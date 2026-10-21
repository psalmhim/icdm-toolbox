function TT = run_unc4d_example(workdir, config_file)
% RUN_UNC4D_EXAMPLE  Configurable longitudinal iCDM workflow with age-specific (UNC-4D) templates.
%
%   TT = run_unc4d_example(workdir)                 synthetic data (make_synthetic_visits) + unc4d_config.json
%   TT = run_unc4d_example(workdir, config_file)    your own data laid out as in the configuration
%
% Steps (guide: docs/NEW_POPULATION.md)
%   1. template selection   : each visit is assigned the age-specific template whose age bin contains it
%   2. common space         : subject->age-template and common->age-template deformations are composed into
%                             one field on the common grid (sub-XX/ses-YY/dartel/y_fs_t1.nii)
%   3. counts               : tractogram + label image -> counts.nii (icdm_tracks_to_counts)
%   4. reference sample     : visits with reference == 1 build ONE covariate-agnostic population prior
%   5. calibration          : alpha_grp by participant-grouped cross-validation (icdm_select_alpha);
%                             trimming, smoothing, floors and thresholds are held at the documented defaults
%   6. inference            : every visit (reference or not) is estimated under the fixed prior
%   7. organisation         : derivatives/icdm/sub-XX/ses-YY/ outputs + trajectory_table.tsv (one row per visit)
% This example demonstrates configuration and data organisation on synthetic inputs. It is not a
% validation of iCDM in infant data.
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root,'core'), fullfile(root,'atlas'), fullfile(root,'examples','minimal_synthetic'), here);
if nargin < 1 || isempty(workdir), workdir = fullfile(tempdir,'icdm_unc4d_example'); end
if nargin < 2
    config_file = fullfile(here,'unc4d_config.json');
    if ~exist(fullfile(workdir,'data','sessions.tsv'),'file'), make_synthetic_visits(workdir); end
end
cfg = jsondecode(fileread(config_file));
spm('defaults','fmri'); spm_jobman('initcfg');
W = @(p) fullfile(workdir, p);
fill = @(pat, sub, ses) strrep(strrep(pat,'{sub}',erase(sub,'sub-')),'{ses}',erase(ses,'ses-'));

T = readtable(W(cfg.sessions_table), 'FileType','text','Delimiter','\t','TextType','char');
n = height(T); targets = cfg.targets.model_channels(:)'; K = numel(targets);
der = W(cfg.derivatives); common = W(cfg.common_space.template);
tmpl = cfg.templates; if ~iscell(tmpl), tmpl = num2cell(tmpl); end
Tsel = cell(n,1); counts = cell(n,1);

for i = 1:n
    sub = T.participant_id{i}; ses = T.session_id{i};
    % 1. template selection by age bin
    b = find(cellfun(@(t) T.age_months(i) >= t.age_months(1) && T.age_months(i) < t.age_months(2), tmpl), 1);
    if isempty(b), error('No template bin contains age %.1f months (%s %s)', T.age_months(i), sub, ses); end
    Tsel{i} = tmpl{b}.name;
    od = fullfile(der, sub, ses); if ~exist(fullfile(od,'dartel'),'dir'), mkdir(fullfile(od,'dartel')); end
    % 2. compose subject->age template with common->age template: one field on the common grid
    job = struct;
    job.comp{1}.def = {W(fullfile(cfg.bids_root, fill(strrep(cfg.subject_to_template_deformation,'{template}',Tsel{i}), sub, ses)))};
    job.comp{2}.def = {W(tmpl{b}.to_common)};
    job.out{1}.savedef.ofname = 'fs_t1';
    job.out{1}.savedef.savedir.saveusr = {fullfile(od,'dartel')};
    spm_deformations(job);
    % 3. counts
    counts{i} = fullfile(od,'counts.nii');
    icdm_tracks_to_counts(W(fullfile(cfg.bids_root, fill(cfg.tractography.tractogram_per_session, sub, ses))), ...
                          W(fullfile(cfg.bids_root, fill(cfg.targets.label_image_per_session, sub, ses))), counts{i});
end

% options: documented defaults; alpha_grp is set after calibration
fx = cfg.calibration.fixed;
opts.mnionly = false; opts.idx_regions = targets; opts.thresh = fx.native_count_threshold;
opts.gamma = 0.35; opts.max_iter = 2; opts.tol = 1e-3; opts.ridge_lambda = 1e-2;
opts.precision = struct('kappa_base',1.0,'prior_eps2',fx.variance_floor,'wbeta_prec',false, ...
                        'alpha_grp',cfg.calibration.alpha_grp_fixed,'kappa_max',fx.precision_cap);
opts.vb  = struct('max_iter',20,'pcg_iter',20,'tol_grad',1e-6,'tol_step',1e-6,'jitter',1e-6,'init','hybrid','parfor',false);
opts.agg = struct('trim_alpha',fx.trim_alpha,'min_kappa',1e-4,'min_subj',3,'min_subjects',3,'verbose',0);
opts.spatial_group = struct('use',true,'lambda',fx.spatial_smoothing_eta,'neighborhood',6);
opts.w_beta = 0; opts.verbose = 0;
pd = fullfile(der,'prior'); if ~exist(fullfile(pd,'subjs'),'dir'), mkdir(fullfile(pd,'subjs')); end
opts.paths.temp_dir = fullfile(pd,'tmp'); if ~exist(opts.paths.temp_dir,'dir'), mkdir(opts.paths.temp_dir); end

ids = strcat(T.participant_id, '_', T.session_id);
for i = 1:n
    si = icdm_compose_subject(counts{i}, {'age_months'}, T.age_months(i), ids{i}, common, fullfile(pd,'subjs'), opts);
    si.idx_regions = targets;
    if i == 1, subjects = repmat(si,n,1); end
    subjects(i) = si;
end
ref = find(T.(cfg.reference.column) == 1);
if ~strcmp(cfg.reference.prior_mode,'pooled')
    error('This example implements prior_mode = "pooled" (age-agnostic prior for trajectory analysis).');
end

% 4. population prior from the reference visits (estimated once, covariate-agnostic)
gopts = icdm_evaluate_group_mask(subjects(ref), fx.group_coverage, fullfile(pd,'group_mask.mat'));
opts.idx_mni = gopts.idx_mni; opts.dim_mni = gopts.dim_mni;
subjects = icdm_build_subject_warps(subjects);
grp0 = icdm_init_group_prior(K, pd, opts);
icdm_population_eb(subjects(ref), K, pd, opts, grp0);
prior_file = fullfile(pd,'group_icdm_iter_002.mat');

% 5. alpha_grp: participant-grouped cross-validated predictive likelihood (no covariate used)
if ischar(cfg.calibration.alpha_grp) && strcmp(cfg.calibration.alpha_grp,'cv')
    A = icdm_select_alpha(fullfile(pd,'iter_002',strcat(ids(ref),'_vb.mat')), prior_file, ...
        struct('nfold',cfg.calibration.cv_folds,'groups',{T.(cfg.calibration.cv_group_by)(ref)}));
    alpha = A.alpha_cv;
else
    alpha = cfg.calibration.alpha_grp_fixed;
end
fid = fopen(fullfile(der,'calibration.json'),'w');
fprintf(fid,'%s',jsonencode(struct('alpha_grp',alpha,'criterion',cfg.calibration.alpha_grp, ...
    'reference_visits',numel(ref),'fixed',fx),'PrettyPrint',true)); fclose(fid);

% 6. every visit under the fixed prior with the calibrated alpha
opts.precision.alpha_grp = alpha;
grp = load(prior_file);
fd = fullfile(der,'final'); if ~exist(fd,'dir'), mkdir(fd); end
icdm_run_subject_vb(subjects, K, grp, opts, fd);

% 7. per-visit outputs and the trajectory table
Vt = spm_vol(common);
rows = cell(n,1);
for i = 1:n
    od = fullfile(der, T.participant_id{i}, T.session_id{i});
    vbf = fullfile(fd,[ids{i} '_vb.mat']); copyfile(vbf, fullfile(od,'icdm_vb.mat'));
    PS = icdm_posterior_support(subjects(i).datafile, vbf, grp, alpha, targets);
    f = {'kappa_data','kappa_post','r_data','valid'};
    for q = 1:numel(f)
        vol = nan(grp.dim_mni(:)'); vol(grp.idx_mni) = double(PS.(f{q}));
        Vo = Vt(1); Vo.fname = fullfile(od,[f{q} '.nii']); Vo.dt = [spm_type('float32') 0];
        spm_write_vol(Vo, vol);
    end
    v = PS.valid;
    rows{i} = {T.participant_id{i}, T.session_id{i}, T.age_months(i), Tsel{i}, T.(cfg.reference.column)(i), ...
               alpha, nnz(v), median(PS.kappa_data(v)), median(PS.kappa_post(v)), median(PS.r_data(v)), ...
               fullfile(T.participant_id{i}, T.session_id{i}, 'icdm_vb.mat')};
end
TT = cell2table(vertcat(rows{:}), 'VariableNames', {'participant_id','session_id','age_months','template', ...
    'reference','alpha_grp','n_valid_voxels','kappa_data_median','kappa_post_median','r_data_median','estimate_file'});
writetable(TT, fullfile(der,'trajectory_table.tsv'), 'FileType','text','Delimiter','\t');
fprintf('[run_unc4d_example] %d visits; alpha_grp = %.3g; table: %s\n', n, alpha, fullfile(der,'trajectory_table.tsv'));
end
