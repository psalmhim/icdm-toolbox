function R = run_minimal_example(workdir, record)
% RUN_MINIMAL_EXAMPLE  Tractogram + atlas -> counts -> population prior -> individualized estimates.
%
%   R = run_minimal_example(workdir)          run and compare with expected_outputs.json
%   R = run_minimal_example(workdir, true)    run and (re)write expected_outputs.json
%
% Steps (input/output contract: docs/CONTRACT.md)
%   0. make_synthetic_cohort : atlas, template, 8 tractograms, identity deformations
%   1. icdm_tracks_to_counts : sub-XX/counts.nii, [X Y Z 7] uint16 streamline--target incidences
%   2. icdm_compose_subject  : native validity mask (sum of counts over all channels > 5) and warps
%   3. icdm_population_eb    : population prior from count-derived representations (targets 1-6),
%                              then subject-level MAP/Laplace inference under the fixed prior
%   4. icdm_posterior_support: N_v, kappa_data, kappa_post, r_data and validity mask in template space
% Requires SPM (SPM12 or SPM25). Runs in about a minute.
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root,'core'), fullfile(root,'atlas'), here);
if nargin < 1 || isempty(workdir), workdir = fullfile(tempdir,'icdm_minimal_example'); end
if nargin < 2, record = false; end
spm('defaults','fmri'); spm_jobman('initcfg');

% 0. synthetic inputs
info = make_synthetic_cohort(fullfile(workdir,'inputs'), 8, 3000, 1);
P = readtable(fullfile(workdir,'inputs','participants.tsv'), 'FileType','text','Delimiter','\t');
names = P.participant_id; age = P.age; S = numel(names);
template = fullfile(workdir,'inputs','template.nii');
atlas = fullfile(workdir,'inputs','atlas.nii');

% 1. counts
counts = cell(S,1);
for s = 1:S
    sd = fullfile(workdir,'inputs',names{s});
    counts{s} = fullfile(sd,'counts.nii');
    icdm_tracks_to_counts(fullfile(sd,'tracts.tck'), atlas, counts{s});
end

% 2-3. prior and individualized estimates (settings of the reported analyses, K = 6 targets)
outdir = fullfile(workdir,'icdm'); if ~exist(fullfile(outdir,'subjs'),'dir'), mkdir(fullfile(outdir,'subjs')); end
targets = info.targets;
opts.mask_thr = 5; opts.mnionly = false; opts.idx_regions = targets;
opts.gamma = 0.35; opts.max_iter = 2; opts.tol = 1e-3; opts.ridge_lambda = 1e-2;
opts.precision = struct('kappa_base',1.0,'prior_eps2',1e-6,'wbeta_prec',false,'alpha_grp',0.25);
opts.vb  = struct('max_iter',20,'pcg_iter',20,'tol_grad',1e-6,'tol_step',1e-6,'jitter',1e-6,'init','hybrid','parfor',false);
opts.agg = struct('trim_alpha',0.20,'min_kappa',1e-4,'min_subj',3,'min_subjects',3,'verbose',0);
opts.w_beta = 0; opts.verbose = 0;
opts.spatial_group = struct('use',true,'lambda',0.5,'neighborhood',6);
opts.paths.temp_dir = fullfile(outdir,'tmp'); if ~exist(opts.paths.temp_dir,'dir'), mkdir(opts.paths.temp_dir); end
for s = 1:S
    si = icdm_compose_subject(counts{s}, {'age'}, age(s), names{s}, template, fullfile(outdir,'subjs'), opts);
    si.idx_regions = targets;
    if s == 1, subjects = repmat(si,S,1); end
    subjects(s) = si;
end
gopts = icdm_evaluate_group_mask(subjects, 0.3, fullfile(outdir,'group_mask.mat'));
opts.idx_mni = gopts.idx_mni; opts.dim_mni = gopts.dim_mni;
subjects = icdm_build_subject_warps(subjects);
grp = icdm_init_group_prior(numel(targets), outdir, opts);
icdm_population_eb(subjects, numel(targets), outdir, opts, grp);

% 4. support indices, written as template-space images
G = load(fullfile(outdir,'group_icdm_iter_002.mat'),'idx_mni','dim_mni','kappa_mni','H','mu_ilr_mni');
od = fullfile(workdir,'outputs'); if ~exist(od,'dir'), mkdir(od); end
Vt = spm_vol(template);
R = struct('id',{},'n_native_valid',{},'n_template_valid',{},'N_total',{},'kappa_data_median',{}, ...
           'kappa_post_median',{},'r_data_median',{});
for s = 1:S
    vbf = fullfile(outdir,'iter_002',[names{s} '_vb.mat']);
    PS = icdm_posterior_support(subjects(s).datafile, vbf, G, 0.25, targets);
    c = load(subjects(s).datafile,'idx_native');
    f = {'kappa_data','kappa_post','r_data','valid'};
    for q = 1:numel(f)
        vol = nan(G.dim_mni(:)'); vol(G.idx_mni) = double(PS.(f{q}));
        Vo = Vt(1); Vo.fname = fullfile(od,sprintf('%s_%s.nii',names{s},f{q})); Vo.dt = [spm_type('float32') 0];
        spm_write_vol(Vo, vol);
    end
    v = PS.valid;
    R(s).id = names{s};
    R(s).n_native_valid = numel(c.idx_native);
    R(s).n_template_valid = nnz(v);
    R(s).N_total = round(sum(PS.N_mni(v)));
    R(s).kappa_data_median = median(PS.kappa_data(v));
    R(s).kappa_post_median = median(PS.kappa_post(v));
    R(s).r_data_median = median(PS.r_data(v));
end
summary.n_group_voxels = numel(G.idx_mni);
summary.subjects = R;
fid = fopen(fullfile(od,'summary.json'),'w'); fprintf(fid,'%s',jsonencode(summary,'PrettyPrint',true)); fclose(fid);

expf = fullfile(here,'expected_outputs.json');
if record
    copyfile(fullfile(od,'summary.json'), expf);
    fprintf('[run_minimal_example] expected outputs recorded in %s\n', expf);
else
    check_outputs(fullfile(od,'summary.json'), expf);
end
end
