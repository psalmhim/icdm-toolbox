function run_icdm_example(icdm_files, names, covariates, template, outdir)
% RUN_ICDM_EXAMPLE  Hierarchical iCDM inference for a cohort (settings used in the paper).
%
%   icdm_files : cell array of per-subject 4-D count images (voxel x target streamline--target
%                incidence counts in native diffusion space)
%   names      : cell array of subject identifiers
%   covariates : struct with fields .names (e.g. {'age','sex'}) and .values ([S x P])
%   template   : path to the DARTEL/MNI template used for spatial normalization
%   outdir     : output folder
%
% Requires SPM on the MATLAB path and the core/ folder of this toolbox.
% The population prior is estimated once from prior-independent count-derived
% representations and held fixed during subject-level inference.

if ~exist(outdir,'dir'), mkdir(outdir); end
S = numel(icdm_files);
idx_regions = 1:68;                                   % Desikan--Killiany cortical targets

% ---- options (defaults used for the reported analyses) ----------------------------
opts.mask_thr    = 5;           % native voxel inclusion: >5 streamline--target incidences
opts.label       = 'desikan68';
opts.mnionly     = false;
opts.idx_regions = idx_regions;
opts.gamma       = 0.35;        % power transform used ONLY for numerical initialization
opts.max_iter    = 2;           % single forward pass (the prior is data-derived)
opts.tol         = 1e-3;
opts.ridge_lambda= 1e-2;
opts.precision   = struct('kappa_base',1.0,'prior_eps2',1e-6,'wbeta_prec',false, ...
                          'alpha_grp',0.25);          % population-transfer strength
opts.vb  = struct('max_iter',20,'pcg_iter',20,'tol_grad',1e-6,'tol_step',1e-6, ...
                  'jitter',1e-6,'init','hybrid','parfor',true);   % Newton--CG settings
opts.agg = struct('trim_alpha',0.20,'min_kappa',1e-4,'min_subj',3,'min_subjects',3,'verbose',1);
opts.w_beta = 0;                                   % covariate-agnostic prior (association-testing mode)
opts.spatial_group = struct('use',true,'lambda',0.5,'neighborhood',6);   % group-prior smoothing (eta)
opts.paths.temp_dir = fullfile(outdir,'tmp'); if ~exist(opts.paths.temp_dir,'dir'), mkdir(opts.paths.temp_dir); end
opts.grp_kappa_thresh = 0;
opts.estimator = 'ols';        % downstream association by ordinary least squares
opts.verbose = 0;

% ---- 1. count-derived subject representations (native -> template) -----------------
V0 = spm_vol(icdm_files{1}); V0 = V0(idx_regions); K = numel(V0);
for i = 1:S
    si = icdm_compose_subject(icdm_files{i}, covariates.names, covariates.values(i,:), ...
                              names{i}, template, fullfile(outdir,'subjs'), opts);
    si.idx_regions = idx_regions;
    if i == 1, subjects = repmat(si,S,1); end
    subjects(i) = si;
end

% ---- 2. template-space analysis domain and warps ------------------------------------
gopts = icdm_evaluate_group_mask(subjects, 0.3, fullfile(outdir,'group_mask.mat'));
opts.idx_mni = gopts.idx_mni; opts.dim_mni = gopts.dim_mni;
subjects = icdm_build_subject_warps(subjects);
opts.design_spec = [{'const'} covariates.names(1)];
opts = icdm_prepare_design_vector(subjects, opts);

% ---- 3. population prior (estimated once) and subject-level Newton--Laplace inference -
grp = icdm_init_group_prior(K, outdir, opts);
icdm_population_eb(subjects, K, outdir, opts, grp);
end
