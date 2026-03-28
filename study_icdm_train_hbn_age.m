% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI", NeuroImage (2026).
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the copyright holder nor the names of its
%      contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES.
%
% REQUIRES: SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

% STUDY_ICDM_TRAIN_HBN_AGE  Train iCDM population model on HBN cohort (age 8-13, DK68).
%
%   End-to-end training script for the iCDM empirical Bayes pipeline on
%   the Healthy Brain Network (HBN) dataset restricted to the age 8-13
%   range.  Uses the Desikan-Killiany 68-parcel atlas with a linear age
%   covariate (design_spec = {'const','age'}).  Steps:
%     1. Scan and filter subjects by demographics
%     2. Compose subject structs and build DARTEL warp operators
%     3. Evaluate group MNI mask (30 % coverage threshold)
%     4. Initialise group prior and run population EB
%     5. Post-hoc prediction and age-effect mapping
%
%   See also ICDM_POPULATION_EB, ICDM_COMPOSE_SUBJECT,
%            ICDM_PREDICT_INDIVIDUAL_CDM, ICDM_BETA_AGE_EFFECT_MAP
% Author: Hae-Jeong Park, Ph.D.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fully aligned with Method section rules:
%
%   1) No FA, ReHo, Jacobian, lesion
%   2) No spatial smoothing (any kind)
%   3) No ILR neighbor smoothing
%   4) Group prior = Trimmed MAD + kappa_base only
%   5) gamma-tempering ONLY: C_gamma = (C+1)^gamma
%   6) Full DARTEL forward/backward warping
%   7) Index-based storage to reduce file size
%   8) Skips existing files automatically
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ===== USER CONFIGURATION (edit before running) =====
% mainPath        = '/path/to/your/hbn_connectome/data';
% excel_file      = fullfile(mainPath, 'codes','HBN_connectome_basic_demos_v2.xlsx');
% dartel_template = '/path/to/your/DARTEL_template/HBN_Temp_6.nii';
% icdm_rel        = 'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii';
% atlaslabel      = 'desikan68';
% outdir          = '/path/to/your/output/ICDM_desikan68';
error('Edit the USER CONFIGURATION section above before running this script.');

%% 0) Paths
mainPath   = '/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/';
excel_file = fullfile(mainPath, 'codes','HBN_connectome_basic_demos_v2.xlsx');
dartel_template = '/remotenas2/projects/hbn_connectome/HBN_DARTEL_TEMP/HBN_Temp_6.nii';

% DK+ ICDM
icdm_rel   = 'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii';
%icdm_rel   = 'WBT_10M_ctx_vep165_to_dwi15_icdm.nii';

atlaslabel='desikan68';
outdir = fullfile('/home/hjpark/data/', ['ICDM_' atlaslabel]);
if ~exist(outdir,'dir'), mkdir(outdir); end

outdirsubj = fullfile(outdir,'subjs');
if ~exist(outdirsubj,'dir'), mkdir(outdirsubj); end

tmpSSD = fullfile(outdir, 'tmp');
if ~exist(tmpSSD,'dir'), mkdir(tmpSSD); end
opts.paths.temp_dir = tmpSSD;

%% 

opts.group_mask='vtemplate.mat';
opts.mask_thr=5;
opts.label = atlaslabel;

idx_regions= 1:68;
%idx_regions=1:162;
a=load(opts.group_mask);
if ~isfield(a,'idx_brain_mni')
    fn=fullfile(mnet_default_path('labels'),'mask_gmwm.nii');
    vtpl=spm_vol(fn); A=spm_read_vols(vtpl);
    idx_brain_mni = uint32(find(A>0.1));
    save(opts.group_mask,'idx_brain_mni','-append');
end

subject_id_column = 'EID'; 
age_column = 'Age'; 
sex_column = 'Sex';

%% 1) Dependencies
if exist('spm','file') ~= 2
    error('SPM12 not on path');
end

%% 2) Scan subjects (UNCHANGED)
lst = dir(fullfile(mainPath, 'NDA*'));
names={}; icdm_all={};

for l=1:numel(lst)
    if lst(l).isdir && numel(lst(l).name) > 3
        f = fullfile(mainPath, lst(l).name, icdm_rel);
        dartel = fullfile(mainPath,lst(l).name,'dartel','y_fs_t1.nii');
        if exist(f,'file') && exist(dartel,'file')
            names{end+1}   = lst(l).name;
            icdm_all{end+1}= f;
        end
    end
end

names=names(:); icdm_all=icdm_all(:);
if isempty(icdm_all), error('No ICDM found'); end

%% 3) Demographics
T = readtable(excel_file);
[lia, locb] = ismember(names, T.(subject_id_column));
names=names(lia);
icdm_all=icdm_all(lia);
locb=locb(lia);

age = double(T.(age_column)(locb));
sex_raw = T.(sex_column)(locb);

sex = nan(size(sex_raw));
if iscell(sex_raw)
    for i=1:numel(sex_raw)
        s = string(sex_raw{i});
        if any(lower(s)==["m","male","1"]), sex(i)=1; continue; end
        if any(lower(s)==["f","female","0","2"]), sex(i)=0; continue; end
        d=str2double(s); if ~isnan(d), sex(i)=d; end
    end
else
    sex=double(sex_raw);
end
sex(sex==2)=0;

valid = ~isnan(age) & ~isnan(sex) & age>8 & age<13;
names=names(valid);
icdm_all=icdm_all(valid);
age=age(valid);
sex=sex(valid);
Nsubj = numel(names);

fprintf('[INFO] Final N = %d subjects\n', Nsubj);

%% 4. Infer K
V0 = spm_vol(icdm_all{1}); V0=V0(idx_regions);
K = numel(V0);
fprintf('[INFO] K = %d\n', K);

opts.mnionly = false;
opts.idx_regions = idx_regions;
%% 7. OPTIONS (final)

opts.gamma = 0.35;

opts.max_iter = 3;
opts.tol = 1e-3;
opts.stop_on_elbo = true;
opts.ridge_lambda = 1e-2;

opts.precision = struct( ...
    'kappa_base', 1.0, ...
    'prior_eps2', 1e-6 ...
);

opts.vb = struct( ...
    'max_iter', 20, ...
    'pcg_iter', 20, ...
    'tol_grad',1e-6, ...
    'tol_step',1e-6, ...
    'jitter',1e-6, ...
    'init','hybrid', ...
    'parfor',true, ...
    'smooth_ilr',false, ...
    'smooth_conn',0, ...
    'smooth_tau',0 ...
);

opts.agg=struct( ...
    'trim_alpha', 0.20, ...    % trimming proportion (two-sided)
    'min_kappa', 1e-4, ...     % minimum kappa for group prior
    'min_subj', 3, ...         % minimum subjects per voxel
    'verbose', 1 ...           % verbosity level
);
opts.w_beta = 0.35;
opts.verbose = 1;
opts.paths.temp_dir = tmpSSD;

opts.spatial = struct( ...
    'use', false, ...        % if true → use spatial prior
    'lambda', 0.5, ...          % weight of spatial penalty
    'neighborhood', 6, ...            % neighbor radius (6/18/26 or voxel-radius)
    'iterations', 1 ...         % number of smoothing passes
);

opts.spatial_group.use         = true;
opts.spatial_group.lambda      = 0.5;   % or whatever you used before
opts.spatial_group.neighborhood = 6;    % or 26, depending on your choice


opts.debug = struct( ...
    'use', true, ...               % master switch for debugging
    'save_figures', true, ...          % save PNG to disk
    'figure_dir', './figures', ...              % will be set in main script
    'show_coronal_slice', true, ...    % MU_nat_4D, KAP_nat_3D before/after VB
    'show_beta_maps', true, ...        % m_beta and fusion effect
    'slice_z', [] ...                  % if empty → middle slice
);

%% 5. Compose subject structs
subjects(Nsubj,1) = struct();
for i=1:Nsubj
    fprintf('[INFO] Composing subject %d/%d: %s\n', i, Nsubj, names{i});
    si = icdm_compose_subject(icdm_all{i},{'age','sex'} ,[age(i), sex(i)], names{i}, dartel_template,outdirsubj,opts);
    si.idx_regions=opts.idx_regions;
    if i==1, subjects=repmat(si,Nsubj,1); end
    subjects(i) = si;
end

%% 8. Group mask
outfile=fullfile(outdir,'group_mask.mat');
gopts = icdm_evaluate_group_mask(subjects, 0.3, outfile);
opts.idx_mni      = gopts.idx_mni;
opts.dim_mni      = gopts.dim_mni;


%% 6. Build warp operators
subjects = icdm_build_subject_warps(subjects);

%% 7. Prepare design vectors
opts.design_spec = {'const','age'};
%opts.design_spec = { 'const', struct('var','age','type','spline','K',3) };
opts = icdm_prepare_design_vector(subjects,opts);


%% 9. INIT group prior (index-based)
grp = icdm_init_group_prior(K, outdir, opts);  % Not needed for standard EB

%% 10. RUN EB
tic
icdm_population_eb(subjects, K, outdir, opts, grp);
toc
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(outdir);
load group_mask.mat
icdm_result_to_nii(outdir);
nid=90;
subj=subjects(nid);
sid  = subj.id;
subjfile = fullfile(outdir,'iter_003',[sid '_vb.mat']);
VBs = load(subjfile);
grp=load(fullfile(outdir,'group_icdm_iter_003.mat'));
VBc = icdm_predict_individual_cdm(subj, VBs, grp, opts);

load vtemplate.mat
R = icdm_2d_to_4d(VBc.pred.pi_grp, idx_mni, [dim_mni K]);
V = vtpl;
V.fname = fullfile(outdir,'s90_pi_grp003.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);

R = icdm_2d_to_4d(VBc.pred.pi_ind, idx_mni, [dim_mni K]);
V = vtpl;
V.fname = fullfile(outdir,'s90_pi_ind.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);

R = icdm_2d_to_4d(VBc.pred.pi_fused, idx_mni, [dim_mni K]);
V = vtpl;
V.fname = fullfile(outdir,'s90_pi_fused.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);

nid=10;
subj=subjects(nid);
sid  = subj.id;
subjfile = fullfile(outdir,'iter_003',[sid '_vb.mat']);
VBs = load(subjfile);
grp=load(fullfile(outdir,'group_icdm_iter_003.mat'));
VBc = icdm_predict_individual_cdm(subj, VBs, grp, opts);

load vtemplate.mat
R = icdm_2d_to_4d(VBc.pred.pi_grp, idx_mni, [dim_mni K]);
V = vtpl;
V.fname = fullfile(outdir,'s10_pi_grp003.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);

R = icdm_2d_to_4d(VBc.pred.pi_ind, idx_mni, [dim_mni K]);
V = vtpl;
V.fname = fullfile(outdir,'s10_pi_ind.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);

R = icdm_2d_to_4d(VBc.pred.pi_fused, idx_mni, [dim_mni K]);
V = vtpl;
V.fname = fullfile(outdir,'s10_pi_fused.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);



R = icdm_predict_pi_for_age(12, grp.Beta, grp.H, opts);
V = vtpl;
V.fname = fullfile(outdir,'P_age12.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);

R = icdm_predict_pi_for_age(8, grp.Beta, grp.H, opts);
V = vtpl;
V.fname = fullfile(outdir,'P_age8.nii');
V.dim = size(R);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, R);


% age_index=2:4;
age_index=2;
E = icdm_beta_age_effect_map(grp.Beta, age_index);
E=icdm_2d_to_4d(E, idx_mni, dim_mni);
V = vtpl;
V.fname = fullfile(outdir,'age_effect.nii');
V.dim = size(E);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, E);


reg_image(fullfile(outdir,['P_age8.nii' '',3']), fullfile(outdir,['P_age12.nii' '',3']));