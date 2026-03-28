% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
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

function icdm_write_group_nifti(group_icdm_mat, group_beta_mat, outdir)
% ICDM_WRITE_GROUP_NIFTI  Write group iCDM results as SPM-compatible NIfTI.
%
%   icdm_write_group_nifti(group_icdm_mat, group_beta_mat, outdir)
%
%   Converts the compact index-based group iCDM output (ILR mean mu,
%   reliability kappa) and voxelwise beta coefficients into full 3-D and
%   4-D NIfTI volumes for viewing in SPM or any NIfTI viewer.
%
%   Written files:
%     group_mu_ilr.nii          : 4-D ILR mean [X Y Z K1]
%     group_kappa.nii           : 3-D reliability map
%     group_beta_all.nii        : 4-D beta (all P*K1 frames)
%     group_beta_p*_ilr*.nii    : individual 3-D beta volumes
%
%   INPUT
%     group_icdm_mat : MAT file with mu_ilr_mni, kappa_mni, idx_mni, dim_mni
%     group_beta_mat : MAT file with Beta [P x Nmni x K1]
%     outdir         : output directory for NIfTI files
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_RESULT_TO_NII, ICDM_POPULATION_EB

if nargin < 3
    error('Usage: icdm_write_group_nifti(group_icdm_mat, group_beta_mat, outdir)');
end

if ~exist(outdir,'dir'), mkdir(outdir); end

%% -------------------------------------------------------------
%% Load both files
%% -------------------------------------------------------------
G = load(group_icdm_mat);
B = load(group_beta_mat);

mu     = G.mu_ilr_mni;       % [Nmni × K1]
kappa  = G.kappa_mni;        % [Nmni × 1]
idx    = G.idx_mni(:);
dim    = double(G.dim_mni(:)');
K1     = size(mu,2);

Beta   = B.Beta;             % [P × Nmni × K1]
P      = size(Beta,1);

fprintf('[WRITE] dim_mni = [%d %d %d]\n', dim);

%% Helper: empty volume
make_vec = @(dim) zeros(prod(dim),1,'single');

%% -------------------------------------------------------------
%% 1) Write μ (ILR mean) as 4D NIfTI
%% -------------------------------------------------------------
mu_4d = zeros(dim(1),dim(2),dim(3),K1,'single');
vec = make_vec(dim);

for d = 1:K1
    vec(:) = 0;
    vec(idx) = mu(:,d);
    mu_4d(:,:,:,d) = reshape(vec,dim);
end

fname_mu = fullfile(outdir,'group_mu_ilr.nii');
write_4d_spm(mu_4d, fname_mu);
fprintf('[WRITE] %s\n', fname_mu);

%% -------------------------------------------------------------
%% 2) Write κ
%% -------------------------------------------------------------
vec(:) = 0;
vec(idx) = kappa;
kap3d = reshape(vec,dim);

fname_k = fullfile(outdir,'group_kappa.nii');
write_3d_spm(kap3d, fname_k);
fprintf('[WRITE] %s\n', fname_k);

%% -------------------------------------------------------------
%% 3) Write β
%% -------------------------------------------------------------
beta_all = zeros(dim(1),dim(2),dim(3),P*K1,'single');
cnt = 1;

for p = 1:P
    for d = 1:K1
        vec(:) = 0;
        vec(idx) = Beta(p,:,d);
        vol = reshape(vec,dim);

        beta_all(:,:,:,cnt) = vol;

        % individual file
        fname_pd = fullfile(outdir, sprintf('group_beta_p%d_ilr%d.nii',p,d));
        write_3d_spm(vol, fname_pd);
        cnt = cnt + 1;
    end
end

fname_all = fullfile(outdir,'group_beta_all.nii');
write_4d_spm(beta_all, fname_all);
fprintf('[WRITE] %s\n', fname_all);

end


%% =======================================================================
%% SPM WRITERS
%% =======================================================================

function write_3d_spm(vol3d, fname)
V = struct();
V.fname = fname;
V.dim   = size(vol3d);           % [X Y Z]
V.dt    = [spm_type('float32') 0];
V.mat   = eye(4);
spm_write_vol(V, vol3d);
end

function write_4d_spm(vol4d, fname)
sz = size(vol4d);    % X Y Z T

for t = 1:sz(4)
    V = struct();
    V.fname = fname;
    V.dim   = sz(1:3);          % MUST be 3 numbers only
    V.dt    = [spm_type('float32') 0];
    V.mat   = eye(4);
    V.n     = [t 1];            % 4D indexing

    spm_write_vol(V, vol4d(:,:,:,t));
end
end

