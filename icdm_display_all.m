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

function icdm_display_all(group_icdm_file, beta_field_file)
% ICDM_DISPLAY_ALL  Comprehensive visualization of group-level iCDM results.
%
%   icdm_display_all(group_icdm_file, beta_field_file)
%
%   Loads the group iCDM MAT file (mu_ilr_mni, idx_mni, dim_mni, H) and
%   the beta-field MAT file (Beta), reconstructs full 4-D volumes for the
%   ILR mean, compositional probabilities (pi), and regression
%   coefficients, and generates four diagnostic figures:
%     1. ILR mosaic (axial slices for each ILR dimension)
%     2. Pi mosaic  (axial slices for each target parcel)
%     3. Beta fields (mid-axial slice for each predictor)
%     4. 3-panel view (T1 | ILR dim | pi parcel)
%
%   Inputs
%     group_icdm_file : path to group_icdm_iter_*.mat.
%     beta_field_file : path to group_beta_fields.mat or similar.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PLOT_3PANEL, ICDM_PLOT_MOSAIC, ICDM_PLOT_ILR_MOSAIC

%% ============================================================
%% 0. LOAD GROUP FILES
%% ============================================================
G = load(group_icdm_file);
B = load(beta_field_file);

idx = G.idx_mni;
dim = G.dim_mni;
H   = G.H;

K1 = size(G.mu_ilr_mni,2);
K = K1 + 1;
Nm = numel(idx);

fprintf('Loaded group: %s  (Nmni=%d, K=%d)\n', group_icdm_file, Nm, K);
fprintf('Loaded beta fields: %s\n', beta_field_file);

%% ============================================================
%% 1. RECONSTRUCT ILR VOLUMES
%% ============================================================
ilr_vols = zeros([dim K1],'single');
for d = 1:K1
    tmp = zeros(prod(dim),1,'single');
    tmp(idx) = G.mu_ilr_mni(:,d);
    ilr_vols(:,:,:,d) = reshape(tmp, dim);
end

%% ============================================================
%% 2. RECONSTRUCT π VOLUMES
%% ============================================================
% convert ILR mean → π distribution
Z = G.mu_ilr_mni * H.';       % linear map
Z = Z - max(Z,[],2);          % stability
A = exp(Z);
PI = A ./ sum(A,2);

pi_vols = zeros([dim K],'single');
for k = 1:K
    tmp = zeros(prod(dim),1,'single');
    tmp(idx) = PI(:,k);
    pi_vols(:,:,:,k) = reshape(tmp, dim);
end

%% ============================================================
%% 3. RECONSTRUCT β FIELDS
%% ============================================================
P = size(B.Beta,1);     % number of covariates including intercept

beta_vols = zeros([dim P K1],'single');
for p = 1:P
    for d = 1:K1
        tmp = zeros(prod(dim),1,'single');
        tmp(idx) = B.Beta(p,:,d);
        beta_vols(:,:,: ,p,d) = reshape(tmp, dim);
    end
end

%% ============================================================
%% 4. LOAD SPM canonical brain
%% ============================================================
Vt1 = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
t1  = spm_read_vols(Vt1);
mid = round(size(t1,3)/2);

%% ============================================================
%% === DISPLAY: ILR MOSAIC ===
%% ============================================================
figure('Name','ILR Mosaic','Color','w','Position',[50 50 1200 1200]);
ns = 25;
for d = 1:min(K1,9)
    sl = round(linspace(1,dim(3),ns));
    subplot(min(K1,9), ns, (d-1)*ns + 1);
    title(sprintf('ILR dim %d', d));
    for i = 1:ns
        subplot(min(K1,9), ns, (d-1)*ns + i);
        imagesc(rot90(ilr_vols(:,:,sl(i),d)));
        axis off equal;
    end
end
colormap turbo;

%% ============================================================
%% === DISPLAY: π MOSAIC ===
%% ============================================================
figure('Name','Pi mosaic','Color','w','Position',[50 50 1200 1200]);
ns = 20;
for k = 1:min(K,12)
    sl = round(linspace(1,dim(3),ns));
    subplot(min(K,12), ns, (k-1)*ns + 1);
    title(sprintf('Parcel %d', k));
    for i = 1:ns
        subplot(min(K,12), ns, (k-1)*ns + i);
        imagesc(rot90(pi_vols(:,:,sl(i),k))); 
        axis off;
    end
end
colormap turbo;

%% ============================================================
%% === DISPLAY: β FIELDS ===
%% ============================================================
figure('Name','Beta fields','Color','w','Position',[100 100 1400 700]);
for p = 1:P
    subplot(1,P,p);
    imagesc(rot90(beta_vols(:,:,mid,p,1)));   % show slope for ILR dim 1
    axis off;
    title(sprintf('Beta p=%d (ILR dim1)', p));
end
colormap turbo; colorbar;

%% ============================================================
%% === DISPLAY: 3-PANEL (T1 | ILR(d) | π(k)) ===
%% ============================================================
example_ilr = 4;       % ILR dimension to show
example_parcel = 22;   % π parcel to show

figure('Name','3-Panel View','Color','w','Position',[100 100 1600 500]);

subplot(1,3,1);
imagesc(rot90(t1(:,:,mid))); axis off;
title('T1 reference');

subplot(1,3,2);
imagesc(rot90(ilr_vols(:,:,mid,example_ilr)));
axis off; title(sprintf('ILR dim %d', example_ilr));
colormap(gca,'turbo'); colorbar;

subplot(1,3,3);
imagesc(rot90(pi_vols(:,:,mid,example_parcel)));
axis off; title(sprintf('\\pi parcel %d', example_parcel));
colormap(gca,'turbo'); colorbar;

fprintf('\n[Display Complete]\n');

end

