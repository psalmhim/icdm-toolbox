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

function icdm_result_to_nii(matfile)
% ICDM_RESULT_TO_NII  Convert a group iCDM iteration MAT file to NIfTI volumes.
%
%   icdm_result_to_nii(matfile)
%
%   Loads the group-level iteration result (mu_ilr_mni, Beta, tau2,
%   kappa_mni, idx_mni, dim_mni, H) from matfile and writes the
%   following NIfTI volumes to the same directory:
%     - group_kappa.nii      (reliability map)
%     - group_mu_ilr.nii     (ILR mean, 4-D)
%     - group_pi.nii         (compositional probabilities, 4-D)
%     - group_tau2_ilr.nii   (ILR variance, 4-D)
%     - group_beta1/2.nii    (beta fields per predictor)
%
%   If matfile is a directory, the default filename
%   group_icdm_iter_003.mat is appended.
%
%   Input
%     matfile : path to the group iteration MAT file, or its parent dir.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_WRITE_GROUP_NIFTI, ICDM_CHECK_GROUP_ITER
if isfolder(matfile)
    matfile=fullfile(matfile,'group_icdm_iter_003.mat');
end

load vtemplate.mat;
[outdir,f,e]=fileparts(matfile);
id=find(f=='_');
itername=f(id(end)+1:end);
%% ---------------------------------------------------------
% Load group ICDM results
%% ---------------------------------------------------------
fprintf('Loading %s..\n', matfile);
if ~exist(matfile, 'file')
    error('Matfile %s does not exist.', matfile);
end

load(matfile);
required_fields = {'mu_ilr_mni', 'Beta', 'tau2', 'kappa_mni', 'idx_mni', 'dim_mni', 'H'};
for i = 1:numel(required_fields)
    field = required_fields{i};
    if ~exist(field, 'var')
        error('Missing required variable ''%s'' in %s.', field, source_name);
    end
end

K1 = size(mu_ilr_mni,2);K  = K1 + 1;

%% ---------------------------------------------------------
% Reconstruct 4D volumes
%% ---------------------------------------------------------
fprintf('Reconstructing full 4D volume..\n');
% ILR mean
MU_ilr_mni = icdm_2d_to_4d(mu_ilr_mni, idx_mni, [dim_mni K1]);

% ILR variance
TAU_mni = icdm_2d_to_4d(tau2, idx_mni, [dim_mni K1]);

% Reliability (kappa)
KAPPA_mni = icdm_2d_to_4d(kappa_mni, idx_mni, dim_mni);

%% ---------------------------------------------------------
% Convert ILR mean back to softmax connectivity
%% ---------------------------------------------------------
fprintf('Back-transforming ILR to pi...\n');
PI_mni = icdm_2d_to_4d(ilr_inverse(mu_ilr_mni,H), idx_mni, [dim_mni K]);

Beta1_mni = icdm_2d_to_4d(squeeze(Beta(1,:,:)), idx_mni, [dim_mni size(Beta,3)]);
Beta2_mni = icdm_2d_to_4d(squeeze(Beta(2,:,:)), idx_mni, [dim_mni size(Beta,3)]);
%% ---------------------------------------------------------
% (Optional) Save as NIfTI
%% ---------------------------------------------------------

fprintf('Saving NIfTI file..\n');
% Save kappa
V = vtpl;
V.fname = fullfile(outdir,[itername '_group_kappa.nii']);
V.dim = size(KAPPA_mni);
V.dt = [spm_type('float32') 0];
vol_write_vols(V, KAPPA_mni);

V.fname=fullfile(outdir,[itername '_group_mu_ilr.nii']);
V.dim = size(MU_ilr_mni);
vol_write_vols(V,MU_ilr_mni);

V.fname=fullfile(outdir,[itername '_group_pi.nii']);
V.dim = size(PI_mni);
vol_write_vols(V,PI_mni);

V.fname=fullfile(outdir,[itername '_group_tau2_ilr.nii']);
V.dim = size(TAU_mni);
vol_write_vols(V,TAU_mni);

V.fname=fullfile(outdir,[itername '_group_beta1.nii']);
V.dim = size(Beta1_mni);
vol_write_vols(V,Beta1_mni);

V.fname=fullfile(outdir,[itername '_group_beta2.nii']);
V.dim = size(Beta2_mni);
vol_write_vols(V,Beta2_mni);
end
