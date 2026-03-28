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

function subj=icdm_features_to_mni(subj, opts, featfile)
% ICDM_FEATURES_TO_MNI  Warp native WM mask to MNI and save MNI mask indices.
%
%   subj = icdm_features_to_mni(subj, opts, featfile)
%
%   Warps the subject's native-space white-matter mask (subj.mask_3d) into
%   MNI template space via DARTEL, thresholds the warped mask, and saves
%   the resulting MNI voxel dimension and mask indices for efficient
%   group-level coverage evaluation and index-based computation.
%
%   INPUT
%     subj     : subject struct with .mask_3d, .dartel_flow, .template
%     opts     : options struct (passed to icdm_warp_4d)
%     featfile : output MAT path for dim_mni and idx_mni_mask
%
%   OUTPUT
%     subj : subject struct (with mask_3d populated if absent)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_EVALUATE_GROUP_MASK, ICDM_WARP_4D

if nargin < 3
    error('Usage: icdm_features_to_mni(subj, opts, featfile)');
end

if exist(featfile,'file')
    % Skip if already created
    fprintf('[MNI features] %s exists. Skip.\n', featfile);
    %return;
end

if ~isfield(subj,'mask_3d') || isempty(subj.mask_3d)
    % -------------------------------------------------------------------------
    % WM mask from ICDM sum > thr
    % -------------------------------------------------------------------------
    thr = 5;
    mask_3d=fullfile(fileparts(subj.icdm_4d), sprintf('icdm_sum_mask_thr%d.nii',thr));
    if exist(mask_3d,'file')
        subj.mask_3d = mask_3d;
    else
        V=spm_vol(subj.icdm_4d);
        C = spm_read_vols(V);
        sum_map = sum(C,4);
        mask = sum_map > thr;
        Vm = V(1);
        Vm.fname = mask_3d;
        Vm.dt(1) = spm_type('uint8');
        Vm.pinfo = [1;0;0];
        Vm.n     = [1 1];
        Vm.descrip = 'ICDM WM mask (sum>thr)';
        spm_write_vol(Vm, uint8(mask));
        subj.mask_3d = Vm.fname;
    end
end
if ~isfield(subj,'dartel_flow') || isempty(subj.dartel_flow)
    error('subj.dartel_flow (y_fs_t1) missing');
end
if ~isfield(subj,'template') || isempty(subj.template)
    error('subj.template (DARTEL template) missing');
end

Vm = spm_vol(subj.mask_3d);
M  = spm_read_vols(Vm);
M  = single(M);

% native -> MNI (nearest neighbor)
Yout = icdm_warp_4d(M, Vm, subj.dartel_flow, subj.template, +1, 0, opts);

dim_mni = size(Yout);
if numel(dim_mni) < 3
    dim_mni(3) = 1;
end

mask_mni = Yout > 0.5;
idx_mni_mask = find(mask_mni(:));

save(featfile,'dim_mni','idx_mni_mask','-v7.3');

fprintf('[MNI features] Saved %s (Nmask=%d)\n', featfile, numel(idx_mni_mask));

end

