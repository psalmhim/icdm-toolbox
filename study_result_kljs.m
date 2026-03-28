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

% STUDY_RESULT_KLJS  Example script for split-half JS/KL divergence computation.
%
%   Demonstrates the use of icdm_cdm_kl_js to compute voxelwise
%   KL and Jensen-Shannon divergence between two split-half ICDM
%   NIfTI files (raw streamline counts) with an aligned brain mask.
%
%   See also ICDM_CDM_KL_JS, ICDM_KL_JS
% Author: Hae-Jeong Park, Ph.D.

%% Example: split-half JS/KL computation

% input NIfTI files (raw counts)
nii_fold1 = ...
'/Volumes/eHDD/DATA/epilepsy/NDARDK983BDA/...
WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold1.nii';

nii_fold2 = ...
'/Volumes/eHDD/DATA/epilepsy/NDARDK983BDA/...
WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold2.nii';

% aligned brain mask
mask_file = ...
'/Volumes/eHDD/DATA/epilepsy/NDARDK983BDA/rfs_t1_brain.nii';

% run KL / JS computation
[outKL, outJS] = icdm_cdm_kl_js( ...
    nii_fold1, ...
    nii_fold2, ...
    mask_file, ...
    'split_half');

fprintf('Saved:\n%s\n%s\n', outKL, outJS);