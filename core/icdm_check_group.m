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

function icdm_check_group(group_path)
% ICDM_CHECK_GROUP  Materialise group beta coefficients as NIfTI volumes.
%
%   icdm_check_group(group_path)
%
%   Loads the compact index-based beta file (beta_ilr.mat) from group_path
%   and writes per-covariate 4-D NIfTI volumes of ILR beta coefficients
%   for visual inspection.
%
%   INPUT
%     group_path : directory containing beta_ilr.mat
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_CHECK_GROUP_ITER, ICDM_ESTIMATE_BETA

% materialize_beta_nii(beta_idx_mat, beta_nii)
beta_idx_mat = fullfile(group_path, 'beta_ilr.mat');

load(beta_idx_mat);
materialize_beta_nii_from_idx(Beta, idx_mni, dim_mni, Vtpl, beta_nii);

end


function materialize_beta_nii_from_idx(Beta, idx_mni, dim_mni, Vtpl, out_path)
[P, ~, Km1] = size(Beta);
F = P*Km1;
Vout = repmat(Vtpl, F, 1);
[p1,f1,e1]=fileparts(out_path);
for p=1:P
    out_path1 = fullfile(p1,sprintf('%s_beta%d.nii',f1,p));
    for f=1:Km1
        Vout(f).fname   = out_path1;
        Vout(f).n       = [f 1];
        Vout(f).dt      = [spm_type('float32') 0];
        Vout(f).pinfo   = [1;0;0];
        Vout(f).descrip = sprintf('β (ILR) p=%d, d=%d', p, d);
        vol = zeros(prod(dim_mni),1,'single');
        vol(idx_mni) = single(Beta(p,:,d));
        spm_write_vol(Vout(f), reshape(vol,dim_mni));
    end
end
end
