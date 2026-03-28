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

function Lout = warp_label_trilinear(Lin, dim_in, idx_in, subj, grp, K)
% WARP_LABEL_TRILINEAR  Warp a discrete label map via one-hot trilinear argmax.
%
%   Lout = warp_label_trilinear(Lin, dim_in, idx_in, subj, grp, K)
%
%   Converts a native-space label vector to one-hot channels, warps each
%   channel to MNI space with trilinear interpolation, and reconstructs
%   the warped label map by argmax voting.
%
%   INPUT
%     Lin    : [N x 1] label vector (0-based class indices)
%     dim_in : [1 x 3] native volume dimensions
%     idx_in : [N x 1] linear indices of labelled voxels
%     subj   : subject struct with warp operator
%     grp    : group struct with .dim_mni
%     K      : number of classes
%
%   OUTPUT
%     Lout : [prod(dim_mni) x 1] warped label map in MNI space
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_WARP_TO_MNI, ICDM_APPLY_WARP

% Convert label to one-hot channels
Lin_full = zeros(prod(dim_in),1,'uint16');
Lin_full(idx_in) = Lin;
Lin3d = reshape(Lin_full, dim_in);

P = zeros([dim_in K],'single');
for k = 1:K
    P(:,:,:,k) = single(Lin3d == (k-1));
end

% Warp each channel using trilinear
P_mni = zeros([grp.dim_mni K],'single');
for k = 1:K
    Ptmp = warp_to_mni(P(:,:,:,k), subj, grp);  % your operator
    P_mni(:,:,:,k) = Ptmp;
end

% Argmax
[~, Lout3] = max(P_mni, [], 4);
Lout = Lout3(:);

end
