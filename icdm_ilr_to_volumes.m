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

function ilr_vols = icdm_ilr_to_volumes(mu_ilr, idx_mni, dim_mni)
% ICDM_ILR_TO_VOLUMES  Expand index-form ILR data into 4-D volume arrays.
%
%   ilr_vols = icdm_ilr_to_volumes(mu_ilr, idx_mni, dim_mni)
%
%   Converts compact index-based ILR coordinates [Nmni x (K-1)] into a
%   full 4-D volume array [X x Y x Z x (K-1)] by placing values at the
%   linear indices specified by idx_mni.
%
%   INPUT
%     mu_ilr  : [Nmni x (K-1)] ILR coordinates
%     idx_mni : [Nmni x 1] linear indices into the MNI volume
%     dim_mni : [1 x 3] MNI volume dimensions [X Y Z]
%
%   OUTPUT
%     ilr_vols : [X x Y x Z x (K-1)] 4-D volume of ILR coordinates
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_2D_TO_4D, ICDM_DISPLAY_ALL

K1 = size(mu_ilr,2);
V = prod(dim_mni);

ilr_vols = zeros([dim_mni K1],'single');

for d = 1:K1
    tmp = zeros(V,1,'single');
    tmp(idx_mni) = mu_ilr(:,d);
    ilr_vols(:,:,:,d) = reshape(tmp, dim_mni);
end
end
