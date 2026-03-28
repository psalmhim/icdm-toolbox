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

function Yout = icdm_apply_warp(Yin, warp)
% ICDM_APPLY_WARP  Apply a prebuilt trilinear warp operator to a 4-D array.
%
%   Yout = icdm_apply_warp(Yin, warp)
%
%   Performs fast spatial warping of a 4-D volume using the 8-neighbor
%   trilinear interpolation indices and weights stored in the warp
%   structure.  Each 3-D frame is warped independently.
%
%   This function is the core engine behind icdm_warp_to_mni and
%   icdm_warp_to_native, replacing SPM-based spm_deformations calls
%   for speed during iterative subject-VB and group-EB loops.
%
%   Inputs
%     Yin  : [d1 x d2 x d3 x D] single array in the source space.
%            Spatial dimensions must match warp.dim_in.
%     warp : struct with fields:
%              .index   [N_out x 8] uint32 neighbor indices into source.
%              .weight  [N_out x 8] single trilinear weights.
%              .dim_in  [1 x 3] source volume dimensions.
%              .dim_out [1 x 3] target volume dimensions.
%
%   Output
%     Yout : [dim_out(1) x dim_out(2) x dim_out(3) x D] warped single array.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_BUILD_WARP_TO_MNI, ICDM_BUILD_WARP_TO_NATIVE,
%            ICDM_TRILINEAR_INDEX_WEIGHT

dim_in  = warp.dim_in;
dim_out = warp.dim_out;

if ~isequal(size(Yin,1:3), dim_in)
    error('icdm_apply_warp: input dim mismatch.');
end

D = size(Yin,4);

Yin_flat = reshape(Yin,[],D);  % [N_in × D]
index8  = warp.index;
weight8 = warp.weight;

Yout = zeros([dim_out D],'single');
for d = 1:D
    src = Yin_flat(:,d);
    neigh = src(index8);        % [Nout × 8]
    dst = sum(neigh .* weight8,2);
    data=reshape(dst, dim_out);
    Yout(:,:,:,d) = data;
end
end
