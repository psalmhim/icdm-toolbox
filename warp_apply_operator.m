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

function Yout = warp_apply_operator(Yin, index8, weight8, dim_out)
% WARP_APPLY_OPERATOR  Apply trilinear warp using precomputed indices.
%
%   Yout = warp_apply_operator(Yin, index8, weight8, dim_out)
%
%   Applies a precomputed 8-neighbour trilinear interpolation warp to
%   flat source data Yin and reshapes the result to dim_out.
%
%   INPUT
%     Yin     : [Nsrc x D] flattened source data
%     index8  : [Nout x 8] uint32 linear indices of 8 neighbours
%     weight8 : [Nout x 8] single trilinear weights
%     dim_out : [1 x 3] output volume dimensions
%
%   OUTPUT
%     Yout : [dim_out x D] warped output volume
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_APPLY_WARP, WARP_APPLY, ICDM_TRILINEAR_INDEX_WEIGHT

D = size(Yin,2);

Yout = zeros(prod(dim_out), D, 'single');

for d = 1:D
    src = Yin(:,d);
    neigh = src(index8);
    val = sum(neigh .* weight8, 2);
    Yout(:,d) = val;
end

Yout = reshape(Yout, [dim_out, D]);
end
