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

function [index8, weight8] = icdm_trilinear_index_weight(x,y,z,dim)
% ICDM_TRILINEAR_INDEX_WEIGHT  Compute 8-neighbour trilinear indices and weights.
%
%   [index8, weight8] = icdm_trilinear_index_weight(x, y, z, dim)
%
%   For each query point (x,y,z) computes the eight corner linear indices
%   and corresponding trilinear interpolation weights within a volume of
%   size dim.  Out-of-bound coordinates are clipped to the volume edges.
%
%   INPUT
%     x, y, z : [N x 1] continuous voxel coordinates (1-based)
%     dim     : [1 x 3] volume dimensions [X Y Z]
%
%   OUTPUT
%     index8  : [N x 8] uint32 linear indices of the 8 corners
%     weight8 : [N x 8] single trilinear weights (rows sum to 1)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_BUILD_WARP_TO_MNI, ICDM_BUILD_WARP_TO_NATIVE

X = dim(1); Y = dim(2); Z = dim(3);

x0=floor(x); x1=x0+1;
y0=floor(y); y1=y0+1;
z0=floor(z); z1=z0+1;

% clipping
x0 = max(1,min(X,x0)); x1=max(1,min(X,x1));
y0 = max(1,min(Y,y0)); y1=max(1,min(Y,y1));
z0 = max(1,min(Z,z0)); z1=max(1,min(Z,z1));

% *** FIXED: row concatenation to get [Nout × 8] ***
index8 = uint32([
    sub2ind(dim,x0,y0,z0), ...
    sub2ind(dim,x1,y0,z0), ...
    sub2ind(dim,x0,y1,z0), ...
    sub2ind(dim,x1,y1,z0), ...
    sub2ind(dim,x0,y0,z1), ...
    sub2ind(dim,x1,y0,z1), ...
    sub2ind(dim,x0,y1,z1), ...
    sub2ind(dim,x1,y1,z1)
]);

% weights also [Nout × 8]
wx = x-x0; wx0=1-wx;
wy = y-y0; wy0=1-wy;
wz = z-z0; wz0=1-wz;

weight8 = single([
    wx0.*wy0.*wz0, ...
    wx.*wy0.*wz0, ...
    wx0.*wy.*wz0, ...
    wx.*wy.*wz0, ...
    wx0.*wy0.*wz, ...
    wx.*wy0.*wz, ...
    wx0.*wy.*wz, ...
    wx.*wy.*wz
]);
end
