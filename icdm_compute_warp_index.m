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

function [index8, weight8] = compute_warp_index(def_field, templatefile, dim_out, direction)
% COMPUTE_WARP_INDEX  Compute 8-neighbor trilinear interpolation indices and weights.
%
%   [index8, weight8] = compute_warp_index(def_field, templatefile, dim_out, direction)
%
%   Evaluates the SPM/DARTEL deformation field on the output grid and
%   returns, for each output voxel, the eight nearest-neighbor linear
%   indices and trilinear interpolation weights in the source grid.
%   This is an alternative implementation to icdm_trilinear_index_weight
%   that operates directly on the deformation field rather than on
%   pre-warped coordinate arrays.
%
%   Inputs
%     def_field    : path to DARTEL deformation field (y_*.nii).
%     templatefile : path to template NIfTI for world-coordinate mapping.
%     dim_out      : [1 x 3] output grid dimensions.
%     direction    : 'forward' (MNI to native) or 'inverse' (native to MNI).
%
%   Outputs
%     index8  : [N_out x 8] uint32 linear indices of the 8 neighbors.
%     weight8 : [N_out x 8] single trilinear weights.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_TRILINEAR_INDEX_WEIGHT, ICDM_BUILD_WARP_TO_MNI

    % Load template deformation flow (SPM/DARTEL)
    Td = spm_vol(templatefile);

    % Output grid coordinates
    [X,Y,Z] = ndgrid(1:dim_out(1), 1:dim_out(2), 1:dim_out(3));

    % Convert to world coordinates
    XYZ = [X(:)'; Y(:)'; Z(:)'; ones(1,numel(X))];
    XYZ_world = Td.mat * XYZ;

    % Deform
    d = spm_deformations(struct('comp',{def_field}, 'space',templatefile));

    if strcmp(direction,'forward')
        XYZ_def = d.mat * XYZ_world;    % forward
    else
        XYZ_def = d.invmat * XYZ_world; % inverse warp
    end

    % Convert to voxel coordinates in source grid
    M = inv(Td.mat);
    XYZ_src = M * [XYZ_def; ones(1,size(XYZ_def,2))];

    % Compute tri-linear interpolation weights
    x = XYZ_src(1,:)';  y = XYZ_src(2,:)';  z = XYZ_src(3,:)';

    x0 = floor(x); x1 = x0+1;
    y0 = floor(y); y1 = y0+1;
    z0 = floor(z); z1 = z0+1;

    % Clip out-of-bound
    x0 = max(1,min(dim_out(1),x0));
    y0 = max(1,min(dim_out(2),y0));
    z0 = max(1,min(dim_out(3),z0));

    x1 = max(1,min(dim_out(1),x1));
    y1 = max(1,min(dim_out(2),y1));
    z1 = max(1,min(dim_out(3),z1));

    % 8 neighbor indices
    index8 = uint32([
        sub2ind(dim_out,x0,y0,z0), ...
        sub2ind(dim_out,x1,y0,z0), ...
        sub2ind(dim_out,x0,y1,z0), ...
        sub2ind(dim_out,x1,y1,z0), ...
        sub2ind(dim_out,x0,y0,z1), ...
        sub2ind(dim_out,x1,y0,z1), ...
        sub2ind(dim_out,x0,y1,z1), ...
        sub2ind(dim_out,x1,y1,z1)
    ]);

    % weights
    wx = x - x0;  wx0 = 1-wx;  wx1 = wx;
    wy = y - y0;  wy0 = 1-wy;  wy1 = wy;
    wz = z - z0;  wz0 = 1-wz;  wz1 = wz;

    weight8 = single([
        wx0.*wy0.*wz0, ...
        wx1.*wy0.*wz0, ...
        wx0.*wy1.*wz0, ...
        wx1.*wy1.*wz0, ...
        wx0.*wy0.*wz1, ...
        wx1.*wy0.*wz1, ...
        wx0.*wy1.*wz1, ...
        wx1.*wy1.*wz1
    ]);

end
