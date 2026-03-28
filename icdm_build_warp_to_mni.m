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

function warp = icdm_build_warp_to_mni(nativefile, def_field,templatefile)
% ICDM_BUILD_WARP_TO_MNI  Build a trilinear warp operator from native to MNI space.
%
%   warp = icdm_build_warp_to_mni(nativefile, def_field, templatefile)
%
%   Creates a precomputed warp operator that maps data from native (DWI)
%   space to MNI (template) space using trilinear interpolation.  An
%   identity coordinate grid in native space is warped forward through the
%   DARTEL deformation field via icdm_warp_4d, and the resulting
%   continuous MNI-space coordinates are decomposed into 8-neighbor
%   indices and trilinear weights by icdm_trilinear_index_weight.
%
%   Inputs
%     nativefile   : path to a NIfTI in native space (for geometry).
%     def_field    : path to DARTEL deformation field (y_*.nii).
%     templatefile : path to DARTEL template NIfTI (defines MNI grid).
%
%   Output
%     warp : struct with fields
%              .index   [N_mni x 8] uint32 neighbor indices in native grid.
%              .weight  [N_mni x 8] single trilinear weights.
%              .dim_in  [1 x 3] native volume dimensions.
%              .dim_out [1 x 3] MNI volume dimensions.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_BUILD_WARP_TO_NATIVE, ICDM_TRILINEAR_INDEX_WEIGHT,
%            ICDM_APPLY_WARP
Vn = spm_vol(nativefile); Vn=Vn(1);
dim_nat = Vn.dim(1:3);

[X,Y,Z] = ndgrid(1:dim_nat(1),1:dim_nat(2),1:dim_nat(3));

Yin = zeros([dim_nat 3],'single');
Yin(:,:,:,1) = X;
Yin(:,:,:,2) = Y;
Yin(:,:,:,3) = Z;

% ground truth warp
Ycoord = icdm_warp_4d(Yin, nativefile, def_field, templatefile, +1, 1);

Vtpl = spm_vol(templatefile);
dim_mni= Vtpl(1).dim(1:3);

if ~isequal(size(Ycoord,1:3), dim_mni)
    error('icdm_build_warp_to_mni: Ycoord size mismatch.');
end

x = Ycoord(:,:,:,1); x = x(:);
y = Ycoord(:,:,:,2); y = y(:);
z = Ycoord(:,:,:,3); z = z(:);

[index8, weight8] = icdm_trilinear_index_weight(x,y,z,dim_nat);

warp.index   = index8;
warp.weight  = weight8;
warp.dim_in  = dim_nat;
warp.dim_out = dim_mni;
end
