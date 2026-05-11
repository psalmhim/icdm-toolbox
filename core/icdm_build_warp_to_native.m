function warp = icdm_build_warp_to_native(nativefile, def_field,templatefile)
% ICDM_BUILD_WARP_TO_NATIVE  Build a trilinear warp operator from MNI to native space.
%
%   warp = icdm_build_warp_to_native(nativefile, def_field, templatefile)
%
%   Creates a precomputed warp operator that maps data from MNI (template)
%   space back to native (DWI) space using trilinear interpolation.  An
%   identity coordinate grid in MNI space is warped inversely through the
%   DARTEL deformation field via icdm_warp_4d, and the resulting
%   continuous native-space coordinates are decomposed into 8-neighbor
%   indices and trilinear weights by icdm_trilinear_index_weight.
%
%   Inputs
%     nativefile   : path to a NIfTI in native space (for geometry).
%     def_field    : path to DARTEL deformation field (y_*.nii).
%     templatefile : path to DARTEL template NIfTI (defines MNI grid).
%
%   Output
%     warp : struct with fields
%              .index   [N_nat x 8] uint32 neighbor indices in MNI grid.
%              .weight  [N_nat x 8] single trilinear weights.
%              .dim_in  [1 x 3] MNI volume dimensions (input space).
%              .dim_out [1 x 3] native volume dimensions (output space).
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_BUILD_WARP_TO_MNI, ICDM_TRILINEAR_INDEX_WEIGHT,
%            ICDM_APPLY_WARP

% MNI grid (based on template affine)
Vtpl = spm_vol(templatefile);
dim_mni= Vtpl(1).dim(1:3);
[X,Y,Z] = ndgrid(1:dim_mni(1),1:dim_mni(2),1:dim_mni(3));

Yin = zeros([dim_mni 3],'single');
Yin(:,:,:,1) = X;
Yin(:,:,:,2) = Y;
Yin(:,:,:,3) = Z;

% ground-truth warp (SPM)
Ycoord = icdm_warp_4d(Yin, nativefile, def_field, templatefile, -1, 1);

Vn = spm_vol(nativefile);Vn=Vn(1);
dim_nat = Vn.dim(1:3);

if ~isequal(size(Ycoord,1:3), dim_nat)
    error('icdm_build_warp_to_native: Ycoord size mismatch.');
end

x = Ycoord(:,:,:,1); x = x(:);
y = Ycoord(:,:,:,2); y = y(:);
z = Ycoord(:,:,:,3); z = z(:);

[index8, weight8] = icdm_trilinear_index_weight(x,y,z,dim_mni);

warp.index   = index8;
warp.weight  = weight8;
warp.dim_in  = dim_mni;
warp.dim_out = dim_nat;
end
