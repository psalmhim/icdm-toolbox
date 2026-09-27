function warp = icdm_build_warp_to_native(nativefile, def_field,templatefile)

% MNI grid (template affine 기준)
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
