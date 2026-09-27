function warp = icdm_build_warp_to_mni(nativefile, def_field,templatefile)
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
