function make_identity_deformation(template_file, out_file)
% MAKE_IDENTITY_DEFORMATION  SPM deformation field y_*.nii that maps the template grid onto itself.
% A deformation stores, for every template voxel, the world (mm) coordinate it samples in subject space;
% the identity field therefore holds each template voxel's own world coordinate.
V = spm_vol(template_file); V = V(1); dim = V.dim(1:3);
[i,j,k] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
xyz = V.mat * [i(:)'; j(:)'; k(:)'; ones(1, numel(i))];
N = nifti; N.dat = file_array(out_file, [dim 1 3], 'float32-le', 0, 1, 0);
N.mat = V.mat; N.mat0 = V.mat; N.mat_intent = 'Aligned'; N.mat0_intent = 'Aligned';
N.descrip = 'Deformation (identity)'; create(N);
for d = 1:3, N.dat(:,:,:,1,d) = reshape(single(xyz(d,:)), dim); end
end
