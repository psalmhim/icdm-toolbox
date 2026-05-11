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
