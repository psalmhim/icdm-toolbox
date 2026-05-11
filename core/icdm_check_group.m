function icdm_check_group(group_path)
% ICDM_CHECK_GROUP  Materialise group beta coefficients as NIfTI volumes.
%
%   icdm_check_group(group_path)
%
%   Loads the compact index-based beta file (beta_ilr.mat) from group_path
%   and writes per-covariate 4-D NIfTI volumes of ILR beta coefficients
%   for visual inspection.
%
%   INPUT
%     group_path : directory containing beta_ilr.mat
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_CHECK_GROUP_ITER, ICDM_ESTIMATE_BETA

% materialize_beta_nii(beta_idx_mat, beta_nii)
beta_idx_mat = fullfile(group_path, 'beta_ilr.mat');

load(beta_idx_mat);
materialize_beta_nii_from_idx(Beta, idx_mni, dim_mni, Vtpl, beta_nii);

end


function materialize_beta_nii_from_idx(Beta, idx_mni, dim_mni, Vtpl, out_path)
[P, ~, Km1] = size(Beta);
F = P*Km1;
Vout = repmat(Vtpl, F, 1);
[p1,f1,e1]=fileparts(out_path);
for p=1:P
    out_path1 = fullfile(p1,sprintf('%s_beta%d.nii',f1,p));
    for f=1:Km1
        Vout(f).fname   = out_path1;
        Vout(f).n       = [f 1];
        Vout(f).dt      = [spm_type('float32') 0];
        Vout(f).pinfo   = [1;0;0];
        Vout(f).descrip = sprintf('β (ILR) p=%d, d=%d', p, d);
        vol = zeros(prod(dim_mni),1,'single');
        vol(idx_mni) = single(Beta(p,:,d));
        spm_write_vol(Vout(f), reshape(vol,dim_mni));
    end
end
end
