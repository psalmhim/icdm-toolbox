function icdm_check_group_iter(group_path)
% ICDM_CHECK_GROUP_ITER  Materialize group-iteration results as NIfTI files for QC.
%
%   icdm_check_group_iter(group_path)
%
%   Loads beta_ilr.mat and yilr_group_stat.mat from the specified
%   iteration directory and writes the following NIfTI volumes for
%   visual quality-control inspection:
%     - Per-predictor beta maps (beta_ilr_beta*.nii)
%     - Group ILR mean (yilr_group_mean.nii)
%     - Group ILR std  (yilr_group_std.nii)
%     - Group kappa    (yilr_group_kappa.nii)
%     - Group pi       (group_pi.nii)  via ILR inverse transform
%
%   Input
%     group_path : path to the iteration folder (default pwd).
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_RESULT_TO_NII, ICDM_CHECK_GROUP
    if nargin<2, group_path=pwd;end
    load('vtemplate.mat')
    beta_idx_mat = fullfile(group_path, 'beta_ilr.mat');
    load(beta_idx_mat);
    beta_idx_nii = fullfile(group_path, 'beta_ilr.nii');
    materialize_beta_nii_from_idx(vtpl, Beta, idx_mni, beta_idx_nii);
    
    yilr_group_mat = fullfile(group_path, 'yilr_group_stat.mat');
    load(yilr_group_mat);
    
    yilr_group_mean=fullfile(group_path,'yilr_group_mean.nii');
    mat_to_nii(vtpl, y_mean, idx_mni,yilr_group_mean);
    yilr_group_std=fullfile(group_path,'yilr_group_std.nii');
    mat_to_nii(vtpl, y_std, idx_mni,yilr_group_std);
    yilr_group_kappa=fullfile(group_path,'yilr_group_kappa.nii');
    mat_to_nii(vtpl, kappa, idx_mni,yilr_group_kappa);

    Km1=size(y_mean,2);
    pi=ilr_inverse(y_mean,Km1+1);
    group_pi=fullfile(group_path,'group_pi.nii');
    mat_to_nii(vtpl, pi, idx_mni,group_pi);

    group_kappa=fullfile(group_path,'group_kappa.nii');
    mat_to_nii(vtpl, sum(pi,2), idx_mni,group_kappa);
end


function mat_to_nii(Vtpl, data, idx_mni,out_path,descrip)
if nargin<5
    [p,descrip,e]=fileparts(out_path);
end

[m, n] = size(data); Km1=n; if n>m, data=data'; Km1=m;end
Vout = repmat(Vtpl,Km1, 1);
for f=1:Km1
    Vout(f).fname   = out_path;
    Vout(f).n       = [f 1];
    Vout(f).dt      = [spm_type('float32') 0];
    Vout(f).pinfo   = [1;0;0];
    Vout(f).descrip = descrip;
    vol = zeros(Vout(f).dim(1:3),'single');
    vol(idx_mni) = single(data(:,f));
    spm_write_vol(Vout(f), reshape(vol,Vout(f).dim(1:3)));
end
end

function materialize_beta_nii_from_idx(Vtpl, Beta, idx_mni,out_path)
[P, ~, Km1] = size(Beta);
[p1,f1,e1]=fileparts(out_path);
for p=1:P
    out_path1 = fullfile(p1,sprintf('%s_beta%d.nii',f1,p));
    Vout = repmat(Vtpl,Km1, 1);
    for f=1:Km1
        Vout(f).fname   = out_path1;
        Vout(f).n       = [f 1];
        Vout(f).dt      = [spm_type('float32') 0];
        Vout(f).pinfo   = [1;0;0];
        Vout(f).descrip = sprintf('β (ILR) p=%d, d=%d', p, f);
        vol = zeros(Vout(f).dim(1:3),'single');
        vol(idx_mni) = single(Beta(p,:,f));
        spm_write_vol(Vout(f), reshape(vol,Vout(f).dim(1:3)));
    end
end
end
