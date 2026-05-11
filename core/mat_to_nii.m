function vo=mat_to_nii(matfile,data,fout)
% MAT_TO_NII  Write index-based data to a NIfTI volume using a template.
%
%   vo = mat_to_nii(matfile, data, fout)
%
%   Loads MNI mask indices (idx_mni) and volume header (v) from matfile,
%   fills a zero volume at those indices with data, and writes the result
%   as a NIfTI file.  Opens the image in the SPM display viewer.
%
%   INPUT
%     matfile : MAT file containing idx_mni and v (SPM volume header)
%     data    : vector of values to place at idx_mni locations
%     fout    : (optional) output NIfTI path; defaults to matfile with .nii
%
%   OUTPUT
%     vo : SPM volume header of the written file
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also IDX_TO_NII, ICDM_RESULT_TO_NII
load('vtemplate.mat');

if nargin<3,
    [p,f,e]=fileparts(matfile); fout=fullfile(p,[f '.nii']);
end

load(matfile);
R=zeros(v.dim);
R(idx_mni)=data;
vo=v;
vo.fname=fout;
vo=spm_write_vol(vo,R);
spm_image('init',vo.fname);
end
