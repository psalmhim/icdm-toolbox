function vo=idx_to_nii(dim_mni,idx_mni,data,v,fout)
% IDX_TO_NII  Write index-based data into a NIfTI volume and display.
%
%   vo = idx_to_nii(dim_mni, idx_mni, data, v, fout)
%
%   Places data at the linear indices idx_mni within a zero volume of
%   size dim_mni, writes it as a NIfTI file using the SPM volume header
%   v, and opens the result in the SPM image viewer.
%
%   INPUT
%     dim_mni : [1 x 3] volume dimensions
%     idx_mni : [N x 1] linear indices
%     data    : [N x 1] values to write
%     v       : SPM volume header (template)
%     fout    : output NIfTI filename
%
%   OUTPUT
%     vo : SPM volume header of the written file
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also MAT_TO_NII, ICDM_RESULT_TO_NII
R=zeros(dim_mni);
R(idx_mni)=data;
vo=v;
vo.fname=fout;
vo=spm_write_vol(vo,R);
spm_image('init',vo.fname);
end
