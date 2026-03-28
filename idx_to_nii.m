% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI", NeuroImage (2026).
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
