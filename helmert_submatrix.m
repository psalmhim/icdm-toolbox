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

function H = helmert_submatrix(K)
% HELMERT_SUBMATRIX  Orthonormal Helmert contrast matrix for ILR transform.
%
%   H = helmert_submatrix(K)
%
%   Constructs a K x (K-1) matrix H such that H'*H = I and H'*1 = 0,
%   mapping K-part compositions to (K-1)-dimensional ILR coordinates.
%   The ILR forward transform is  y = H' * clr(pi)  and the inverse is
%   pi = softmax(H * y).  (Manuscript Eqs. ilr, softmax)
%
%   INPUT
%     K : number of compositional components (e.g. 68 cortical targets)
%
%   OUTPUT
%     H : [K x (K-1)] orthonormal Helmert submatrix
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ILR_INVERSE, ICDM_SUBJECT_VB

H = zeros(K, K-1);
for i=1:(K-1)
    H(1:i, i) =  1 / sqrt(i*(i+1));
    H(i+1, i) = -i / sqrt(i*(i+1));
end
end

