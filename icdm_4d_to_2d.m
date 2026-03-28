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

function R2d=icdm_4d_to_2d(R4d,idx)
% ICDM_4D_TO_2D  Extract index-form rows from a 4-D (or 3-D) volume.
%
%   R2d = icdm_4d_to_2d(R4d, idx)
%
%   Reshapes R4d into [prod(spatial dims) x D] and extracts the rows at
%   linear indices idx, yielding a compact [numel(idx) x D] matrix.
%   This is the standard iCDM operation for converting between full-volume
%   and index-based representations.
%
%   Inputs
%     R4d : array of size [X Y Z D] (or [X Y Z] treated as D=Z).
%     idx : linear indices into the spatial grid.
%
%   Output
%     R2d : [numel(idx) x D] extracted data.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_2D_TO_4D
sz=size(R4d); D=sz(end);
R4d=reshape(R4d,[],D);
R2d=R4d(idx,:);
end
