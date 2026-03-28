% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
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

function icdm_plot_ilr_vs_pi(ilr_vols, d, H)
% ICDM_PLOT_ILR_VS_PI  Display a single ILR basis dimension at mid-axial slice.
%
%   icdm_plot_ilr_vs_pi(ilr_vols, d, H)
%
%   Shows the d-th ILR coordinate volume at the mid-axial slice.  Note
%   that ILR basis dimensions do not correspond one-to-one to parcels.
%
%   INPUT
%     ilr_vols : [X x Y x Z x (K-1)] 4-D ILR volume array
%     d        : ILR dimension to display
%     H        : [K x (K-1)] Helmert submatrix (reserved for future use)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PLOT_ILR_MOSAIC, ICDM_ILR_TO_VOLUMES

% Show ILR(d) and all pi parcels separately
Z = ilr_vols(:,:,:,d);   % a single ILR basis
K1 = size(H,2);
K = K1+1;

% Just display ILR(d)
figure; imagesc(rot90(Z(:,:,round(end/2)))); axis off; 
title(sprintf('ILR dimension %d', d));

% Not showing π(d) since ILR basis doesn’t correspond 1:1 to parcels.
% But I can compute π from the entire mu_ilr if needed.
