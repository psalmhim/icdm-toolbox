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

function icdm_plot_mosaic(vol, title_str, nslices)
% ICDM_PLOT_MOSAIC  Axial mosaic display of a 3-D volume.
%
%   icdm_plot_mosaic(vol, title_str, nslices)
%
%   Displays equally-spaced axial slices of a 3-D volume (e.g. pi or beta)
%   in a tiled mosaic layout with a shared colorbar.
%
%   INPUT
%     vol       : 3-D volume array
%     title_str : figure super-title string
%     nslices   : number of slices to show (default 20)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PLOT_ILR_MOSAIC, ICDM_DISPLAY_ALL

if nargin < 3
    nslices = 20;
end

dims = size(vol);
sl = round(linspace(1,dims(3),nslices));

figure('Color','w','Position',[50 50 1000 800]);
tiledlayout(ceil(nslices/5),5,'Padding','compact','TileSpacing','compact');

for i = 1:nslices
    nexttile;
    imagesc(rot90(vol(:,:,sl(i)))); axis off equal;
    title(sprintf('z=%d',sl(i)));
end

sgtitle(title_str);
colormap parula; colorbar;
end

