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

function icdm_plot_ci_heatmap(Ymap, Qdiag, H, T1, mask, out_png)
% ICDM_PLOT_CI_HEATMAP  Overlay posterior uncertainty heatmap on anatomy.
%
%   icdm_plot_ci_heatmap(Ymap, Qdiag, H, T1, mask, out_png)
%
%   Computes the mean pi-space standard deviation at each voxel by
%   propagating ILR posterior precision through the softmax Jacobian,
%   then overlays the result on a background T1 slice and saves to PNG.
%
%   INPUT
%     Ymap    : [nV x (K-1)] ILR MAP estimates
%     Qdiag   : [nV x (K-1)] diagonal posterior precision
%     H       : [K x (K-1)] Helmert submatrix
%     T1      : [H x W] background anatomical slice
%     mask    : [H x W] logical voxel mask
%     out_png : output PNG file path
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_CONFIDENCE_INTERVAL, ICDM_CI_RECOVERY

[nV, Km1] = size(Ymap);
K = Km1 + 1;
Hmat = H;

std_map = zeros(nV,1,'single');

for v = 1:nV
    y = double(Ymap(v,:)');
    prec = double(Qdiag(v,:)');
    Cov_y = diag(1./max(prec,1e-12));

    z = Hmat*y;
    pi = softmax(z);

    J = (diag(pi) - pi*pi') * Hmat;
    Cov_pi = J * Cov_y * J';
    std_pi = sqrt(max(diag(Cov_pi),1e-12));
    std_map(v) = mean(std_pi);
end

% reshape back to 2D
std2 = zeros(size(mask),'single');
std2(mask) = std_map;

figure('Color','w','Position',[100 100 800 600]);
imagesc(rot90(T1,1)); colormap gray; hold on;
h = imagesc(rot90(std2,1)); 
set(h,'AlphaData',rot90(mask>0,1)*0.85);
title('Posterior Uncertainty (π-space STD)');
axis image off;
colorbar;

exportgraphics(gcf, out_png, 'Resolution',150);
close(gcf);

end

function pi = softmax(z)
    a = exp(z - max(z));
    pi = a / sum(a);
end
