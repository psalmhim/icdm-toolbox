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
