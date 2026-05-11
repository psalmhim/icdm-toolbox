function icdm_plot_ilr_mosaic(ilr_vols, d, nslices)
% ICDM_PLOT_ILR_MOSAIC  Axial mosaic display of one ILR dimension.
%
%   icdm_plot_ilr_mosaic(ilr_vols, d, nslices)
%
%   Displays equally-spaced axial slices of the d-th ILR basis volume
%   from a 4-D ILR volume array in a tiled mosaic layout.
%
%   INPUT
%     ilr_vols : [X x Y x Z x (K-1)] 4-D ILR volume array
%     d        : ILR dimension to display (1-based)
%     nslices  : number of axial slices (default 20)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PLOT_MOSAIC, ICDM_ILR_TO_VOLUMES

if nargin < 3
    nslices = 20;
end

vol = ilr_vols(:,:,:,d);
dims = size(vol);
sl = round(linspace(1,dims(3),nslices));

figure('Color','w','Position',[100 50 900 900]);
tiledlayout(ceil(nslices/5),5,'Padding','compact','TileSpacing','compact');

for i = 1:nslices
    nexttile;
    imagesc(rot90(vol(:,:,sl(i))));
    axis off equal;
    title(sprintf('z=%d',sl(i)));
end
sgtitle(sprintf('ILR dimension %d', d));
colormap turbo;
colorbar;
end
