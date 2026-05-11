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

