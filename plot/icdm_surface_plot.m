function icdm_surface_plot(pi_vol, hemi, surf_path, out_png)
% ICDM_SURFACE_PLOT  Map a volumetric pi map onto a cortical surface mesh.
%
%   icdm_surface_plot(pi_vol, hemi, surf_path, out_png)
%
%   Projects a 3-D volumetric compositional proportion map onto a GIfTI
%   cortical surface mesh by nearest-voxel lookup, then renders the result
%   as a 3-D surface patch with Gouraud lighting.
%
%   INPUT
%     pi_vol    : 3-D volume of compositional proportions
%     hemi      : hemisphere label ('lh' or 'rh')
%     surf_path : path to GIfTI surface file (.surf.gii)
%     out_png   : (optional) output PNG file path
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_DISPLAY_ALL, ICDM_PLOT_MOSAIC

G = gifti(surf_path);
coords = G.vertices;      % Nx3
faces  = G.faces;

% MNI voxel coordinates
[X,Y,Z] = ndgrid(1:size(pi_vol,1), 1:size(pi_vol,2), 1:size(pi_vol,3));
Vxyz = [X(:), Y(:), Z(:)];

% convert voxel->world using SPM affine
Vtmp = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
M = Vtmp.mat;
Wxyz = (M * [Vxyz'; ones(1,size(Vxyz,1))])';  % world coords

% map surface → nearest voxel
pi_flat = pi_vol(:);
surf_pi = zeros(size(coords,1),1);

for v = 1:size(coords,1)
    [~,idx] = min(sum((Wxyz(:,1:3)-coords(v,:)).^2,2));
    surf_pi(v) = pi_flat(idx);
end

% 3D plot
figure('Color','w');
patch('Faces',faces,'Vertices',coords,...
      'FaceVertexCData',surf_pi,...
      'FaceColor','interp','EdgeColor','none');
axis equal off;
camlight headlight; lighting gouraud;
title(sprintf('%s surface π',hemi));

if nargin > 3
    exportgraphics(gcf,out_png,'Resolution',300);
end
end

