function icdm_plot_ilr_3panel(group_icdm_file, parcel_idx, ilr_idx)
% ICDM_PLOT_ILR_3PANEL  Three-panel view: T1, ILR dimension, and pi parcel.
%
%   icdm_plot_ilr_3panel(group_icdm_file, parcel_idx, ilr_idx)
%
%   Loads group iCDM results and displays a three-panel figure at the
%   mid-axial slice showing:
%     (1) SPM canonical T1 anatomy
%     (2) Selected ILR basis dimension
%     (3) Compositional proportion pi for the specified parcel
%
%   INPUT
%     group_icdm_file : MAT file with group iCDM results
%     parcel_idx      : parcel index for pi display
%     ilr_idx         : ILR dimension index for ILR display
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PLOT_3PANEL, ICDM_ILR_TO_VOLUMES

G = load(group_icdm_file);

% ILR → full vol
ilr_vols = icdm_ilr_to_volumes(G.mu_ilr_mni, G.idx_mni, G.dim_mni);

% π
H = G.H;
K1 = size(G.mu_ilr_mni,2);
K = K1+1;

Z = G.mu_ilr_mni * H.';
Z = Z - max(Z,[],2);
A = exp(Z);
PI = A ./ sum(A,2);

VOL = zeros(prod(G.dim_mni),1,'single');
VOL(G.idx_mni) = PI(:,parcel_idx);
pi_vol = reshape(VOL,G.dim_mni);

% ILR dimension d
ilr_vol = ilr_vols(:,:,:,ilr_idx);

% T1
Vt1 = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
t1 = spm_read_vols(Vt1);

mid = round(size(t1,3)/2);

figure('Color','w','Position',[100 100 1200 400]);

subplot(1,3,1);
imagesc(rot90(t1(:,:,mid))); axis off;
title('T1');

subplot(1,3,2);
imagesc(rot90(ilr_vol(:,:,mid))); axis off;
title(sprintf('ILR dim %d', ilr_idx));
colormap(gca,'turbo'); colorbar;

subplot(1,3,3);
imagesc(rot90(pi_vol(:,:,mid))); axis off;
title(sprintf('\\pi parcel %d', parcel_idx));
colormap(gca,'turbo'); colorbar;

end
