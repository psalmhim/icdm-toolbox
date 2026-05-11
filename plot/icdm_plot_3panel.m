function icdm_plot_3panel(group_icdm_file, group_beta_file, parcel_idx, cov_idx, ilr_idx)
% ICDM_PLOT_3PANEL  Three-panel display: T1 anatomy, pi map, and beta map.
%
%   icdm_plot_3panel(group_icdm_file, group_beta_file, parcel_idx, cov_idx, ilr_idx)
%
%   Displays a three-panel figure at the mid-axial slice showing:
%     (1) SPM canonical T1 anatomy
%     (2) Compositional proportion pi for the specified parcel
%     (3) Beta coefficient for the specified covariate and ILR dimension
%
%   INPUT
%     group_icdm_file : MAT file with group iCDM results (mu_ilr_mni, H, ...)
%     group_beta_file : MAT file with Beta field [P x Nmni x K1]
%     parcel_idx      : parcel index for pi display
%     cov_idx         : covariate index for beta display
%     ilr_idx         : ILR dimension index for beta display
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PLOT_ILR_3PANEL, ICDM_DISPLAY_ALL

cg = load(group_icdm_file);
cb = load(group_beta_file);

% -------------------------------------------------------------------------
% 1. ILR → π
% -------------------------------------------------------------------------
mu = cg.mu_ilr_mni;          % Nmni × (K-1)
K1 = size(mu,2);
K  = K1+1;
H  = cg.H;

Z = mu * H.';                % Nmni × K
Z = Z - max(Z,[],2);
A = exp(Z);
PI = A ./ sum(A,2);

% build π(k) volume
vol = zeros(prod(cg.dim_mni),1,'single');
vol(cg.idx_mni) = PI(:,parcel_idx);
pi_vol = reshape(vol, cg.dim_mni);

% build β(cov_idx,ilr_idx)
Beta = cb.Beta;              % [P × Nmni × K1]
beta_vec = squeeze(Beta(cov_idx,:,ilr_idx));
vol2 = zeros(prod(cg.dim_mni),1,'single');
vol2(cg.idx_mni) = beta_vec;
beta_vol = reshape(vol2, cg.dim_mni);

% canonical
Vt1 = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
t1  = spm_read_vols(Vt1);

% -------------------------------------------------------------------------
% 2. FIGURE
% -------------------------------------------------------------------------
figure('Color','w','Position',[100 100 1200 400]);

subplot(1,3,1);
imagesc(rot90(t1(:,:,round(end/2)))); axis off equal
title('T1 Anatomy');

subplot(1,3,2);
imagesc(rot90(pi_vol(:,:,round(end/2)))); axis off equal
title(sprintf('\\pi parcel %d', parcel_idx));
colorbar

subplot(1,3,3);
imagesc(rot90(beta_vol(:,:,round(end/2)))); axis off equal
title(sprintf('\\beta cov=%d, ilr=%d', cov_idx, ilr_idx));
colorbar

end
