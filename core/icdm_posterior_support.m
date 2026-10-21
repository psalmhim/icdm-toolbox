function PS = icdm_posterior_support(subj_datafile, vb, grp, alpha_grp, targets)
% ICDM_POSTERIOR_SUPPORT  Template-space support indices for one subject.
%
%   PS = icdm_posterior_support(subj_datafile, vb, grp, alpha_grp, targets)
%
% Inputs
%   subj_datafile : subject cache written by icdm_compose_subject (icdm2d, idx_native, warp)
%   vb            : subject output of icdm_subject_vb (struct or <id>_vb.mat path) with
%                   y_ilr_mni, kappa_data_mni, kappa_mni, valid_mni
%   grp           : population prior (struct or group_icdm_iter_XXX.mat path) with idx_mni, kappa_mni, H
%                   (grp.kappa_mni is lambda_grp; the prior precision applied is alpha_grp * lambda_grp)
%   alpha_grp     : population-transfer strength used for inference (e.g. 0.25)
%   targets       : target channels entering the model (e.g. 1:68)
%
% Output (all [N_mni x 1], in grp.idx_mni order)
%   PS.N_mni      : streamline--target incidences N_v (sum over targets), warped with the subject's operator
%   PS.kappa_data : prior-independent data information (likelihood curvature at the relative frequencies)
%   PS.kappa_post : posterior certainty (mean posterior curvature at the MAP)
%   PS.r_data     : posterior data fraction tr(Q_data)/tr(Q_post) at the MAP,
%                   N c(pi_hat) / (N c(pi_hat) + alpha_grp * lambda_grp); 0 where N_v = 0
%   PS.valid      : subject validity mask in template space (valid_mni)
if ischar(vb), vb = load(vb, 'y_ilr_mni','kappa_data_mni','kappa_mni','valid_mni'); end
if ischar(grp), grp = load(grp, 'idx_mni','kappa_mni','H'); end
c = load(subj_datafile, 'icdm2d','warp','idx_native');
n = double(c.icdm2d(:, targets)); n(~isfinite(n) | n < 0) = 0;
N_mni = double(icdm_warp_to_mni(single(sum(n,2)), c.warp, grp.idx_mni, c.idx_native));
N_mni = N_mni(:);

H = double(grp.H); if size(H,1) < size(H,2), H = H.'; end
Y = double(vb.y_ilr_mni); if size(Y,2) ~= size(H,2), Y = Y.'; end
P = exp(Y * H.'); P = P ./ sum(P, 2);
cpi = mean(P * (H.^2) - (P * H).^2, 2);          % mean diagonal of H' (diag(pi) - pi pi') H

lam = double(grp.kappa_mni(:));
num = N_mni .* cpi;
r = num ./ max(num + alpha_grp * lam, 1e-30);
r(N_mni <= 0) = 0;
valid = logical(vb.valid_mni(:));
r(~valid) = NaN;

PS.N_mni      = N_mni;
PS.kappa_data = double(vb.kappa_data_mni(:));
PS.kappa_post = double(vb.kappa_mni(:));
PS.r_data     = r;
PS.valid      = valid;
end
