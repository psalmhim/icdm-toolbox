function [ci_low, ci_high] = icdm_confidence_interval(Ymap, Qdiag, H)
% ICDM_CONFIDENCE_INTERVAL  Voxelwise 95% confidence interval on pi-space.
%
%   [ci_low, ci_high] = icdm_confidence_interval(Ymap, Qdiag, H)
%
%   Propagates posterior ILR uncertainty through the softmax Jacobian to
%   obtain approximate 95% Wald confidence intervals for each compositional
%   proportion pi_k at every voxel.
%
%   INPUT
%     Ymap  : [nV x (K-1)] ILR MAP estimates
%     Qdiag : [nV x (K-1)] diagonal posterior precision per voxel
%     H     : [K x (K-1)] Helmert submatrix
%
%   OUTPUT
%     ci_low  : [nV x K] lower bound of 95% CI for pi
%     ci_high : [nV x K] upper bound of 95% CI for pi
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_GROUP_CI, ICDM_PLOT_CI_HEATMAP, ICDM_CI_RECOVERY

[nV,Km1] = size(Ymap);
K = Km1 + 1;

ci_low  = zeros(nV,K,'single');
ci_high = zeros(nV,K,'single');

for v = 1:nV
    y = double(Ymap(v,:)');
    prec = double(Qdiag(v,:)');
    Cov_y = diag(1 ./ max(prec,1e-12));   % posterior covariance

    z = H * y;
    pi = softmax(z);

    J = (diag(pi) - pi*pi') * H;  % K×(K−1)
    Cov_pi = J * Cov_y * J';      % K×K

    std_pi = sqrt(max(diag(Cov_pi), 1e-12));

    ci_low(v,:)  = single(pi' - 1.96 * std_pi');
    ci_high(v,:) = single(pi' + 1.96 * std_pi');
end
end

function pi = softmax(z)
    mz = max(z);
    a = exp(z - mz);
    pi = a / sum(a);
end
