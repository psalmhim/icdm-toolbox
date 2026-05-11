function [CI_low, CI_high, std_pi] = icdm_group_ci(y_mean, y_std, H)
% ICDM_GROUP_CI  Group-level 95% confidence interval on pi-space.
%
%   [CI_low, CI_high, std_pi] = icdm_group_ci(y_mean, y_std, H)
%
%   Propagates group ILR variance (from trimmed standard deviations)
%   through the softmax Jacobian to compute approximate 95% Wald
%   confidence intervals for the compositional proportions pi_k.
%
%   INPUT
%     y_mean : [nV x (K-1)] group mean ILR coordinates
%     y_std  : [nV x (K-1)] group trimmed standard deviations
%     H      : [K x (K-1)] Helmert submatrix
%
%   OUTPUT
%     CI_low  : [nV x K] lower bound of 95% CI for pi
%     CI_high : [nV x K] upper bound of 95% CI for pi
%     std_pi  : [nV x K] standard deviation of pi
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_CONFIDENCE_INTERVAL, ICDM_UPDATE_GROUP_PRIOR

[nV,Km1] = size(y_mean);
K = Km1 + 1;

CI_low  = zeros(nV,K,'single');
CI_high = zeros(nV,K,'single');
std_pi  = zeros(nV,K,'single');

for v = 1:nV
    mu = y_mean(v,:)';
    s  = y_std(v,:)';
    Cov_y = diag(s.^2);

    z = H*mu;
    pi = softmax(z);

    J = (diag(pi) - pi*pi') * H;
    Cov_pi = J * Cov_y * J';

    stdp = sqrt(max(diag(Cov_pi),1e-12));
    CI_low(v,:)  = pi' - 1.96 * stdp';
    CI_high(v,:) = pi' + 1.96 * stdp';
    std_pi(v,:)  = stdp';
end
end
