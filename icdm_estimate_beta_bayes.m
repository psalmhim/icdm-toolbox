% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI", NeuroImage (2026).
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

function [beta_mat, Beta_field, Var_field, Z_field, PP_field] = ...
    icdm_estimate_beta_bayes(subjects, OUTS, K, opts)
% ICDM_ESTIMATE_BETA_BAYES  Bayesian voxelwise ridge regression for ILR beta fields.
%
%   [beta_mat, Beta_field, Var_field, Z_field, PP_field] = ...
%       icdm_estimate_beta_bayes(subjects, OUTS, K, opts)
%
%   Performs Bayesian ridge regression at every MNI voxel, estimating
%   posterior mean and variance of the regression coefficients linking
%   subject covariates (design matrix) to the ILR-transformed
%   connectivity.  Weights combine per-subject posterior precision
%   (kappa^alpha) with voxel coverage (coverage^gamma).
%
%   The posterior is:
%     Q_v  = X' W_v X + lambda * I     (posterior precision)
%     mu_v = Q_v^{-1} X' W_v y         (posterior mean)
%     Sigma_v = Q_v^{-1}               (posterior covariance)
%
%   Posterior z-scores and probabilities Phi(z) are also computed for
%   each predictor and ILR dimension (manuscript Eq. for beta inference).
%
%   Inputs
%     subjects : [S x 1] struct array with property fields.
%     OUTS     : {S x 1} cell array of subject VB outputs, each with
%                y_ilr_mni [Nmni x K1] and kappa_mni [Nmni x 1].
%     K        : number of compositional targets (parcels).
%     opts     : struct with design_spec, idx_mni, ridge_lambda,
%                beta_alpha, beta_gamma.
%
%   Outputs
%     beta_mat   : [K1 x (P-1)] group-mean beta (excluding intercept).
%     Beta_field : [P x Nmni x K1] posterior means.
%     Var_field  : [P x Nmni x K1] posterior variances.
%     Z_field    : [P x Nmni x K1] posterior z-scores.
%     PP_field   : [P x Nmni x K1] posterior probabilities Phi(z).
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_ESTIMATE_BETA, ICDM_DESIGN_VECTOR
%
% -------------------------------------------------------------------------
fprintf('=== Estimating voxelwise β (BAYESIAN) ===\n');

%% Dimensions
S    = numel(subjects);
K1   = K - 1;
Nmni = numel(opts.idx_mni);
P    = numel(opts.design_spec);

%% -----------------------------------------------------------------------
%% 1. DESIGN MATRIX  (S × P)
%% -----------------------------------------------------------------------
X = zeros(S, P, 'single');
for s = 1:S
    X(s,:) = single(icdm_design_vector(subjects(s), opts));
end
Xd = double(X);
clear X

%% -----------------------------------------------------------------------
%% 2. STACK Y, KAPPA
%% -----------------------------------------------------------------------
Yall = zeros(S, Nmni, K1, 'single');
Kap  = zeros(S, Nmni, 'single');

for s = 1:S
    Y = OUTS{s}.y_ilr_mni;
    Y(~isfinite(Y)) = 0;
    Yall(s,:,:) = single(Y);
    Kap(s,:)    = single(OUTS{s}.kappa_mni);
end

%% -----------------------------------------------------------------------
%% 3. COVERAGE
%% -----------------------------------------------------------------------
coverage = mean(Kap > 0, 1);     % [1 × Nmni]

%% -----------------------------------------------------------------------
%% 4. WEIGHT PARAMETERS
%% -----------------------------------------------------------------------
alpha  = getfield_default(opts,'beta_alpha',0.5);
gamma  = getfield_default(opts,'beta_gamma',0.2);
lambda = opts.ridge_lambda;

Ireg = lambda * eye(P);

%% -----------------------------------------------------------------------
%% 5. PRECOMPUTE POWERED WEIGHTS
%% -----------------------------------------------------------------------
Kap_alpha = max(double(Kap), eps).^alpha;
w_cov     = max(coverage, eps).^gamma;

%% -----------------------------------------------------------------------
%% 6. ALLOCATE RESULTS
%% -----------------------------------------------------------------------
Beta_field = zeros(P, Nmni, K1, 'single');
Var_field  = zeros(P, Nmni, K1, 'single');
Z_field    = zeros(P, Nmni, K1, 'single');
PP_field   = zeros(P, Nmni, K1, 'single');

%% -----------------------------------------------------------------------
%% 7. MAIN LOOP — Bayesian voxelwise regression
%% -----------------------------------------------------------------------
fprintf('  -- Running Bayesian regression for %d voxels\n', Nmni);

parfor v = 1:Nmni

    % Weights
    w  = Kap_alpha(:,v) * w_cov(v);
    sw = sqrt(w);

    % Weighted design
    Xsw = Xd .* sw;
    Xt  = Xsw';
    XtX = Xt*Xsw + Ireg;

    % Posterior covariance (Σ = Q^{-1})
    Sigma = inv(XtX);
    varB  = diag(Sigma);

    % Extract ILR responses
    Yv  = reshape(Yall(:,v,:), [S K1]);
    Ysw = bsxfun(@times, Yv, sw);    % S × K1

    % Posterior mean for all ILR dims
    B = XtX \ (Xt * Ysw);            % P × K1

    Beta_field(:,v,:) = single(B);
    Var_field(:,v,:)  = single(varB);

    % Compute z & posterior probability
    for d = 1:K1
        z  = B(:,d) ./ sqrt(varB + eps);
        Z_field(:,v,d)  = single(z);
        PP_field(:,v,d) = single(0.5 * (1 + erf(z ./ sqrt(2))));
    end
end

%% -----------------------------------------------------------------------
%% 8. GROUP-LEVEL SUMMARY  (expected by your EB pipeline)
%% -----------------------------------------------------------------------
beta_mat = zeros(K1, P-1, 'single');
for d = 1:K1
    for p = 2:P
        beta_mat(d,p-1) = mean(Beta_field(p,:,d), 'omitnan');
    end
end

fprintf('=== Bayesian β estimation completed ===\n');

end
