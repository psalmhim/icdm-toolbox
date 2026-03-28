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

function Beta = icdm_estimate_beta(subjects, OUTS, K, opts)
% ICDM_ESTIMATE_BETA  Voxelwise weighted ridge regression for ILR beta fields.
%
%   Beta = icdm_estimate_beta(subjects, OUTS, K, opts)
%
%   Core group-level M-step of the iCDM EB algorithm.  At each MNI
%   analysis voxel, a weighted ridge regression is solved simultaneously
%   for all (K-1) ILR dimensions:
%
%     Beta(:,v,:) = (X'W_v X + lambda*I) \ (X'W_v Y_v)
%
%   Weights W_v combine kappa^alpha (posterior precision of each subject)
%   with coverage^gamma (fraction of subjects covering the voxel).
%   Corresponds to the regression model in the iCDM manuscript
%   (population-level covariate modelling section).
%
%   Inputs
%     subjects : [S x 1] struct array with property fields for design.
%     OUTS     : {S x 1} cell of subject VB outputs, each containing
%                y_ilr_mni [Nmni x (K-1)] and kappa_mni [Nmni x 1].
%     K        : number of compositional targets.
%     opts     : struct with design_spec, idx_mni, ridge_lambda,
%                beta_alpha, beta_gamma.
%
%   Output
%     Beta : [P x Nmni x (K-1)] voxelwise regression coefficients.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_ESTIMATE_BETA_BAYES, ICDM_DESIGN_VECTOR,
%            ICDM_POPULATION_EB
%
% -------------------------------------------------------------------------
fprintf('=== Estimating voxelwise β-field (FAST, weighted ridge) ===\n');

K1   = K - 1;               % ILR dims
Nmni = numel(opts.idx_mni);          % voxels in MNI analysis space

%% 1. BUILD DESIGN MATRIX (X)
X_first = icdm_design_vector(subjects(1), opts);
P       = numel(X_first);
S       = numel(subjects);

X = zeros(S, P, 'single');
X(1,:) = single(X_first);

for s = 2:S
    Xs = icdm_design_vector(subjects(s), opts);
    if numel(Xs) ~= P
        error('Design vector length mismatch at subject %d', s);
    end
    X(s,:) = single(Xs);
end
Xd = double(X);      % design in double precision
clear X Xs

%% 2. STACK SUBJECT ILR (Y) & KAPPA
Yall = zeros(S, Nmni, K1, 'single');
Kap  = zeros(S, Nmni, 'single');

for s = 1:S
    ytmp = OUTS{s}.y_ilr_mni;        % [Nmni × K1]
    ytmp(~isfinite(ytmp)) = 0;

    Yall(s,:,:) = single(ytmp);
    Kap(s,:)    = single(OUTS{s}.kappa_mni);
end
clear OUTS;
%% 3. COVERAGE
coverage = mean(Kap > 0, 1);        % [1 × Nmni]

%% 4. Weighted ridge parameters
alpha  = getfield_default(opts,'beta_alpha',0.5);
gamma  = getfield_default(opts,'beta_gamma',0.2);
lambda = opts.ridge_lambda;

% ridge regularization
Ireg = lambda * eye(P);

%% 5. PRECOMPUTE Kap^alpha AND coverage^gamma

fprintf('  -- precomputing Kap^alpha and coverage^gamma...\n');
Kap_alpha = max(double(Kap), eps('double')).^alpha;     % [S × Nmni]
w_cov_all = max(coverage, eps).^gamma;                  % [1 × Nmni]

clear Kap

%% 6. RESULT ALLOCATION
Beta = zeros(P, Nmni, K1, 'single');

%% 7. MAIN LOOP (per voxel)
fprintf('  -- running weighted ridge for %d voxels...\n', Nmni);
for v = 1:Nmni
    % ------------------------------------------------------
    % Compute weights for voxel v
    % ------------------------------------------------------
    w = Kap_alpha(:,v) * w_cov_all(v);    % [S × 1]
    sw = sqrt(w);                         % [S × 1]

    % Weighted design matrix Xsw = X .* sw
    Xsw = Xd .* sw;                       % [S × P]

    Xt  = Xsw.';                          % [P × S]
    XtX = Xt * Xsw + Ireg;                % [P × P]

    % ------------------------------------------------------
    % Weighted responses for ALL ILR dims simultaneously
    % ------------------------------------------------------
    % Extract Yall(:,v,:) as an [S × K1] matrix
    Yv = double( reshape(Yall(:,v,:), [S K1]) );
    Ysw = bsxfun(@times, Yv, sw);         % [S × K1]

    % ------------------------------------------------------
    % Solve K1 regressions in one shot:
    % B = (XtX) \ (Xt * Ysw)
    % Result: [P × K1]
    % ------------------------------------------------------
    B = XtX \ (Xt * Ysw);

    % store
    Beta(:,v,:) = single(B);
end

fprintf('=== β-field estimation completed ===\n');

end
