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
