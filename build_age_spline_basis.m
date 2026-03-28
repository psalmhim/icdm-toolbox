% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
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

function [B, age_z, knots] = build_age_spline_basis(age, K)
% BUILD_AGE_SPLINE_BASIS  Construct a natural cubic spline basis for age.
%
%   [B, age_z, knots] = build_age_spline_basis(age, K)
%
%   Builds a natural cubic spline basis matrix using K internal knots
%   placed at equally-spaced quantiles of the z-scored age distribution.
%   The basis follows the Hastie & Tibshirani truncated-power formulation
%   and is used in the iCDM population-level regression model to capture
%   nonlinear age effects on compositional connectivity.
%
%   Inputs
%     age   : [S x 1] vector of subject ages.
%     K     : number of internal knots (e.g., 3).
%
%   Outputs
%     B     : [S x (K+1)] spline basis matrix (first column is linear term).
%     age_z : [S x 1] z-scored age vector used internally.
%     knots : [1 x K] knot positions in z-scored age units.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PREPARE_DESIGN_VECTOR, ICDM_DESIGN_VECTOR

    age = age(:);
    S = numel(age);

    % z-score age
    mu = mean(age);
    sd = std(age); if sd == 0, sd = 1; end
    age_z = (age - mu)/sd;

    % choose K internal knots at equally spaced quantiles
    qs = linspace(0,1,K+2); 
    qs = qs(2:end-1);  % drop 0 and 1
    knots = quantile(age_z, qs);

    % build natural cubic spline basis (Hastie & Tibshirani style)
    % columns: [age_z, d1, d2, ... dK], so total K+1 columns
    B = zeros(S, K+1);

    % first basis = linear term
    B(:,1) = age_z;

    % helper for truncated cubic
    function g = hfun(x, k)
        g = max(x - k, 0).^3;
    end

    for j = 1:K
        k_j = knots(j);
        g_j = hfun(age_z, k_j);
        g_K = hfun(age_z, knots(end));
        g_1 = hfun(age_z, knots(1));
        % natural spline transform
        B(:,j+1) = (g_j - g_K)/(knots(end) - k_j) ...
                 - (g_1 - g_K)/(knots(end) - knots(1));
    end
end
