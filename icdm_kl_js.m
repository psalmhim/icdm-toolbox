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

function [KL12, KL21, JS] = icdm_kl_js(c1, c2)
% ICDM_KL_JS  Compute KL and Jensen-Shannon divergence between two connectivity profiles.
%
%   [KL12, KL21, JS] = icdm_kl_js(c1, c2)
%
%   Computes the Kullback-Leibler divergence in both directions and the
%   symmetric Jensen-Shannon divergence between two connectivity profiles.
%   Inputs may be raw streamline counts or probability vectors; internal
%   normalization to the probability simplex is applied automatically.
%   A numerical floor (eps = 1e-12) prevents log(0).
%
%   Inputs
%     c1, c2 : [K x 1] or [1 x K] connectivity profiles (counts or
%              probabilities).
%
%   Outputs
%     KL12 : KL(c1 || c2).
%     KL21 : KL(c2 || c1).
%     JS   : Jensen-Shannon divergence (symmetric, bounded).
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_CDM_KL_JS

    % numerical floor
    EPS = 1e-12;

    % ensure column vectors
    c1 = double(c1(:));
    c2 = double(c2(:));

    % check dimension
    if length(c1) ~= length(c2)
        error('c1 and c2 must have the same length');
    end

    % handle empty / zero vectors
    if sum(c1) == 0 || sum(c2) == 0
        KL12 = NaN;
        KL21 = NaN;
        JS   = NaN;
        return;
    end

    % normalize to probability simplex
    p = c1 / sum(c1);
    q = c2 / sum(c2);

    % numerical stabilization
    p = max(p, EPS);
    q = max(q, EPS);

    % renormalize after clipping
    p = p / sum(p);
    q = q / sum(q);

    % midpoint distribution
    m = 0.5 * (p + q);

    % KL divergences
    KL12 = sum(p .* log(p ./ q));
    KL21 = sum(q .* log(q ./ p));

    % Jensen?Shannon divergence
    JS = 0.5 * sum(p .* log(p ./ m)) + ...
         0.5 * sum(q .* log(q ./ m));
end