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

function icdm_ci_recovery(pi_true, ci_low, ci_high, mask, out_png)
% ICDM_CI_RECOVERY  Evaluate confidence interval coverage and generate a recovery map.
%
%   icdm_ci_recovery(pi_true, ci_low, ci_high, mask, out_png)
%
%   Computes the empirical coverage probability, i.e., the fraction of
%   masked voxels where the true compositional probability pi_true falls
%   within the 95 % credible interval [ci_low, ci_high].  A binary
%   recovery map (green = inside CI, red = outside) is saved as PNG.
%
%   This is used to validate the posterior uncertainty estimates from the
%   iCDM variational Bayes inference (see icdm_confidence_interval).
%
%   Inputs
%     pi_true  : ground-truth compositional probability (same size as mask).
%     ci_low   : lower bound of the 95 % credible interval.
%     ci_high  : upper bound of the 95 % credible interval.
%     mask     : logical array selecting valid voxels.
%     out_png  : output file path for the recovery-map PNG.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_CONFIDENCE_INTERVAL, ICDM_GROUP_CI

nV = numel(pi_true);
inside = (pi_true >= ci_low) & (pi_true <= ci_high);
coverage = mean(inside(mask));

fprintf('CI Coverage = %.3f\n', coverage);

% Visualize
heat = zeros(size(mask),'single');
heat(mask) = inside(mask);

figure('Color','w');
imagesc(rot90(heat,1)); axis image off;
colormap([1 0.4 0.4; 0.2 0.8 0.2]); 
title(sprintf('CI Recovery Map (Coverage=%.3f)', coverage));

exportgraphics(gcf, out_png,'Resolution',150);
close(gcf);
end
