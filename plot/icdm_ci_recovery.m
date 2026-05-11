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
