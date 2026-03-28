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

function cmap = coolwarm(n)
% COOLWARM  Generate a diverging blue-white-red colormap.
%
%   cmap = coolwarm(n)
%
%   Returns an [n x 3] colormap that transitions linearly from blue
%   ([59 76 192]/255) through white to red ([180 4 38]/255).  Useful for
%   visualizing signed quantities such as beta coefficients or ILR
%   differences where zero is a natural midpoint.
%
%   Input
%     n : number of colormap entries (default 256).
%
%   Output
%     cmap : [n x 3] RGB colormap matrix with values in [0, 1].
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also PARULA, TURBO, JET

if nargin < 1
    n = 256;
end

cmap = zeros(n,3);

% coolwarm endpoints
c0 = [59,76,192] / 255;   % blue
c1 = [180,4,38] / 255;    % red

% midpoint white
white = [1 1 1];

% interpolate from blue ? white ? red
half = floor(n/2);

for i = 1:half
    t = (i-1)/(half-1);
    cmap(i,:) = (1-t)*c0 + t*white;
end

for i = half+1:n
    t = (i-half)/(n-half);
    cmap(i,:) = (1-t)*white + t*c1;
end
end