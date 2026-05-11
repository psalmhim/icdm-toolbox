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