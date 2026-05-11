function pi = ilr_inverse(y_ilr, H)
% ILR_INVERSE  Inverse ILR transform: map ILR coordinates back to the simplex.
%
%   pi = ilr_inverse(y_ilr, H)
%
%   Converts (K-1)-dimensional ILR coordinates to K-part compositional
%   proportions via  pi = softmax(H * y).  This is the inverse of the
%   forward ILR transform  y = H' * clr(pi).
%
%   INPUT
%     y_ilr : [N x (K-1)] matrix of ILR coordinates (e.g. mu or ymap)
%     H     : [K x (K-1)] Helmert submatrix (columns orthonormal), or
%             scalar K (number of components; Helmert matrix built internally)
%
%   OUTPUT
%     pi    : [N x K] matrix of compositional proportions, each row sums to 1
%
%   Example
%     H  = helmert_submatrix(K);
%     pi = ilr_inverse(mu_group, H);
%
%   Reference: Egozcue et al., Math Geol. (2003) 35:279-300.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also HELMERT_SUBMATRIX, ICDM_SUBJECT_VB

if nargin < 2
    error('ilr_inverse requires y_ilr and Helmert matrix H.');
end

[N, D] = size(y_ilr);
if numel(H)==1
    K = H;
    H = helmert_submatrix(K);  % [K × (K-1)]
end

K = size(H,1);
if size(H,2) ~= D
    error('Helmert matrix size mismatch: expect %d columns, got %d.', D, size(H,2));
end

% Compute log-ratio basis projection
Z = y_ilr * H';      % [N × K]

% Exponentiate and normalize (softmax-like)
pi = exp(Z);
pi = bsxfun(@rdivide, pi, sum(pi,2) + eps('single'));

end


function H = helmert_submatrix(K)
% helmert_submatrix  Generate orthonormal Helmert submatrix (K × (K-1))
H = zeros(K, K-1);
for i=1:(K-1)
    H(1:i, i) =  1 / sqrt(i*(i+1));
    H(i+1,i) = -i / sqrt(i*(i+1));
end
end
