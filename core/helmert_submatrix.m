function H = helmert_submatrix(K)
% HELMERT_SUBMATRIX  Orthonormal Helmert contrast matrix for ILR transform.
%
%   H = helmert_submatrix(K)
%
%   Constructs a K x (K-1) matrix H such that H'*H = I and H'*1 = 0,
%   mapping K-part compositions to (K-1)-dimensional ILR coordinates.
%   The ILR forward transform is  y = H' * clr(pi)  and the inverse is
%   pi = softmax(H * y).  (Manuscript Eqs. ilr, softmax)
%
%   INPUT
%     K : number of compositional components (e.g. 68 cortical targets)
%
%   OUTPUT
%     H : [K x (K-1)] orthonormal Helmert submatrix
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ILR_INVERSE, ICDM_SUBJECT_VB

H = zeros(K, K-1);
for i=1:(K-1)
    H(1:i, i) =  1 / sqrt(i*(i+1));
    H(i+1, i) = -i / sqrt(i*(i+1));
end
end

