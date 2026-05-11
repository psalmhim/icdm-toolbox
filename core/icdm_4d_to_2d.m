function R2d=icdm_4d_to_2d(R4d,idx)
% ICDM_4D_TO_2D  Extract index-form rows from a 4-D (or 3-D) volume.
%
%   R2d = icdm_4d_to_2d(R4d, idx)
%
%   Reshapes R4d into [prod(spatial dims) x D] and extracts the rows at
%   linear indices idx, yielding a compact [numel(idx) x D] matrix.
%   This is the standard iCDM operation for converting between full-volume
%   and index-based representations.
%
%   Inputs
%     R4d : array of size [X Y Z D] (or [X Y Z] treated as D=Z).
%     idx : linear indices into the spatial grid.
%
%   Output
%     R2d : [numel(idx) x D] extracted data.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_2D_TO_4D
sz=size(R4d); D=sz(end);
R4d=reshape(R4d,[],D);
R2d=R4d(idx,:);
end
