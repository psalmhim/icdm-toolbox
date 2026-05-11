function R4d=icdm_2d_to_4d(R2d,idx,dim)
% ICDM_2D_TO_4D  Scatter index-form data back into a full 3-D or 4-D volume.
%
%   R4d = icdm_2d_to_4d(R2d, idx, dim)
%
%   Inverse of icdm_4d_to_2d.  Places the rows of R2d into the voxel
%   locations given by linear indices idx, producing a volume of size dim.
%   When dim has 3 elements, R2d is treated as a single-column vector and
%   the output is 3-D; when dim has 4 elements, R2d is [numel(idx) x D]
%   and the output is 4-D.
%
%   Inputs
%     R2d : [Nidx x D] or [Nidx x 1] data in index form.
%     idx : linear voxel indices into the output grid.
%     dim : [X Y Z] for 3-D output, or [X Y Z D] for 4-D output.
%
%   Output
%     R4d : single array of size dim with zeros at non-indexed locations.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_4D_TO_2D
if length(dim)==3
    R4d=zeros(dim,'single');
    R4d(idx)=R2d;
else
    D=dim(end);
    R4d=zeros(prod(dim(1:end-1)),D,'single');
    R4d(idx,:)=R2d;
    R4d=reshape(R4d,dim);
end
end
