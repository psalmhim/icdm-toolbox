function Yout = warp_apply_operator(Yin, index8, weight8, dim_out)
% WARP_APPLY_OPERATOR  Apply trilinear warp using precomputed indices.
%
%   Yout = warp_apply_operator(Yin, index8, weight8, dim_out)
%
%   Applies a precomputed 8-neighbour trilinear interpolation warp to
%   flat source data Yin and reshapes the result to dim_out.
%
%   INPUT
%     Yin     : [Nsrc x D] flattened source data
%     index8  : [Nout x 8] uint32 linear indices of 8 neighbours
%     weight8 : [Nout x 8] single trilinear weights
%     dim_out : [1 x 3] output volume dimensions
%
%   OUTPUT
%     Yout : [dim_out x D] warped output volume
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_APPLY_WARP, WARP_APPLY, ICDM_TRILINEAR_INDEX_WEIGHT

D = size(Yin,2);

Yout = zeros(prod(dim_out), D, 'single');

for d = 1:D
    src = Yin(:,d);
    neigh = src(index8);
    val = sum(neigh .* weight8, 2);
    Yout(:,d) = val;
end

Yout = reshape(Yout, [dim_out, D]);
end
