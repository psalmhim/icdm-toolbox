function Yout = warp_apply(Ysrc, src_index, warp, dst_index)
% WARP_APPLY  Generic index-based warp with trilinear interpolation.
%
%   Yout = warp_apply(Ysrc, src_index, warp, dst_index)
%
%   Applies a precomputed trilinear warp operator to sparse index-form
%   source data and returns the warped values at the requested destination
%   indices.  Uses chunked processing to limit memory consumption.
%
%   INPUT
%     Ysrc      : [Nsrc x K1] source data (index form)
%     src_index : [Nsrc x 1] linear indices in the source grid
%     warp      : warp struct with .index [Ndst x 8] and .weight [Ndst x 8]
%     dst_index : [Ndst_out x 1] destination voxels to extract
%
%   OUTPUT
%     Yout : [numel(dst_index) x K1] warped data at destination indices
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_APPLY_WARP, ICDM_WARP_TO_MNI, ICDM_WARP_TO_NATIVE

    idx8 = warp.index;      % [Ndst × 8]
    w8   = warp.weight;     % [Ndst × 8]
    Ndst = size(idx8,1);
    [~,K1] = size(Ysrc);

    % Make LUT for sparse mapping
    max_idx = max(idx8(:));
    LUT = zeros(max_idx,1,'uint32');
    LUT(src_index) = uint32(1:size(Ysrc,1));

    Ydst_full = zeros(Ndst,K1,'single');

    chunk=20000;
    nb = ceil(Ndst/chunk);

    for bi=1:nb
        id0 = (bi-1)*chunk+1;
        id1 = min(bi*chunk,Ndst);
        ib  = id0:id1;

        idx_block = idx8(ib,:);
        w_block   = w8(ib,:);

        mapped = LUT(idx_block);

        for k=1:K1
            Yk = Ysrc(:,k);
            V = zeros(numel(ib),8,'single');
            for n=1:8
                rows = mapped(:,n);
                valid = rows>0;
                if any(valid)
                    V(valid,n) = Yk(rows(valid));
                end
            end
            Ydst_full(ib,k) = sum(V.*w_block,2);
        end
    end

    Yout = Ydst_full(dst_index,:);
    Yout(~isfinite(Yout))=0;
end
