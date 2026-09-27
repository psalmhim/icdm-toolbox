function Yout = icdm_apply_warp(Yin, warp)
% Yin  : [d1 d2 d3 D]
% warp : trilinear operator

dim_in  = warp.dim_in;
dim_out = warp.dim_out;

if ~isequal(size(Yin,1:3), dim_in)
    error('icdm_apply_warp: input dim mismatch.');
end

D = size(Yin,4);

Yin_flat = reshape(Yin,[],D);  % [N_in × D]
index8  = warp.index;
weight8 = warp.weight;

Yout = zeros([dim_out D],'single');
for d = 1:D
    src = Yin_flat(:,d);
    neigh = src(index8);        % [Nout × 8]
    dst = sum(neigh .* weight8,2);
    data=reshape(dst, dim_out);
    Yout(:,:,:,d) = data;
end
end
