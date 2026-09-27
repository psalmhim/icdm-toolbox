function tfce_vec = icdm_tfce(stat_vec, idx_mni, dim_mni, E, H, dh, conn)
% ICDM_TFCE  Threshold-free cluster enhancement of a voxelwise statistic.
%
%   tfce_vec = icdm_tfce(stat_vec, idx_mni, dim_mni, E, H, dh, conn)
%
% Standard TFCE (Smith & Nichols, 2009):
%   TFCE(v) = sum_{h=dh:dh:max} extent(v,h)^E * h^H * dh
% where extent(v,h) is the size of the supra-threshold cluster containing v at
% threshold h. Statistic must be non-negative (here T = ||beta_age||^2 studentized).
%
% INPUT
%   stat_vec : [Nmni x 1] statistic on the in-mask voxels idx_mni (>=0)
%   idx_mni  : linear indices of in-mask voxels within a dim_mni volume
%   dim_mni  : [x y z] volume dimensions
%   E,H      : TFCE exponents (defaults 0.5, 2.0 — FSL standard)
%   dh       : FIXED threshold step (must be the SAME for observed and every
%              permutation so the enhancement is one consistent transform)
%   conn     : 3D connectivity (default 26, FSL default)
%
% OUTPUT
%   tfce_vec : [Nmni x 1] TFCE-enhanced values at the same voxels
%
% Notes: cropped to the mask bounding box for speed; parfor-safe (no state).

if nargin<4 || isempty(E),   E=0.5;  end
if nargin<5 || isempty(H),   H=2.0;  end
if nargin<7 || isempty(conn),conn=26; end

stat_vec = double(stat_vec(:));
stat_vec(~isfinite(stat_vec)) = 0;
stat_vec(stat_vec<0) = 0;
tfce_vec = zeros(numel(stat_vec),1);

smax = max(stat_vec);
if smax<=0 || dh<=0, return; end

% place onto 3D grid, crop to bounding box of the mask
V = zeros(dim_mni);
V(idx_mni) = stat_vec;
[ii,jj,kk] = ind2sub(dim_mni, idx_mni(:));
r1=min(ii); r2=max(ii); c1=min(jj); c2=max(jj); s1=min(kk); s2=max(kk);
Vc = V(r1:r2, c1:c2, s1:s2);
Tc = zeros(size(Vc));

hs = dh:dh:smax;
for h = hs
    bw = Vc >= h;
    if ~any(bw(:)), continue; end
    cc = bwconncomp(bw, conn);
    np = cellfun(@numel, cc.PixelIdxList);
    incr = (np.^E) * (h^H) * dh;          % per-cluster increment
    for k = 1:cc.NumObjects
        Tc(cc.PixelIdxList{k}) = Tc(cc.PixelIdxList{k}) + incr(k);
    end
end

Vt = zeros(dim_mni);
Vt(r1:r2, c1:c2, s1:s2) = Tc;
tfce_vec = Vt(idx_mni);
end
