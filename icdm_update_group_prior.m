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

function [mu_ilr, kappa_group, tau2] = icdm_update_group_prior(OUTS, opts, it)
% ICDM_UPDATE_GROUP_PRIOR  M-step: robust group aggregation of ILR posteriors.
%
%   [mu_ilr, kappa_group, tau2] = icdm_update_group_prior(OUTS, opts, it)
%
%   Implements the group-level parameter update (M-step) of the iterative
%   empirical Bayes loop (manuscript Section: Empirical Bayes group priors).
%
%   For each voxel v and ILR dimension d:
%     1. Trimmed mean of subject MAP estimates (Eq. group_mu)
%     2. MAD-based robust dispersion tau2 (Eq. group_dispersion)
%     3. Weighted kappa aggregation (Eq. group-reliability)
%     4. Optional spatial smoothing of mu and kappa (Eq. spatial_smoothing)
%
%   INPUT
%     OUTS : {S x 1} cell array of subject posterior structs (from icdm_subject_vb)
%     opts : options struct (.idx_mni, .agg.trim_alpha, .precision, .spatial_group)
%     it   : current EB iteration number
%
%   OUTPUT
%     mu_ilr      : [Nv x (K-1)] group prior mean in ILR space
%     kappa_group : [Nv x 1] group-level posterior reliability
%     tau2        : [Nv x (K-1)] per-dimension dispersion (MAD^2)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_POPULATION_EB, ICDM_SUBJECT_VB, SPATIAL_SMOOTH_PRIOR

fprintf('[EB] update group prior (iter %03d)\n', it);

S      = numel(OUTS);
idx_mni = opts.idx_mni(:);
Nv      = numel(idx_mni);          % FIXED voxel count
K1      = size(OUTS{find(~cellfun(@isempty,OUTS),1)}.y_ilr_mni, 2);

%% -------------------------------------------------------------
%% Stack ILR and kappa
%% -------------------------------------------------------------
Ystack = zeros(Nv, K1, S, 'single');
Kstack = zeros(Nv, S,    'single');

for s = 1:S
    Os = OUTS{s};
    if isempty(Os)
        error('Subject %d OUTS is empty.', s);
    end

    % SAFETY CHECK: size must match Nv
    if size(Os.y_ilr_mni,1) ~= Nv
        error(['Subject %d y_ilr_mni dimension mismatch. ' ...
               'Expected Nv=%d, got %d'], s, Nv, size(Os.y_ilr_mni,1));
    end

    Ystack(:,:,s) = single(Os.y_ilr_mni);
    Kstack(:,s)   = single(Os.kappa_mni);
end

%% -------------------------------------------------------------
%% 1. Trimmed mean & tau2 (robust)
%% -------------------------------------------------------------
alpha       = opts.agg.trim_alpha;
prior_eps2  = opts.precision.prior_eps2;

mu_ilr = zeros(Nv, K1, 'single');
tau2   = zeros(Nv, K1, 'single');

for v = 1:Nv
    yv = squeeze(Ystack(v,:,:));     % [K1 × S]

    for d = 1:K1

        yd = yv(d,:);
        yd = yd(~isnan(yd));

        if numel(yd) < 3
            mu_ilr(v,d) = 0;
            tau2(v,d)   = prior_eps2;
            continue;
        end

        yd_sorted = sort(yd);
        L = numel(yd_sorted);
        kL = max(1, round(alpha * L));
        kH = min(L, round((1 - alpha) * L));
        yd_trim = yd_sorted(kL:kH);

        % mean
        mu_ilr(v,d) = mean(yd_trim);

        % MAD variance
        mad_val = median(abs(yd_trim - mu_ilr(v,d)));
        tau2(v,d) = max(mad_val.^2, prior_eps2);
    end
end

%% -------------------------------------------------------------
%% 2. Weighted Aggregate Kappa
%% -------------------------------------------------------------
alpha_kappa = getfield_default(opts, 'kappa_alpha', 0.5);
kappa_base  = opts.precision.kappa_base;

kappa_group = zeros(Nv,1,'single');

for v = 1:Nv
    kv = Kstack(v,:);
    kv = kv(~isnan(kv));

    if isempty(kv)
        kappa_group(v) = kappa_base;
        continue;
    end

    % weights = κ^α
    w = kv .^ alpha_kappa;

    kg = sum(w .* kv) / sum(w);

    if ~isfinite(kg) || kg <= 0
        kg = kappa_base;
    end

    kappa_group(v) = single(kg);
end


%% -------------------------------------------------------------
%% 3. OPTIONAL: Group-level Spatial Smoothing (Option B)
%% -------------------------------------------------------------
if isfield(opts,'spatial_group') && opts.spatial_group.use

    fprintf('[EB] Group-level spatial smoothing (Option B)...\n');

    dim = opts.dim_mni;     % [X Y Z]

    mask3d = false(dim);
    mask3d(idx_mni) = true;

    idx_mask = find(mask3d(:));   % full mask voxel ids
    numMask  = numel(idx_mask);

    % LUT: full-index → mask-index
    LUT = zeros(prod(dim),1,'uint32');
    LUT(idx_mask) = 1:numMask;

    % Expand MU and K into full mask-space
    MU_mask = zeros(numMask, K1, 'single');
    K_mask  = zeros(numMask, 1,   'single');

    % SAFETY: LUT(idx_mni) must equal 1:Nv
    MU_mask( LUT(idx_mni), : ) = mu_ilr;
    K_mask( LUT(idx_mni) )    = kappa_group;

    % neighbor list
    nbr = get_neighbors_mask(mask3d, dim, LUT, opts.spatial_group.neighborhood);

    % smooth
    [MU_mask_s, K_mask_s] = spatial_smooth_prior_mask( ...
        MU_mask, K_mask, nbr, opts.spatial_group.lambda);

    % back to WM-index order (Nv)
    mu_ilr      = MU_mask_s(LUT(idx_mni), :);
    kappa_group = K_mask_s(LUT(idx_mni));

else
    fprintf('[EB] Spatial smoothing OFF.\n');
end

fprintf('[EB] group prior updated (iter %03d)\n', it);
end


function nbr = get_neighbors_mask(mask3d, dim, LUT, neighborhood)
% Build neighbor lists for mask-indexed voxels
% mask3d : logical [X Y Z]
% dim    : [X Y Z]
% LUT    : full-index → mask-index (0 = not in mask)
% neighborhood = 6 or 26

idx_mask = find(mask3d(:));
numMask  = numel(idx_mask);

nbr = cell(numMask,1);

X = dim(1); Y = dim(2); Z = dim(3);

if neighborhood == 26
    offsets = [];
    cnt = 1;
    for dx=-1:1
        for dy=-1:1
            for dz=-1:1
                if dx==0 && dy==0 && dz==0, continue; end
                offsets(cnt,:) = [dx dy dz]; %#ok<AGROW>
                cnt = cnt + 1;
            end
        end
    end
else
    offsets = [
        1 0 0
       -1 0 0
        0 1 0
        0 -1 0
        0 0 1
        0 0 -1
    ];
end

for k = 1:numMask
    lin = idx_mask(k);
    [x,y,z] = ind2sub(dim, lin);

    neigh = [];

    for j = 1:size(offsets,1)
        xn = x + offsets(j,1);
        yn = y + offsets(j,2);
        zn = z + offsets(j,3);

        if xn>=1 && xn<=X && yn>=1 && yn<=Y && zn>=1 && zn<=Z
            lin2 = sub2ind(dim, xn,yn,zn);
            mk = LUT(lin2);
            if mk > 0
                neigh(end+1) = mk; %#ok<AGROW>
            end
        end
    end

    nbr{k} = neigh;
end

end


function [MU2, K2] = spatial_smooth_prior_mask(MU, K, nbr, lambda)
% One-pass spatial smoothing of group-level priors
% MU, K defined on mask-index space
% nbr is {v} neighbor index list
% lambda is smoothing strength

if lambda <= 0
    MU2 = MU;
    K2  = K;
    return;
end

[numV, K1] = size(MU);
MU2 = MU;
K2  = K;

for v = 1:numV
    nv = nbr{v};
    if isempty(nv), continue; end

    deg = numel(nv);

    num = K(v)*MU(v,:) + lambda * sum(MU(nv,:),1);
    den = K(v) + lambda * deg;

    MU2(v,:) = num ./ den;
    K2(v)    = den;
end
end

