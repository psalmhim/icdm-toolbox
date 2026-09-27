function [mu_ilr, kappa_group, tau2, beta_pred] = icdm_update_group_prior(OUTS, opts, it)
% ======================================================================
% icdm_update_group_prior.m (robust population prior; optional EB beta^pred)
%
% Group prior update:
%   • Robust symmetric trimmed mean and dispersion
%   • Population-prior precision derived from between-subject dispersion
%   • Optional spatial smoothing (Option B)
%   • Optional EB-shrinkage predictive covariate coefficient beta^pred
%     (4th output; only when opts.beta_pred_design is supplied), estimated from
%     the SAME prior-independent data-only representations (Sec. 2.5.3, general mode).
%
% IMPORTANT:
%   Nv MUST be fixed = numel(opts.idx_mni)
%   OUTS(s).y_data_mni MUST be [Nv × K1]
% ======================================================================

fprintf('[EB] update group prior (iter %03d)\n', it);

S      = numel(OUTS);
idx_mni = opts.idx_mni(:);
Nv      = numel(idx_mni);          % FIXED voxel count
first_valid = find(~cellfun(@isempty,OUTS),1);
if isempty(first_valid), error('No nonempty subject output was supplied.'); end
if ~isfield(OUTS{first_valid},'y_data_mni')
    error('First nonempty subject output lacks y_data_mni.');
end
K1 = size(OUTS{first_valid}.y_data_mni,2);

%% -------------------------------------------------------------
%% Stack prior-independent ILR estimates
%% -------------------------------------------------------------
Ystack = zeros(Nv, K1, S, 'single');

for s = 1:S
    Os = OUTS{s};
    if isempty(Os)
        error('Subject %d OUTS is empty.', s);
    end

    % SAFETY CHECK: size must match Nv
    required = {'y_data_mni'};
    for j = 1:numel(required)
        if ~isfield(Os,required{j})
            error('Subject %d is missing required field %s.',s,required{j});
        end
    end
    if ~isequal(size(Os.y_data_mni),[Nv K1])
        error(['Subject %d y_data_mni dimension mismatch. ' ...
               'Expected [%d %d], got [%d %d]'], s, Nv, K1, ...
               size(Os.y_data_mni,1),size(Os.y_data_mni,2));
    end
    Ystack(:,:,s) = single(Os.y_data_mni);
end

%% -------------------------------------------------------------
%% 1. Trimmed mean & tau2 (robust)
%% -------------------------------------------------------------
alpha       = opts.agg.trim_alpha;
prior_eps2  = opts.precision.prior_eps2;
if ~isscalar(alpha) || alpha < 0 || alpha >= 0.5
    error('opts.agg.trim_alpha must satisfy 0 <= alpha < 0.5.');
end
min_subjects = getfield_default(opts.agg,'min_subjects',3);

mu_ilr = zeros(Nv, K1, 'single');
tau2   = zeros(Nv, K1, 'single');

for v = 1:Nv
    yv = squeeze(Ystack(v,:,:));     % [K1 × S]

    for d = 1:K1

        yd = yv(d,:);
        yd = double(yd(isfinite(yd)));

        if numel(yd) < min_subjects
            mu_ilr(v,d) = NaN;
            tau2(v,d)   = NaN;
            continue;
        end

        yd_sorted = sort(yd);
        L = numel(yd_sorted);
        ntrim = floor(alpha * L);
        yd_trim = yd_sorted((ntrim+1):(L-ntrim));

        % mean
        mu_ilr(v,d) = mean(yd_trim);

        % Robust MAD-based variance
        med_val = median(yd_trim);
        robust_sd = 1.4826 * median(abs(yd_trim - med_val));
        if robust_sd <= 10*eps(max(1,abs(med_val)))
            % MAD can be exactly zero in small/discrete samples.  Fall back
            % to the ordinary variance rather than claiming infinite precision.
            robust_var = var(yd_trim,0);
        else
            robust_var = robust_sd.^2;
        end
        tau2(v,d) = max(robust_var, prior_eps2);
    end
end

%% -------------------------------------------------------------
%% 2. Population prior precision from robust between-subject variance
%% -------------------------------------------------------------
kappa_base  = opts.precision.kappa_base;
kappa_min   = getfield_default(opts.precision,'kappa_min',1e-3);
kappa_max   = getfield_default(opts.precision,'kappa_max',1e4);

kappa_group = zeros(Nv,1,'single');

for v = 1:Nv
    tv = double(tau2(v,:));
    tv = tv(isfinite(tv) & tv > 0);
    if isempty(tv)
        kappa_group(v) = kappa_base;
        continue;
    end
    % The subject solver currently accepts an isotropic ILR prior, so reduce
    % dimension-specific population variances by their arithmetic mean.
    % This is conservative relative to averaging dimension-wise precisions.
    kg = 1 / max(mean(tv),prior_eps2);
    kappa_group(v) = single(min(max(kg,kappa_min),kappa_max));
end

% Voxels with insufficient coverage carry a neutral mean and baseline precision.
bad_mu = ~isfinite(mu_ilr);
mu_ilr(bad_mu) = 0;


%% -------------------------------------------------------------
%% 2b. OPTIONAL: EB predictive covariate coefficient beta^pred (Sec. 2.5.3).
%%     Estimated HERE, at the group level, by EB-shrinkage ridge from the SAME
%%     prior-independent data-only representations (deviations about mu_ilr).
%%     beta^pred is the deviation-form predictive coefficient used only in the
%%     covariate-informed (individualized) prior mean; it is NOT the downstream
%%     scientific coefficient beta^assoc, and is not used in covariate-agnostic tests.
%% -------------------------------------------------------------
beta_pred = [];
if nargout >= 4 && isfield(opts,'beta_pred_design') && ~isempty(opts.beta_pred_design)
    Xd = double(opts.beta_pred_design);           % [S x P] covariates (centered internally)
    Xd = Xd - mean(Xd,1);
    P  = size(Xd,2);
    beta_pred = zeros(Nv, K1, P, 'single');
    if P == 1
        x = Xd(:); sx = x' * x;                    % scalar X'X (age-only case)
        for d = 1:K1
            Rd = double(squeeze(Ystack(:,d,:)));   % [Nv x S] data-only reps
            Rd = Rd - double(mu_ilr(:,d));         % deviations about group mean
            Rd(~isfinite(Rd)) = 0;
            Xtr = Rd * x;  RtR = sum(Rd.^2, 2);    % X'r and ||r||^2 per voxel
            lam = ones(Nv,1);                      % vectorized voxelwise MacKay REML
            for it2 = 1:30
                b    = Xtr ./ (sx + lam);
                geff = sx  ./ (sx + lam);
                RSS  = max(RtR - 2*b.*Xtr + b.^2.*sx, 0);
                tau2b = max(b.^2 ./ max(geff,1e-9),   1e-12);
                s2b   = max(RSS  ./ max(S - geff,1e-9), 1e-12);
                lam   = s2b ./ tau2b;
            end
            beta_pred(:,d,1) = single(Xtr ./ (sx + lam));
        end
    else
        XtX = Xd' * Xd;                            % general P: per-voxel ridge REML
        for v = 1:Nv
            for d = 1:K1
                r = double(squeeze(Ystack(v,d,:))) - double(mu_ilr(v,d));
                r(~isfinite(r)) = 0;
                beta_pred(v,d,:) = single(reml_ridge_local(Xd, XtX, r, S, P));
            end
        end
    end
    fprintf('[EB] beta^pred (EB-shrinkage ridge, P=%d) estimated from data-only reps.\n', P);
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
    % Normalized (scale-preserving) reliability smoothing: a convex combination
    % of the voxel's own precision and its neighbors', so smoothing introduces
    % spatial coherence WITHOUT inflating the overall precision scale (an additive
    % update K + lambda*deg would spuriously grow precision with neighbor count).
    K2(v)    = (K(v) + lambda * sum(K(nv))) ./ (1 + lambda * deg);
end
end


function b = reml_ridge_local(X, XtX, r, S, P)
% EB-shrinkage ridge (MacKay evidence fixed point) for one voxel/dimension,
% y = X b + e, b~N(0,tau2 I), e~N(0,s2 I). Penalty s2/tau2 from marginal likelihood.
Xtr = X' * r;  s2 = var(r); if ~isfinite(s2)||s2<=0, s2 = 1e-6; end
tau2 = s2; b = zeros(P,1);
for it = 1:50
    lam = s2 / max(tau2,1e-12);
    A   = XtX + lam*eye(P);
    b   = A \ Xtr;
    geff = min(max(P - lam*trace(inv(A)), 1e-6), P);
    tau2 = max((b'*b)/geff, 1e-12);
    resid = r - X*b;
    s2   = max((resid'*resid)/max(S-geff,1e-6), 1e-12);
end
b = (XtX + (s2/max(tau2,1e-12))*eye(P)) \ Xtr;
end
