function [Beta, stats, info] = icdm_estimate_beta(design_input, group, K, opts)
% ICDM_ESTIMATE_BETA  (FULL IMPLEMENTATION / SINGLE FILE)
% ========================================================================
% Voxelwise regression in ILR space (K1=K-1) + inference options.
%
% KEY DESIGN GOAL:
%   - Correct + publication-ready defaults
%   - Memory-safe permutation mode by default (avoids MATLAB crashes)
%   - Optional "fast" permutation mode (uses more RAM; faster)
%   - No nested functions inside parfor (parfor-safe)
%
% ------------------------------------------------------------------------
% opts.estimator  (default 'wls-ridge'):
%   'ols'        : ordinary least squares
%   'wls'        : weighted least squares (weights from kappa + coverage)
%   'wls-ridge'  : wls + ridge penalty
%
% opts.inference (default 'signflip'):
%   'wald'         : analytic chi-square Wald (often liberal if var underestimated)
%   'signflip'     : permutation sign-flip on FULL-model residuals (max-T FWER)
%   'freedman-lane': permutation Freedman–Lane (permute residuals from REDUCED model)
%
% PERMUTATION MODES:
%   opts.perm_mode = 'safe' (default)  : recompute per-voxel residuals on-the-fly (low RAM, slower)
%                   'fast'            : precompute residual tensors (high RAM, faster; may crash)
%
% MULTIVARIATE STATISTIC:
%   T(v) = || beta_age(v,:) ||_2^2   in ILR space (K1 dims)
%
% MULTIPLE COMPARISON OUTPUTS (both reported):
%   - stats.p_map      : voxelwise p-values (permutation => max-T FWER; Wald => uncorrected)
%   - stats.q_map      : BH-FDR q-values (reporting; not required if using max-T)
%   - stats.sig_mask   : BH significant mask
%   - stats.p_crit     : BH critical p threshold
%
% CLUSTER (optional):
%   opts.do_cluster (default false)
%   Cluster correction here is **cluster-size FWER** using permutation null
%   of maximum cluster size (computed during permutation).
%
% Required group fields:
%   group{s}.y_ilr_mni : [Nmni x (K-1)]
%   group{s}.kappa_mni : [Nmni x 1]
%
% Required opts fields:
%   opts.idx_mni       : linear indices into a 3D MNI volume (Nmni entries)
%   opts.dim_mni       : [X Y Z] of MNI volume
%   opts.design_spec   : cellstr like {'const','age', ...}
%   opts.design.mu.age, opts.design.sd.age  (if 'age' used)
%
% ========================================================================

fprintf('\n=== ICDM: Estimation + Inference ===\n');

% ----------------------------- defaults
opts.estimator        = getfield_default(opts,'estimator','wls-ridge');
opts.inference        = getfield_default(opts,'inference','signflip');
opts.ridge_lambda     = getfield_default(opts,'ridge_lambda',1e-3);
opts.beta_alpha       = getfield_default(opts,'beta_alpha',0.5);
opts.beta_gamma       = getfield_default(opts,'beta_gamma',0.2);
opts.n_perm           = getfield_default(opts,'n_perm',2000);
opts.batch_size       = getfield_default(opts,'batch_size',256);
opts.perm_seed        = getfield_default(opts,'perm_seed',[]);
opts.perm_mode        = getfield_default(opts,'perm_mode','safe'); % 'safe'|'fast'
opts.usePar           = getfield_default(opts,'usePar',[]);
opts.fdr_q            = getfield_default(opts,'fdr_q',0.05);

opts.do_cluster        = getfield_default(opts,'do_cluster',false);
opts.cluster_alpha     = getfield_default(opts,'cluster_alpha',0.05); % cluster-FWER
opts.cluster_form_q    = getfield_default(opts,'cluster_form_q',0.99); % voxelwise forming threshold quantile on null-maxT
opts.cluster_min_size  = getfield_default(opts,'cluster_min_size',20);
opts.use_pca         = getfield_default(opts,'use_pca',false);
opts.pca_mode        = getfield_default(opts,'pca_mode','fixed');   % 'fixed' or 'auto'
opts.pca_ncomp       = getfield_default(opts,'pca_ncomp',10);
opts.pca_var_ratio   = getfield_default(opts,'pca_var_ratio',0.90);
opts.pca_max_rank    = getfield_default(opts,'pca_max_rank',20);
opts.pca_min_rank    = getfield_default(opts,'pca_min_rank',3);
opts.pca_center      = getfield_default(opts,'pca_center',false);   % ? ??? ??
% sanity
if ~isfield(opts,'idx_mni') || ~isfield(opts,'dim_mni')
    error('opts.idx_mni and opts.dim_mni are required for cluster mapping.');
end

% ----------------------------- 1) design
[X,S,P] = icdm_build_design(design_input, opts);
age_idx = find(strcmp(opts.design_spec,'age'));
if isempty(age_idx)
    error('Age predictor not found in opts.design_spec.');
end

% ----------------------------- 2) stack data
[Yall, Kap] = icdm_stack_data(group, S);
coverage = mean(Kap > 0, 1);
[S2,Nmni_full,K1] = size(Yall);
if S2 ~= S, error('Internal size mismatch.'); end
if K1 ~= (K-1)
    warning('K1 from data is %d but K-1 is %d. Using K1 from data.', K1, K-1);
end

% ----------------------------- 2b) kappa threshold mask (default 50)
grp_kappa_thresh = getfield_default(opts, 'grp_kappa_thresh', 50);
idx_mni_full = opts.idx_mni;   % save full idx_mni before any masking
wm_mask = [];
if grp_kappa_thresh > 0 && isfield(opts,'grp_kappa_mni') && ~isempty(opts.grp_kappa_mni)
    kappa_grp = opts.grp_kappa_mni(:);
    wm_mask   = kappa_grp > grp_kappa_thresh;
    if sum(wm_mask) < Nmni_full
        Yall     = Yall(:, wm_mask, :);
        Kap      = Kap(:,  wm_mask);
        coverage = coverage(wm_mask);
        % update idx_mni in opts so cluster correction maps correctly
        opts.idx_mni = opts.idx_mni(wm_mask);
        fprintf('  -- kappa>%.0f mask: %d / %d voxels (%.1f%%)\n', ...
            grp_kappa_thresh, sum(wm_mask), Nmni_full, 100*sum(wm_mask)/Nmni_full);
    end
end
[S2, Nmni, K1] = size(Yall);

fprintf('  S=%d subjects | P=%d predictors | Nmni=%d voxels | K1=%d ILR dims\n', S,P,Nmni,K1);
fprintf('  estimator=%s | inference=%s | perm_mode=%s\n', opts.estimator, opts.inference, opts.perm_mode);

% ======================= PCA BLOCK (AUTO RANK) =======================
if opts.use_pca
    fprintf('  -- Applying GLOBAL PCA\n');
    % ---- 1) Global covariance 계산
    C = zeros(K1,K1);
    for s = 1:S
        Ys = double(squeeze(Yall(s,:,:)));

        if opts.pca_center
            Ys = Ys - mean(Ys,1);
        end
        C = C + (Ys' * Ys);
    end
    C = C / (S * Nmni);
    % ---- 2) Eigen decomposition
    [V,D] = eig(C);
    evals = diag(D);
    [evals,idx] = sort(evals,'descend');
    V = V(:,idx);

    total_var = sum(evals);
    cumvar = cumsum(evals) / total_var;

    % ---- 3) Rank selection
    if strcmpi(opts.pca_mode,'auto')
        r = find(cumvar >= opts.pca_var_ratio,1);
        if isempty(r)
            r = length(evals);
        end
        % 안전장치
        r = max(r, opts.pca_min_rank);
        r = min(r, opts.pca_max_rank);
        r = min(r, S-5);   % S 대비 과도한 rank 방지
        fprintf('     auto-selected rank = %d (%.2f%% variance explained)\n', ...
                r, 100*cumvar(r));
    else
        r = min(opts.pca_ncomp, K1);
        fprintf('     fixed rank = %d\n', r);
    end

    % ---- 4) Projection
    V = V(:,1:r); nvtx=size(Yall,2);
    % Optional: save the PCA ILR basis (K1_orig x r) so the age effect can be
    % mapped back to composition/target space offline (H * V * beta_pca).
    if isfield(opts,'save_pca_basis') && ~isempty(opts.save_pca_basis)
        try, save(opts.save_pca_basis, 'V', '-v7'); fprintf('  -- saved PCA basis -> %s\n', opts.save_pca_basis); catch, end
    end
    Yall1=zeros(S,nvtx,r,'single');
    for s = 1:S
        Ys = double(squeeze(Yall(s,:,:)));
        Yall1(s,:,:) = single(Ys * V);
    end
    Yall=Yall1;Yall1=[];
    K1 = r;
    fprintf('  -- PCA reduced ILR dims to %d\n', K1);
end
% =====================================================================

% ----------------------------- 3) estimation
Beta = icdm_beta_estimator(X, Yall, Kap, coverage, opts);

% ----------------------------- 4) inference
switch lower(opts.inference)
    case 'none'
        % point estimate only — skip all testing (intermediate EB iterations)
        stats = struct('p_map',[],'T_obs',[],'T_null_max',[],...
                       'q_map',[],'sig_mask',[],'p_crit',[]);
        fprintf('  -- Inference skipped (point estimate only)\n');
        % expand Beta back and return immediately
        if ~isempty(wm_mask) && sum(wm_mask) < Nmni_full
            Beta_full = zeros(P, Nmni_full, K1, 'single');
            Beta_full(:, wm_mask, :) = Beta;
            Beta = Beta_full;
        end
        info.estimator = opts.estimator; info.inference = 'none';
        info.S = S; info.P = P; info.Nmni = Nmni_full; info.K1 = K1;
        fprintf('=== Done ===\n\n');
        return;

    case 'wald'
        stats = icdm_wald_test(X, Yall, Kap, coverage, Beta, age_idx, opts);

    case 'signflip'
        stats = icdm_perm_signflip_wls(X, Yall, Kap, coverage, Beta, age_idx, opts);

    case {'freedman-lane','freedmanlane'}
        stats = icdm_perm_freedman_lane_wls(X, Yall, Kap, coverage, Beta, age_idx, opts);

    otherwise
        error('Unknown inference: %s', opts.inference);
end

% ----------------------------- 5) BH-FDR (reporting)
if isfield(stats,'p_uncorr') && all(isfinite(stats.p_uncorr))
    [stats.q_map, stats.sig_mask, stats.p_crit] = bh_fdr(stats.p_uncorr, opts.fdr_q);
else
    [stats.q_map, stats.sig_mask, stats.p_crit] = bh_fdr(stats.p_map, opts.fdr_q);
end

% ----------------------------- 6) cluster correction (optional)
if opts.do_cluster
    if ~isfield(stats,'max_cluster_null') || isempty(stats.max_cluster_null)
        warning('Cluster requested but max_cluster_null is missing. Cluster skipped.');
    else
        stats.cluster_map = icdm_cluster_correction(stats, opts);
    end
end

% ----------------------------- 7) expand Beta + stats back to Nmni_full
if ~isempty(wm_mask) && sum(wm_mask) < Nmni_full
    Beta_full = zeros(P, Nmni_full, K1, 'single');
    Beta_full(:, wm_mask, :) = Beta;
    Beta = Beta_full;

    stat_fields = {'p_map','p_uncorr','T_obs','q_map','sig_mask','tfce_obs','tfce_p_fwer'};
    for fi = 1:numel(stat_fields)
        fn = stat_fields{fi};
        if isfield(stats,fn) && ~isempty(stats.(fn))
            v_full = nan(Nmni_full, 1);
            v_full(wm_mask) = stats.(fn)(:);
            stats.(fn) = v_full;
        end
    end
    % restore full idx_mni so downstream code (saving, NIfTI write) works correctly
    opts.idx_mni = idx_mni_full;
    fprintf('  -- Expanded results back to Nmni_full=%d\n', Nmni_full);
end

% ----------------------------- info
info = struct();
info.estimator        = opts.estimator;
info.inference        = opts.inference;
info.perm_mode        = opts.perm_mode;
info.grp_kappa_thresh = grp_kappa_thresh;
info.S = S; info.P = P; info.Nmni = Nmni_full; info.K1 = K1;

fprintf('=== Done ===\n\n');
end

% ========================================================================
%  ESTIMATOR
% ========================================================================
function Beta = icdm_beta_estimator(X, Yall, Kap, coverage, opts)

[S,Nmni,K1] = size(Yall);
P = size(X,2);

alpha  = opts.beta_alpha;
gamma  = opts.beta_gamma;
lambda = opts.ridge_lambda;

Kap_alpha = max(Kap, eps).^alpha;        % [S×Nmni]
w_cov_all = max(coverage, eps).^gamma;   % [1×Nmni]

Beta = zeros(P,Nmni,K1,'single');

XtX_ols = X' * X;   % reuse

for v = 1:Nmni
    Yv = double(reshape(Yall(:,v,:), [S K1]));  % [S×K1]

    switch lower(opts.estimator)
        case 'ols'
            B = XtX_ols \ (X' * Yv);

        case 'wls'
            w  = double(Kap_alpha(:,v) * w_cov_all(v));
            sw = sqrt(w);
            Xsw = X .* sw;
            Ysw = Yv .* sw;
            B = (Xsw' * Xsw) \ (Xsw' * Ysw);

        case {'wls-ridge','wlsridge'}
            w  = double(Kap_alpha(:,v) * w_cov_all(v));
            sw = sqrt(w);
            Xsw = X .* sw;
            Ysw = Yv .* sw;
            XtX = (Xsw' * Xsw) + lambda * eye(P);
            B = XtX \ (Xsw' * Ysw);

        otherwise
            error('Unknown estimator: %s', opts.estimator);
    end

    Beta(:,v,:) = single(B);
end
end

% ========================================================================
%  WALD TEST (analytic)  -- often liberal if sigma2 underestimated
% ========================================================================
function stats = icdm_wald_test(X, Yall, Kap, coverage, Beta, age_idx, opts)

fprintf('  -- Wald multivariate test\n');

[S,Nmni,K1] = size(Yall);
P = size(X,2);

alpha  = opts.beta_alpha;
gamma  = opts.beta_gamma;
lambda = opts.ridge_lambda;

Kap_alpha = max(Kap, eps).^alpha;
w_cov_all = max(coverage, eps).^gamma;

dof = max(S - P, 1);

p_map = zeros(Nmni,1);
T_obs = zeros(Nmni,1);

XtX_ols = X' * X;

for v = 1:Nmni

    Yv = double(reshape(Yall(:,v,:), [S K1]));    % [S×K1]
    B  = double(reshape(Beta(:,v,:), [P K1]));    % [P×K1]

    if strcmpi(opts.estimator,'ols')
        w = ones(S,1);
        XtX_inv = XtX_ols \ eye(P);
    else
        w  = double(Kap_alpha(:,v) * w_cov_all(v));
        sw = sqrt(w);
        Xsw = X .* sw;
        XtX = (Xsw' * Xsw) + lambda * eye(P);
        XtX_inv = XtX \ eye(P);
    end

    R = Yv - X * B;                              % [S×K1]
    WRSS = sum((R.^2) .* w, 1);                  % [1×K1]
    sigma2 = WRSS ./ dof;                        % [1×K1]

    var_age  = sigma2(:) .* XtX_inv(age_idx,age_idx);  % [K1×1]
    beta_age = B(age_idx,:).';                         % [K1×1]

    Z2 = sum(beta_age.^2 ./ var_age);            % chi2 approx
    T_obs(v) = Z2;
    p_map(v) = 1 - chi2cdf(Z2, K1);
end

stats = struct();
stats.p_map = p_map;
stats.T_obs = T_obs;
end

% ========================================================================
%  SIGN-FLIP PERMUTATION (WLS/WLS+RIDGE CONSISTENT)
%  - max-T FWER p-values
%  - cluster-size FWER optional null (max cluster size per perm)
%
%  IMPORTANT:
%   - For WLS/WLS+ridge, we must refit per voxel for each perm.
%   - Default perm_mode='safe' avoids storing gigantic residual tensors.
% ========================================================================
function stats = icdm_perm_signflip_wls(X,Yall,Kap,coverage,Beta,age_idx,opts)

fprintf('  -- Sign-flip permutation (WLS-consistent) | n_perm=%d | mode=%s\n', opts.n_perm, opts.perm_mode);

if ~isempty(opts.perm_seed)
    rng(opts.perm_seed);
end

[S,Nmni,K1] = size(Yall);
P = size(X,2);

alpha  = opts.beta_alpha;
gamma  = opts.beta_gamma;
lambda = opts.ridge_lambda;
nperm  = opts.n_perm;
Bsz    = opts.batch_size;

Kap_alpha = max(Kap, eps).^alpha;
w_cov_all = max(coverage, eps).^gamma;

dof = max(S - P, 1);

% -----------------------------
% 1) Observed statistic (studentized)
% -----------------------------
[T_obs] = icdm_compute_Tobs_studentized(X,Yall,Kap_alpha,w_cov_all,Beta,age_idx,lambda,dof);

% -----------------------------
% 2) permutation loop
%    - max-T null (FWER)
%    - voxelwise null counts (for uncorrected p -> BH-FDR)
% -----------------------------
T_null_max = zeros(nperm,1,'single');
exceed_cnt = zeros(Nmni,1,'uint32');   % voxelwise exceedance counts (uncorrected p)
doClusterNull = isfield(opts,'do_cluster') && opts.do_cluster;
max_cluster_null = zeros(nperm,1,'single');

% NOTE: voxelwise counting + cluster null은 reduction이 필요해서
%       parfor로 깔끔하게 하려면 더 복잡해짐.
%       결과를 "정확히" 얻는 목적이면 여기서는 FOR 권장.
usePar = false; % 안전하게 꺼두는 것을 권장 (메모리/리덕션 이슈)
if isfield(opts,'usePar') && ~isempty(opts.usePar)
    usePar = opts.usePar;
end

if usePar
    warning('Signflip: voxelwise permutation p(FDR)까지 하려면 for-loop가 가장 안전합니다. 현재는 max-T만 병렬로 수행합니다.');
    parfor m = 1:nperm
        [Tmax_m, ~, ~] = icdm_signflip_oneperm_stats( ...
            X, Yall, Kap_alpha, w_cov_all, Beta, age_idx, lambda, dof, ...
            Bsz, false, opts.dim_mni, opts.idx_mni, opts.cluster_form_q);
        T_null_max(m) = single(Tmax_m);
    end
else
    for m = 1:nperm
        if rem(m,50)==0, fprintf('    perm %d/%d\n', m, nperm); end

        sgn = sign(randn(S,1));
        sgn(sgn==0) = 1;

        [Tmax_m, Tperm_all, maxCl_m] = icdm_signflip_oneperm_stats( ...
            X, Yall, Kap_alpha, w_cov_all, Beta, age_idx, lambda, dof, ...
            Bsz, doClusterNull, opts.dim_mni, opts.idx_mni, opts.cluster_form_q, sgn);

        T_null_max(m) = single(Tmax_m);

        % voxelwise exceed count (for permutation p_uncorr)
        exceed_cnt = exceed_cnt + uint32(Tperm_all >= T_obs);

        if doClusterNull
            max_cluster_null(m) = single(maxCl_m);
        end
    end
end

% -----------------------------
% 3) p-values
% -----------------------------
% (A) max-T FWER p-values
p_fwer = zeros(Nmni,1);
for v = 1:Nmni
    p_fwer(v) = (1 + sum(T_null_max >= T_obs(v))) / (nperm + 1);
end

% (B) voxelwise permutation p-values (uncorrected) for BH-FDR
if ~usePar
    p_uncorr = (double(exceed_cnt) + 1) / (nperm + 1);
else
    p_uncorr = nan(Nmni,1); % not computed in parfor branch
end

stats = struct();
stats.T_obs = T_obs;
stats.T_null_max = T_null_max;

stats.p_map = p_fwer;       % 기존 인터페이스 유지: p_map = FWER max-T
stats.p_fwer = p_fwer;
stats.p_uncorr = p_uncorr;  % BH-FDR용

if doClusterNull && ~usePar
    stats.max_cluster_null = max_cluster_null;
else
    stats.max_cluster_null = [];
end

fprintf('  -- Sign-flip done\n');
end


function [Tmax, Tperm_all, maxCluster] = icdm_signflip_oneperm_stats( ...
    X, Yall, Kap_alpha, w_cov_all, Beta, age_idx, lambda, dof, ...
    batch_size, doClusterNull, dim_mni, idx_mni, cluster_form_q, sgn)

[S,Nmni,K1] = size(Yall);
P = size(X,2);

if nargin < 14 || isempty(sgn)
    sgn = sign(randn(S,1));
    sgn(sgn==0) = 1;
end

Tmax = 0;
maxCluster = 0;

Tperm_all = zeros(Nmni,1,'single');  % needed for voxelwise p + optional cluster

nbatch = ceil(Nmni / batch_size);

for b = 1:nbatch
    v_start = (b-1)*batch_size + 1;
    v_end   = min(b*batch_size, Nmni);

    for v = v_start:v_end

        % voxel weights
        w  = double(Kap_alpha(:,v) * w_cov_all(v));
        sw = sqrt(w);

        % observed fit (from Beta) and residual
        Bv   = double(reshape(Beta(:,v,:), [P K1]));   % [P×K1]
        Yv   = double(reshape(Yall(:,v,:), [S K1]));   % [S×K1]
        Yhat = X * Bv;                                 % [S×K1]
        Rv   = Yv - Yhat;                              % [S×K1]

        % sign-flip residuals
        Yperm = Yhat + (Rv .* sgn);

        % refit full model under weights
        Xsw = X .* sw;
        Ysw = Yperm .* sw;
        XtX = (Xsw' * Xsw) + lambda * eye(P);
        Bp  = XtX \ (Xsw' * Ysw);

        % studentized T
        R = Yperm - X*Bp;
        WRSS = sum((R.^2) .* w, 1);
        sigma2 = WRSS ./ dof;

        XtX_inv = XtX \ eye(P);
        var_age = sigma2(:) .* XtX_inv(age_idx,age_idx);

        beta_age = Bp(age_idx,:).';
        Tperm = sum(beta_age.^2 ./ var_age);

        Tperm_all(v) = single(Tperm);

        if Tperm > Tmax
            Tmax = Tperm;
        end
    end
end

% optional cluster null (max cluster size)
if doClusterNull
    thr = quantile(double(Tperm_all), cluster_form_q);
    Tvol = nan(dim_mni,'single');
    Tvol(idx_mni) = Tperm_all;

    mask = Tvol > thr;
    CC = bwconncomp(mask,26);
    if CC.NumObjects > 0
        sizes = cellfun(@numel, CC.PixelIdxList);
        maxCluster = max(sizes);
    else
        maxCluster = 0;
    end
end
end


function T_obs = icdm_compute_Tobs_studentized(X,Yall,Kap_alpha,w_cov_all,Beta,age_idx,lambda,dof)

[S,Nmni,K1] = size(Yall);
P = size(X,2);

T_obs = zeros(Nmni,1);

for v = 1:Nmni
    Yv = double(reshape(Yall(:,v,:),[S K1]));
    Bv = double(reshape(Beta(:,v,:),[P K1]));

    w  = double(Kap_alpha(:,v) * w_cov_all(v));
    sw = sqrt(w);

    Xsw = X .* sw;
    XtX = (Xsw'*Xsw) + lambda*eye(P);
    XtX_inv = XtX \ eye(P);

    R = Yv - X*Bv;
    WRSS = sum((R.^2).*w,1);
    sigma2 = WRSS ./ dof;

    var_age = sigma2(:) .* XtX_inv(age_idx,age_idx);
    beta_age = Bv(age_idx,:).';

    T_obs(v) = sum(beta_age.^2 ./ var_age);
end
end

% ========================================================================
%  FREEDMAN–LANE PERMUTATION (WLS/WLS+RIDGE CONSISTENT)
%  - Permute residuals from reduced model (X without age) and refit full model
%  - max-T FWER p-values
%  - cluster-size null optional
%
%  Memory-safe: does NOT store full residual tensors by default; it computes
%  reduced-model residuals per voxel inside each permutation iteration.
%  This is slower but avoids huge RAM spikes/crashes.
% ========================================================================
function stats = icdm_perm_freedman_lane_wls(X,Yall,Kap,coverage,Beta,age_idx,opts)

fprintf('  -- Freedman–Lane permutation (WLS-consistent) | n_perm=%d | mode=%s\n', opts.n_perm, opts.perm_mode);

if ~isempty(opts.perm_seed)
    rng(opts.perm_seed);
end

[S,Nmni,K1] = size(Yall);
P = size(X,2);

alpha  = opts.beta_alpha;
gamma  = opts.beta_gamma;
lambda = opts.ridge_lambda;
nperm  = opts.n_perm;
Bsz    = opts.batch_size;

Kap_alpha = max(Kap, eps).^alpha;
w_cov_all = max(coverage, eps).^gamma;

dof = max(S - P, 1);

% reduced design (remove age)
X_red = X;
X_red(:,age_idx) = [];
Pr = size(X_red,2);

% -----------------------------
% 1) Observed statistic (studentized)
% -----------------------------
T_obs = icdm_compute_Tobs_studentized(X,Yall,Kap_alpha,w_cov_all,Beta,age_idx,lambda,dof);

% -----------------------------
% 2) Fit reduced model once per voxel: Yhat_red, R_red
%    (이건 Freedman–Lane에 필수)
% -----------------------------
Yhat_red = zeros(S,Nmni,K1,'single');
R_red    = zeros(S,Nmni,K1,'single');

fprintf('    fitting reduced model once...\n');
for v = 1:Nmni
    Yv = double(reshape(Yall(:,v,:),[S K1]));

    w  = double(Kap_alpha(:,v) * w_cov_all(v));
    sw = sqrt(w);

    Xsw = X_red .* sw;
    Ysw = Yv .* sw;

    XtX = (Xsw'*Xsw) + lambda*eye(Pr);
    Br  = XtX \ (Xsw'*Ysw);

    Yh = X_red * Br;

    Yhat_red(:,v,:) = single(Yh);
    R_red(:,v,:)    = single(Yv - Yh);
end

% -----------------------------
% 2a') SAVE-RESIDUALS checkpoint: dump the reduced-model residuals + all params
%      needed to run permutations externally (multi-PROCESS perm-split), then
%      return a stub. Used to bypass the MATLAB-internal parfor/MVM crash on
%      large cohorts: independent serial worker processes run perm slices.
% -----------------------------
if isfield(opts,'save_residuals') && ~isempty(opts.save_residuals)
    FLR = struct();
    FLR.X=X; FLR.Yhat_red=Yhat_red; FLR.R_red=R_red; FLR.Kap_alpha=Kap_alpha;
    FLR.w_cov_all=w_cov_all; FLR.age_idx=age_idx; FLR.lambda=lambda; FLR.dof=dof;
    FLR.Bsz=Bsz; FLR.T_obs=T_obs; FLR.Nmni=Nmni; FLR.nperm=nperm;
    FLR.dim_mni=opts.dim_mni; FLR.idx_mni=opts.idx_mni; FLR.cluster_form_q=opts.cluster_form_q;
    FLR.perm_seed=opts.perm_seed;
    FLR.tfce=struct('E',getfield_default(opts,'tfce_E',0.5),'H',getfield_default(opts,'tfce_H',2.0),...
                    'nsteps',getfield_default(opts,'tfce_nsteps',50),'conn',getfield_default(opts,'tfce_conn',26));
    save(opts.save_residuals,'-struct','FLR','-v7.3');
    fprintf('  -- saved FL residuals -> %s (run perms externally, then combine)\n', opts.save_residuals);
    stats=struct('T_obs',T_obs,'residuals_saved',opts.save_residuals);
    return;
end

% -----------------------------
% 2b) TFCE branch (opt-in): max-T and TFCE nulls on the SAME permutations
%     (seed preserved from rng() above; reduced fit uses no RNG). Returns early.
% -----------------------------
if isfield(opts,'do_tfce') && opts.do_tfce
    tfceE = getfield_default(opts,'tfce_E',0.5);
    tfceH = getfield_default(opts,'tfce_H',2.0);
    tfceN = getfield_default(opts,'tfce_nsteps',50);
    tfceC = getfield_default(opts,'tfce_conn',26);
    idxm  = opts.idx_mni;  dimm = opts.dim_mni;  cfq = opts.cluster_form_q;
    tfce_dh = double(max(T_obs)) / tfceN;
    doClusterNull = isfield(opts,'do_cluster') && opts.do_cluster;
    nWk = getfield_default(opts,'par_workers',0);   % >0 => parallel on nWk workers
    if nWk>0, modestr=sprintf('PARALLEL %d workers',nWk); else modestr='SERIAL'; end
    fprintf('  -- FL COMPLETE inference (max-T + FDR + cluster + TFCE) | E=%.2f H=%.2f dh=%.4g | n_perm=%d | %s\n',...
            tfceE,tfceH,tfce_dh,nperm, modestr);
    if ~isempty(opts.perm_seed), rng(opts.perm_seed); end
    PERM = cell(nperm,1); for m=1:nperm, PERM{m}=randperm(S); end   % same perms as serial
    T_null_max       = zeros(nperm,1,'single');
    tfce_null_max    = zeros(nperm,1,'single');
    max_cluster_null = zeros(nperm,1,'single');
    exceed_cnt       = zeros(Nmni,1,'uint32');   % voxelwise exceedances -> BH-FDR

    if nWk>0
        clear Yall Kap coverage Beta;             % free large client arrays before pool
        delete(gcp('nocreate')); parpool('Processes', nWk);
        Cyh=parallel.pool.Constant(Yhat_red); Crr=parallel.pool.Constant(R_red); Ckap=parallel.pool.Constant(Kap_alpha);
        Tobs_c=parallel.pool.Constant(T_obs);
        parfor m = 1:nperm
            [Tm, Tperm_all, maxCl] = icdm_fl_oneperm( ...
                X, Cyh.Value, Crr.Value, Ckap.Value, w_cov_all, age_idx, lambda, dof, ...
                Bsz, doClusterNull, dimm, idxm, cfq, PERM{m});
            T_null_max(m)       = single(Tm);
            max_cluster_null(m) = single(maxCl);
            tv                  = icdm_tfce(Tperm_all, idxm, dimm, tfceE, tfceH, tfce_dh, tfceC);
            tfce_null_max(m)    = single(max(tv));
            exceed_cnt          = exceed_cnt + uint32(Tperm_all >= Tobs_c.Value);   % reduction
        end
    else
        for m = 1:nperm
            if rem(m,50)==0, fprintf('    perm %d/%d\n', m, nperm); end
            [Tm, Tperm_all, maxCl] = icdm_fl_oneperm( ...
                X, Yhat_red, R_red, Kap_alpha, w_cov_all, age_idx, lambda, dof, ...
                Bsz, doClusterNull, dimm, idxm, cfq, PERM{m});
            T_null_max(m)=single(Tm); max_cluster_null(m)=single(maxCl);
            tv=icdm_tfce(Tperm_all, idxm, dimm, tfceE, tfceH, tfce_dh, tfceC);
            tfce_null_max(m)=single(max(tv));
            exceed_cnt = exceed_cnt + uint32(Tperm_all >= T_obs);
        end
    end

    tfce_obs = icdm_tfce(T_obs, idxm, dimm, tfceE, tfceH, tfce_dh, tfceC);
    p_fwer = zeros(Nmni,1); tfce_p_fwer = zeros(Nmni,1);
    for v = 1:Nmni
        p_fwer(v)      = (1 + sum(T_null_max    >= T_obs(v)))    / (nperm + 1);
        tfce_p_fwer(v) = (1 + sum(tfce_null_max >= tfce_obs(v))) / (nperm + 1);
    end
    stats = struct();
    stats.T_obs=T_obs; stats.T_null_max=T_null_max;
    stats.p_map=p_fwer; stats.p_fwer=p_fwer;
    stats.p_uncorr = (double(exceed_cnt)+1)/(nperm+1);   % for BH-FDR (section 5)
    stats.tfce_obs=tfce_obs; stats.tfce_null_max=tfce_null_max; stats.tfce_p_fwer=tfce_p_fwer;
    stats.tfce_params=struct('E',tfceE,'H',tfceH,'nsteps',tfceN,'conn',tfceC,'dh',tfce_dh);
    if doClusterNull, stats.max_cluster_null=max_cluster_null; else stats.max_cluster_null=[]; end
    fprintf('  -- FL complete inference done (max-T, FDR, cluster, TFCE on identical perms)\n');
    return;
end

% -----------------------------
% 3) Permutations
% -----------------------------
T_null_max = zeros(nperm,1,'single');
exceed_cnt = zeros(Nmni,1,'uint32'); % voxelwise p_uncorr for BH-FDR
doClusterNull = isfield(opts,'do_cluster') && opts.do_cluster;
max_cluster_null = zeros(nperm,1,'single');

usePar = false; % 정확한 voxelwise p+FDR까지 원하면 for가 안전
if isfield(opts,'usePar') && ~isempty(opts.usePar)
    usePar = opts.usePar;
end

if usePar
    warning('Freedman–Lane: voxelwise permutation p(FDR)까지 하려면 for-loop가 가장 안전합니다. 현재는 max-T만 병렬로 수행합니다.');
    parfor m = 1:nperm
        [Tmax_m, ~, ~] = icdm_freedmanlane_oneperm_stats( ...
            X, Yhat_red, R_red, Kap_alpha, w_cov_all, age_idx, lambda, dof, ...
            Bsz, false, opts.dim_mni, opts.idx_mni, opts.cluster_form_q);
        T_null_max(m) = single(Tmax_m);
    end
else
    for m = 1:nperm
        if rem(m,50)==0, fprintf('    perm %d/%d\n', m, nperm); end

        idx = randperm(S);

        [Tmax_m, Tperm_all, maxCl_m] = icdm_freedmanlane_oneperm_stats( ...
            X, Yhat_red, R_red, Kap_alpha, w_cov_all, age_idx, lambda, dof, ...
            Bsz, doClusterNull, opts.dim_mni, opts.idx_mni, opts.cluster_form_q, idx);

        T_null_max(m) = single(Tmax_m);

        exceed_cnt = exceed_cnt + uint32(Tperm_all >= T_obs);

        if doClusterNull
            max_cluster_null(m) = single(maxCl_m);
        end
    end
end

% -----------------------------
% 4) p-values
% -----------------------------
p_fwer = zeros(Nmni,1);
for v = 1:Nmni
    p_fwer(v) = (1 + sum(T_null_max >= T_obs(v))) / (nperm + 1);
end

if ~usePar
    p_uncorr = (double(exceed_cnt) + 1) / (nperm + 1);
else
    p_uncorr = nan(Nmni,1);
end

stats = struct();
stats.T_obs = T_obs;
stats.T_null_max = T_null_max;

stats.p_map = p_fwer;      % 기존 인터페이스: p_map = max-T FWER
stats.p_fwer = p_fwer;
stats.p_uncorr = p_uncorr;

if doClusterNull && ~usePar
    stats.max_cluster_null = max_cluster_null;
else
    stats.max_cluster_null = [];
end

fprintf('  -- Freedman–Lane done\n');
end


function [Tmax, Tperm_all, maxCluster] = icdm_freedmanlane_oneperm_stats( ...
    X, Yhat_red, R_red, Kap_alpha, w_cov_all, age_idx, lambda, dof, ...
    batch_size, doClusterNull, dim_mni, idx_mni, cluster_form_q, idx)

[S,Nmni,K1] = size(Yhat_red);
P = size(X,2);

if nargin < 14 || isempty(idx)
    idx = randperm(S);
end

Tmax = 0;
maxCluster = 0;

Tperm_all = zeros(Nmni,1,'single');

nbatch = ceil(Nmni / batch_size);

for b = 1:nbatch
    v_start = (b-1)*batch_size + 1;
    v_end   = min(b*batch_size, Nmni);

    for v = v_start:v_end

        Yh = double(squeeze(Yhat_red(:,v,:)));   % [S×K1]
        Rv = double(squeeze(R_red(:,v,:)));      % [S×K1]

        % Freedman–Lane: permute reduced residuals
        Yperm = Yh + Rv(idx,:);

        % voxel weights
        w  = double(Kap_alpha(:,v) * w_cov_all(v));
        sw = sqrt(w);

        % full model refit
        Xsw = X .* sw;
        Ysw = Yperm .* sw;

        XtX = (Xsw'*Xsw) + lambda*eye(P);
        Bp  = XtX \ (Xsw'*Ysw);

        % studentized T
        R = Yperm - X*Bp;
        WRSS = sum((R.^2).*w,1);
        sigma2 = WRSS ./ dof;

        XtX_inv = XtX \ eye(P);
        var_age = sigma2(:) .* XtX_inv(age_idx,age_idx);

        beta_age = Bp(age_idx,:).';
        Tperm = sum(beta_age.^2 ./ var_age);

        Tperm_all(v) = single(Tperm);

        if Tperm > Tmax
            Tmax = Tperm;
        end
    end
end

% optional cluster null
if doClusterNull
    thr = quantile(double(Tperm_all), cluster_form_q);

    Tvol = nan(dim_mni,'single');
    Tvol(idx_mni) = Tperm_all;

    mask = Tvol > thr;
    CC = bwconncomp(mask,26);
    if CC.NumObjects > 0
        sizes = cellfun(@numel, CC.PixelIdxList);
        maxCluster = max(sizes);
    else
        maxCluster = 0;
    end
end
end

% ========================================================================
%  CLUSTER CORRECTION (cluster-size FWER using permutation max-cluster null)
% ========================================================================
function cluster_map = icdm_cluster_correction(stats, opts)

alpha = opts.cluster_alpha;

% reconstruct T volume
Tvol = nan(opts.dim_mni);
Tvol(opts.idx_mni) = stats.T_obs;

% cluster-forming threshold: use quantile of null max-T (conservative),
% OR use stats.T_obs quantile if you prefer. Here we tie it to null max-T.
thr_form = quantile(stats.T_null_max, opts.cluster_form_q);

mask = Tvol > thr_form;
CC = bwconncomp(mask, 26);

cluster_map = zeros(opts.dim_mni);

if CC.NumObjects == 0
    return;
end

for c = 1:CC.NumObjects
    sz = numel(CC.PixelIdxList{c});

    % cluster-size p-value (FWER) using null distribution of max cluster size
    p_cluster = (1 + sum(stats.max_cluster_null >= sz)) / (numel(stats.max_cluster_null) + 1);

    if (p_cluster < alpha) && (sz >= opts.cluster_min_size)
        cluster_map(CC.PixelIdxList{c}) = 1;
    end
end

end

% ========================================================================
%  BUILD DESIGN MATRIX
% ========================================================================
function [X, S, P] = icdm_build_design(design_input, opts)

if isstruct(design_input)
    subjects = design_input;
    S = numel(subjects);

    x1 = icdm_design_vector(subjects(1), opts);
    P  = numel(x1);

    X = zeros(S, P, 'double');
    X(1,:) = double(x1(:)).';

    for s = 2:S
        xs = icdm_design_vector(subjects(s), opts);
        if numel(xs) ~= P
            error('Design vector length mismatch at subject %d.', s);
        end
        X(s,:) = double(xs(:)).';
    end

elseif isnumeric(design_input)
    X = double(design_input);
    [S, P] = size(X);

else
    error('design_input must be subjects struct array or numeric design matrix.');
end

end

% ========================================================================
%  DESIGN VECTOR (minimal default; customize as needed)
% ========================================================================
function x = icdm_design_vector1(subj, opts)
% Supports opts.design_spec entries: 'const','age','sex'
spec = opts.design_spec;
P = numel(spec);
x = zeros(P,1);

for j = 1:P
    name = lower(string(spec{j}));
    switch name
        case "const"
            x(j) = 1;

        case "age"
            if ~isfield(subj,'age')
                error('Subject missing field .age');
            end
            mu = opts.design.mu.age;
            sd = opts.design.sd.age;
            if sd == 0, sd = 1; end
            x(j) = (double(subj.age) - mu) / sd;

        case "sex"
            if isfield(subj,'sex')
                x(j) = double(subj.sex);
            else
                x(j) = 0;
            end

        otherwise
            x(j) = 0; % unknown covariate: user can extend
    end
end

x = double(x(:));
end

% ========================================================================
%  STACK DATA: group -> Yall [S×Nmni×K1], Kap [S×Nmni]
% ========================================================================
function [Yall, Kap] = icdm_stack_data(group, S)

if numel(group) ~= S
    error('Group size mismatch: numel(group)=%d, S=%d', numel(group), S);
end

Y0 = group{1}.y_ilr_mni;
[Nmni, K1] = size(Y0);

Yall = zeros(S, Nmni, K1, 'single');
Kap  = zeros(S, Nmni, 'single');

for s = 1:S
    if ~isfield(group{s},'y_ilr_mni') || ~isfield(group{s},'kappa_mni')
        error('group{%d} missing y_ilr_mni or kappa_mni', s);
    end

    Ys = group{s}.y_ilr_mni;
    Ks = group{s}.kappa_mni;

    if size(Ys,1) ~= Nmni || size(Ys,2) ~= K1
        error('Subject %d: y_ilr_mni must be [Nmni×K1].', s);
    end
    if numel(Ks) ~= Nmni
        error('Subject %d: kappa_mni must have Nmni elements.', s);
    end

    Yall(s,:,:) = single(Ys);
    Kap(s,:)    = single(Ks(:)).';
end

end

% ========================================================================
%  BH-FDR
% ========================================================================
function [q_map, sig_mask, p_crit] = bh_fdr(p_map, q)

if nargin < 2 || isempty(q)
    q = 0.05;
end

orig_size = size(p_map);
p = p_map(:);

valid = isfinite(p);
p_valid = p(valid);

m = numel(p_valid);
if m == 0
    q_map = nan(orig_size);
    sig_mask = zeros(orig_size);
    p_crit = NaN;
    return
end

[ps, sort_idx] = sort(p_valid);
rank = (1:m)';
bh_line = (rank/m) * q;

below = find(ps <= bh_line);

q_map = ones(size(p));
sig_mask = zeros(size(p));

if ~isempty(below)
    max_i = max(below);
    p_crit = ps(max_i);

    sig_mask(valid) = p_valid <= p_crit;

    q_vals = ps .* m ./ rank;
    q_vals = min(1, flipud(cummin(flipud(q_vals))));

    q_temp = zeros(m,1);
    q_temp(sort_idx) = q_vals;
    q_map(valid) = q_temp;
else
    p_crit = 0;
end

q_map = reshape(q_map, orig_size);
sig_mask = reshape(sig_mask, orig_size);

end

% ========================================================================
% util
% ========================================================================
function val = getfield_default(s,field,default)
if isfield(s,field) && ~isempty(s.(field))
    val = s.(field);
else
    val = default;
end
end