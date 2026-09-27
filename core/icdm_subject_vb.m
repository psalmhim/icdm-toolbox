function OUT = icdm_subject_vb(subj, K, grp, opts)
% ========================================================================
% icdm_subject_vb.m  (covariate-free MAP/Laplace implementation)
%
%   - Loads native ICDM counts
%   - Warps group prior μ, κ from MNI -> native
%   - Covariate-free population prior (covariates are downstream only)
%   - Voxel-wise MAP estimation with a Laplace covariance approximation
%   - Warps posterior back to MNI
%   - Debug visualization if opts.verbose
% ========================================================================

%fprintf('[VB] Subject %s\n', subj.id);
H  = grp.H;
K1 = K - 1;

% 1. LOAD NATIVE ICDM COUNTS
if isfield(subj,'datafile') && ~isempty(subj.datafile)
    load(subj.datafile, 'icdm2d', 'warp','idx_native', 'dim_native','V');
    C_2d=icdm2d(:,subj.idx_regions);
    clear icdm2d;
    Vnative=V;
else
    V  = spm_vol(subj.icdm_4d); V=V(subj.idx_regions);
    C4 = spm_read_vols(V);
    thr = getfield_default(opts,'thresh',5);
    sumC = sum(C4,4);
    mask_3d = sumC > thr;
    idx_native = find(mask_3d(:));
    icdm2d = reshape(C4, [], size(C4,4));
    C_2d = icdm2d(idx_native, :);
    Vnative=V(1);   
    dim_native = Vnative.dim;
    warp.to_native = icdm_build_warp_to_native( ...
                C4, subj.dartel_flow,subj.template);
    
    warp.to_mni = icdm_build_warp_to_mni( ...
                C4, subj.dartel_flow, subj.template);
    warp=warp;
    clear C4 icdm2d;
end

Nnat = numel(idx_native);
if size(C_2d,1) ~= Nnat || size(C_2d,2) ~= K
    error('Native count matrix must be [%d x %d], but is [%d x %d].', ...
        Nnat, K, size(C_2d,1), size(C_2d,2));
end
if any(~isfinite(C_2d(:))) || any(C_2d(:) < 0)
    error('ICDM counts must be finite and nonnegative.');
end
gamma = getfield_default(opts,'gamma',1);
if ~isscalar(gamma) || ~isfinite(gamma) || gamma <= 0
    error('opts.gamma must be a positive finite scalar.');
end
Ct_2d    = (C_2d + 1).^gamma;
% Normalise to tempered compositions π^(γ)
Ct_2d = Ct_2d ./ max(sum(Ct_2d,2), eps);
Ct_2d(~isfinite(Ct_2d)) = 0;


% 2. WARP GROUP PRIOR FROM MNI → NATIVE
dim_mni = grp.dim_mni;
idx_mni = grp.idx_mni;
fprintf('  Warp group prior MNI -> native...\n');
% --- μ ---% --- κ ---
is_initial = isempty(grp.mu_ilr_mni) || ...
    (all(isfinite(grp.mu_ilr_mni(:))) && all(grp.mu_ilr_mni(:) == 0));
if is_initial
    MU_native= zeros(Nnat,K1,'single');
    KAP_native = opts.precision.kappa_base * ones(Nnat,1,'single');
else
    MU_native=icdm_warp_to_native(grp.mu_ilr_mni,warp,idx_native,idx_mni);
    KAP_native=icdm_warp_to_native(grp.kappa_mni,warp,idx_native,idx_mni);
    KAP_native(~isfinite(KAP_native) | KAP_native<=0) = opts.precision.kappa_base;
end

% -------------------------------------------------------------------------
% General mode (Sec. 2.5.3): OPTIONAL covariate-informed prior MEAN.
% Default (opts.cov_prior empty): covariate-agnostic testing mode, m_v = mu_grp,
%   used for ALL reported covariate-effect tests -- Beta^assoc estimated downstream.
% If a cross-fitted / reference predictive coefficient beta^pred (deviation form,
% from the data-only representation) is supplied, shift the prior MEAN toward the
% subject's covariate prediction:  m_v = mu_grp + (beta^pred)' x   (Eq. cov_prior).
% MEAN-ONLY: the prior precision KAP_native is UNCHANGED (no precision inflation).
% -------------------------------------------------------------------------
cov = getfield_default(opts,'cov_prior',[]);
if ~is_initial && ~isempty(cov) && isfield(cov,'Beta') && ~isempty(cov.Beta) ...
        && isfield(cov,'x') && ~isempty(cov.x)
    xs   = double(cov.x(:))';                     % [1 x P] centered covariate for this subject
    Beta = cov.Beta;                              % [Nmni x K1 x P] deviation-form beta^pred (grp layout)
    Nmni = numel(grp.idx_mni);
    pred_mni = zeros(Nmni,K1,'single');
    for p = 1:size(Beta,3)
        pred_mni = pred_mni + single( xs(p) * double(Beta(:,:,p)) ); % sum_p x_p * beta^pred_p  (deviation)
    end
    off_native = icdm_warp_to_native(pred_mni,warp,idx_native,idx_mni);
    off_native(~isfinite(off_native)) = 0;
    MU_native = MU_native + off_native;           % m_v = mu_grp + (beta^pred)' x  (mean only)
    fprintf('  [cov-prior] covariate-informed prior MEAN applied (precision unchanged).\n');
end

fprintf('  Voxel-wise MAP/Laplace inference...\n');

% 5. VOXELWISE MAP + LAPLACE CURVATURE
vb_opts = opts.vb;

y_nat   = zeros(Nnat,K1,'single');
y_data_nat = zeros(Nnat,K1,'single');             % data-only ILR (prior-independent) for group mean/dispersion
kappa_v = zeros(Nnat,1,'single');
kappa_data_v = zeros(Nnat,1,'single');            % data-only reliability (prior-independent)
kappa_base=opts.precision.kappa_base;
alpha_grp = getfield_default(opts.precision,'alpha_grp',1);   % global group-transfer scale (default 1)
% Group-prior precision is Lambda_v^(s) = alpha_grp * kappa_v^grp * I  (reliability factor omega == 1;
% the optional phi_Les/phi_N tempering and any baseline elevation are not applied in the reported estimator).
H2col = H.^2;                                     % [K x K1] for data-curvature

parfor v = 1:Nnat
    % ---------------------------------------------------------
    % 1. Likelihood (raw counts only)
    % ---------------------------------------------------------
    n_v = double(C_2d(v,:)');      % RAW streamline counts
    n_v(~isfinite(n_v) | n_v < 0) = 0;
    N_v = sum(n_v);

    % ---------------------------------------------------------
    % 2. Prior at voxel v
    % ---------------------------------------------------------
    m_v = MU_native(v,:)';
    P_v = (alpha_grp*KAP_native(v))*ones(K1,1,'double');

    if ~isfinite(N_v) || N_v<0, N_v=0; end
    m_v(~isfinite(m_v)) = 0;
    P_v(~isfinite(P_v)|P_v<=0) = kappa_base;

    % ---------------------------------------------------------
    % (3) ILR initialisation using γ-tempered composition
    %     (Ct_2d(v,:) comes from (C4+1)^γ and is normalized)
    % ---------------------------------------------------------
    y0 = y0_from(m_v, Ct_2d(v,:), H);

    % ---------------------------------------------------------
    % 4. Newton–Raphson MAP update
    % ---------------------------------------------------------
    
    [y_row, kap] = solve_voxel_newton( ...
        n_v, N_v, m_v, P_v, H, vb_opts, y0);

    if any(~isfinite(y_row))
        y_row = m_v';
    end
    if ~isfinite(kap), kap = mean(P_v); end

    y_nat(v,:) = y_row;
    kappa_v(v) = kap;
    % data-only reliability: kappa_data = tr(N H'S(pi_data)H)/(K-1) at the count-derived composition (NO prior)
    pd = n_v / max(N_v,1e-9);
    hpd = H' * pd;
    data_diag = N_v * max(sum(H2col.*pd,1)' - hpd.^2, 0);
    kappa_data_v(v) = single(mean(data_diag));
    % data-only ILR (prior-independent) for group mean & dispersion (avoids mu/tau feedback)
    p_reg = (n_v+0.5) / (N_v + 0.5*K);
    y_data_nat(v,:) = single((log(p_reg))'*H);
end


% 6. WARP POSTERIOR → MNI
fprintf('  Warp posterior to MNI space...\n');
y_ilr_mni =icdm_warp_to_mni(y_nat,warp,idx_mni,idx_native);
y_ilr_mni (~isfinite(y_ilr_mni )) = 0;
kappa_mni =icdm_warp_to_mni(kappa_v,warp,idx_mni,idx_native);
kappa_mni(~isfinite(kappa_mni) | kappa_mni<=0) = opts.precision.kappa_base;
kappa_data_mni = icdm_warp_to_mni(kappa_data_v,warp,idx_mni,idx_native);
kappa_data_mni(~isfinite(kappa_data_mni) | kappa_data_mni<=0) = opts.precision.kappa_base;
y_data_mni = icdm_warp_to_mni(y_data_nat,warp,idx_mni,idx_native);
support_native = ones(Nnat,1,'single');
support_mni = icdm_warp_to_mni(support_native,warp,idx_mni,idx_native);
support_mni = support_mni(:);
support_thresh = getfield_default(opts,'warp_support_thresh',0.5);
valid_mni = isfinite(support_mni) & support_mni >= support_thresh;
y_data_mni(~isfinite(y_data_mni)) = NaN;
y_data_mni(~repmat(valid_mni,1,K1)) = NaN;
kappa_data_mni(~valid_mni) = NaN;
y_ilr_mni(~repmat(valid_mni,1,K1)) = NaN;
kappa_mni(~valid_mni) = NaN;

% 7. DEBUG VISUALIZATION

% if isfield(opts,'verbose') && opts.verbose
%     debug_subject_plots(subj, ...
%         MU_nat_4D, KAP_nat_3D, ...
%         Y_mni_4D, K_mni_3D);
% end


% 8. RETURN STRUCT
OUT.id            = subj.id;
%OUT.idx_native    = idx_native;
%OUT.y_ilr_native  = y_nat;
%OUT.kappa_native  = kappa_v;
OUT.y_ilr_mni     = y_ilr_mni;
OUT.kappa_mni     = kappa_mni;          % posterior reliability (data + prior) — OUTPUT/reporting only
OUT.kappa_data_mni = kappa_data_mni;    % data-only reliability for QC/reporting (not population homogeneity)
OUT.y_data_mni    = y_data_mni;         % data-only ILR — used for group mean & dispersion (no mu/tau feedback)
OUT.valid_mni     = valid_mni;

end  % -------------------- END MAIN FUNCTION -----------------------------


%% DEBUG PLOTS FOR SUBJECT
function debug_subject_plots(subj, MU4, KAP3,Yfull,KAP_nat)
dimN=size(MU4);
[Xn,Yn,Zn] = deal(dimN(1),dimN(2),dimN(3));
mid = round(Zn/2);
try
    outdir = fileparts(subj.icdm_4d);
    fig = figure('Visible','off','Position',[100 100 1500 900]);

    tiledlayout(3,3);

    % -------- Original MU -----------
    nexttile;
    imagesc(MU4(:,:,mid,1)); axis image off;
    title('MU nat (d=1)');

    % -------- Original Kappa --------
    nexttile;
    imagesc(KAP3(:,:,mid)); axis image off;
    title('Kappa nat');


    % -------- Posterior y (ILR) -----
    nexttile;
    imagesc(Yfull(:,:,mid,1)); axis image off;
    title('Posterior ILR d=1');

    % -------- Posterior y (ILR) -----
    nexttile;
    imagesc(KAP_nat(:,:,mid)); axis image off;
    title('Posterior Kappa');

    saveas(fig, fullfile(outdir, sprintf('%s_debug.png',subj.id)));
    close(fig);

catch ME
    warning('debug plot failed: %s', ME.message);
end

end


function [y_row, kappa_scalar] = solve_voxel_newton(n, N, m, Pdiag, H, vb, y0)

K1 = numel(m);
H = double(H);

Pdiag = max(double(Pdiag),1e-8);
m = double(m(:));
n = double(n(:));

if nargin < 7 || isempty(y0)
    y = m;
else
    y = double(y0(:));
end

for it = 1:vb.max_iter

    z = H*y;
    [~, mu] = logsumexp_softmax(z);
    mu = max(mu,1e-12);

    g = H'*(n - N*mu) - (Pdiag.*(y - m));

    if norm(g,2) < vb.tol_grad
        break;
    end

    pcg_tol = max(vb.tol_grad, min(1e-1, 0.1*norm(g,2)));
    pcg_iter = getfield_default(vb,'pcg_iter',20);
    [step, ~, ~] = pcg_newton_step(H, mu, Pdiag, N, g, pcg_tol, pcg_iter);
    if norm(step,2) < vb.tol_step
        break;
    end

    alpha = 1.0;
    F0 = voxel_logpost(y,n,N,m,Pdiag,H);
    accepted = false;
    for b = 1:10
        y_try = y + alpha*step;
        F_try = voxel_logpost(y_try,n,N,m,Pdiag,H);
        if F_try >= F0
            y = y_try;
            accepted = true;
            break;
        end
        alpha = alpha*0.5;
    end
    if ~accepted, break; end
end

z = H*y;
[~, mu] = logsumexp_softmax(z);
mu = max(mu,1e-12);
Hd2 = sum((H.^2).*mu,1)';
Hmu = H' * mu;
Qdiag = Pdiag + N*max(Hd2 - Hmu.^2,0) + 1e-12;

kappa_scalar = single(mean(Qdiag));
y_row = single(y');
end


%% helper: PCG Newton step

function [step, iters, ok] = pcg_newton_step(H, mu, Pdiag, N, g, tol, maxit)

Km1 = numel(Pdiag);
H = double(H);
mu = double(mu(:));
g  = double(g(:));
Pdiag = double(Pdiag(:));

Htmu = H' * mu;                               % (K-1)x1
Hd2  = sum((H.^2).*mu,1)';                    % diag(H' diag(mu) H)
Mdiag = Pdiag + N*Hd2;
Minv  = 1 ./ max(Mdiag, 1e-12);

Qx = @(x) (Pdiag.*x) + N*( H'*(mu.*(H*x)) - (Htmu*(Htmu'*x)) );

x = zeros(Km1,1);
r = g - Qx(x);
z = Minv .* r;
p = z;

rz_old = r'*z;
ok = false;
iters = 0;

for k=1:maxit
    Ap = Qx(p);
    denom = p'*Ap;
    if denom <= 0, break; end

    alpha = rz_old / denom;

    x = x + alpha*p;
    r = r - alpha*Ap;

    if norm(r,2) <= tol
        ok = true; iters = k;
        break;
    end

    z = Minv .* r;
    rz_new = r'*z;
    beta = rz_new / max(rz_old,1e-20);
    p = z + beta*p;

    rz_old = rz_new;
    iters = k;
end

step = x;
end



%% helper: F(y)

function Fv = voxel_logpost(y, n, N, m, Pdiag, H)

z = H*y;
[lse, ~] = logsumexp_softmax(z);

ll = n'*z - N*lse;
lp = -0.5*(y-m)'*(Pdiag.*(y-m));
Fv = ll + lp;

end



%% helper: logsumexp + softmax

function [lse, mu] = logsumexp_softmax(z)
mz = max(z);
a = exp(z - mz);
s = sum(a);
mu = a / s;
lse = log(s) + mz;
end



%% y0 initialization from empirical proportions

function y0 = y0_from(m, p, H)
p = double(p(:));
pi0 = p / max(sum(p),realmin('double'));

z0 = log(max(pi0, realmin('double')));
z0 = z0 - mean(z0);
y_emp = H' * z0;

y0 = 0.7*y_emp + 0.3*double(m(:));
end
