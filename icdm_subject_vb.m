% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
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

function OUT = icdm_subject_vb(subj, K, grp, opts)
% ICDM_SUBJECT_VB  Subject-level Newton-Laplace MAP inference in ILR space.
%
%   OUT = icdm_subject_vb(subj, K, grp, opts)
%
%   Performs voxelwise MAP estimation of ILR coordinates for a single
%   subject, implementing the E-step of the iterative empirical Bayes
%   loop (manuscript Section 2: Newton-Laplace inference in ILR space).
%
%   For each voxel v, maximises the log-posterior (Eq. log-posterior):
%     L(y) = ell(n|y) - 0.5*(y-m)'*Lambda*(y-m)
%   using Newton-CG with backtracking line search (Eqs. grad-voxel,
%   hess-voxel, newton-linear-system, newton-update).
%
%   Posterior precision is approximated via Laplace (Eq. laplace-posterior)
%   and summarised as kappa = tr(Q)/(K-1) (Eq. def-kappa).
%
%   Pipeline:
%     1. Load native ICDM streamline counts
%     2. Warp group prior mu, kappa from MNI -> native space
%     3. Optional: blend covariate-predicted prior (Eq. subject_prior)
%     4. Voxelwise Newton-Laplace MAP with gamma-tempered initialisation
%        (Eq. gamma_tempered, Eq. initial_y)
%     5. Warp posterior y_ilr, kappa back to MNI space
%
%   INPUT
%     subj : struct with fields .id, .datafile or .icdm_4d, .dartel_flow
%     K    : number of compositional components (e.g. 68)
%     grp  : group prior struct (.mu_ilr_mni, .kappa_mni, .H, .Beta, ...)
%     opts : options struct (.gamma, .vb, .precision, .w_beta, ...)
%
%   OUTPUT
%     OUT  : struct with fields
%       .id          : subject identifier
%       .y_ilr_mni   : [Nmni x (K-1)] MAP ILR coordinates in MNI space
%       .kappa_mni   : [Nmni x 1] posterior reliability kappa in MNI space
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_POPULATION_EB, ICDM_UPDATE_GROUP_PRIOR, HELMERT_SUBMATRIX

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

[Xn, Yn, Zn] = deal(dim_native(1), dim_native(2), dim_native(3));
Nnat = numel(idx_native);
gamma = opts.gamma;
Ct_2d    = (C_2d + 1).^gamma;
% Normalise to tempered compositions π^(γ)
Ct_2d = Ct_2d ./ max(sum(Ct_2d,2), eps);
Ct_2d(~isfinite(Ct_2d)) = 0;


% 2. WARP GROUP PRIOR FROM MNI → NATIVE
dim_mni = grp.dim_mni;
idx_mni = grp.idx_mni;
fprintf('  Warp group prior MNI -> native...\n');
% --- μ ---% --- κ ---
if sum(abs(grp.mu_ilr_mni(:))) == 0 %first iteration
    MU_native= zeros(Nnat,K1,'single');
    KAP_native = opts.precision.kappa_base * ones(Nnat,1,'single');
else
    MU_native=icdm_warp_to_native(grp.mu_ilr_mni,warp,idx_native,idx_mni);
    KAP_native=icdm_warp_to_native(grp.kappa_mni,warp,idx_native,idx_mni);
    KAP_native(~isfinite(KAP_native) | KAP_native<=0) = opts.precision.kappa_base;
end

% 4. β-FUSION (unchanged)
w_beta = getfield_default(opts,'w_beta',0);
if w_beta > 0
    w_beta = min(max(w_beta,0),1);
    x = icdm_design_vector(subj, opts);   % [1 × P]
    Nmni       = numel(grp.idx_mni);
    pred_mni = zeros(Nmni, K1, 'single');
    for d = 1:K1
        Bd = double(squeeze(grp.Beta(:,:,d)));      % [P × Nmni]
        pred_mni(:,d) = single((x * Bd).');
    end
    pMU_wbeta=icdm_warp_to_native(pred_mni,warp,idx_native,idx_mni);
    if ~isempty(pMU_wbeta)
        pMU_wbeta(~isfinite(pMU_wbeta)) = 0;
        rho = w_beta/(1-w_beta);
        KAP_native = (1+rho).*KAP_native;
        MU_native  = (1-w_beta).*MU_native + w_beta.*pMU_wbeta;
    end
end

fprintf('  Voxel-wise VB inference...\n');

% 5. VOXELWISE VB (Newton-Laplace)
vb_opts = opts.vb;

y_nat   = zeros(Nnat,K1,'single');
kappa_v = zeros(Nnat,1,'single');
kappa_base=opts.precision.kappa_base;

parfor v = 1:Nnat
    % ---------------------------------------------------------
    % 1. Likelihood (raw counts only)
    % ---------------------------------------------------------
    n_v = C_2d(v,:)';      % RAW streamline counts
    N_v = sum(n_v);

    % ---------------------------------------------------------
    % 2. Prior at voxel v
    % ---------------------------------------------------------
    m_v = MU_native(v,:)';
    P_v = KAP_native(v)*ones(K1,1,'double');

    n_v(~isfinite(n_v)) = 0;
    if ~isfinite(N_v) || N_v<0, N_v=0; end
    m_v(~isfinite(m_v)) = 0;
    P_v(~isfinite(P_v)|P_v<=0) = kappa_base;

    % ---------------------------------------------------------
    % (3) ILR initialisation using γ-tempered composition
    %     (Ct_2d(v,:) comes from (C4+1)^γ and is normalized)
    % ---------------------------------------------------------
    y0 = y0_from(m_v, Ct_2d(v,:), H);

    % ---------------------------------------------------------
    % 4. Newton–Raphson VB update
    % ---------------------------------------------------------
    
    [y_row, kap] = solve_voxel_newton( ...
        n_v, N_v, m_v, P_v, H, vb_opts, y0);

    if any(~isfinite(y_row))
        y_row = m_v';
    end
    if ~isfinite(kap), kap = mean(P_v); end

    y_nat(v,:) = y_row;
    kappa_v(v) = kap;
end


% 6. WARP POSTERIOR → MNI
fprintf('  Warp posterior to MNI space...\n');
y_ilr_mni =icdm_warp_to_mni(y_nat,warp,idx_mni,idx_native);
y_ilr_mni (~isfinite(y_ilr_mni )) = 0;
kappa_mni =icdm_warp_to_mni(kappa_v,warp,idx_mni,idx_native);
kappa_mni(~isfinite(kappa_mni) | kappa_mni<=0) = opts.precision.kappa_base;

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
OUT.kappa_mni     = kappa_mni;

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
    F0 = voxel_obj(y,n,N,m,Pdiag,H);
    for b = 1:10
        y_try = y + alpha*step;
        F_try = voxel_obj(y_try,n,N,m,Pdiag,H);
        if F_try >= F0
            y = y_try;
            break;
        end
        alpha = alpha*0.5;
    end
end

z = H*y;
[~, mu] = logsumexp_softmax(z);
mu = max(mu,1e-12);
Hd2 = sum((H.^2).*mu,1)';

Qdiag = Pdiag + N*Hd2 + 1e-12;

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

function Fv = voxel_obj(y, n, N, m, Pdiag, H)

z = H*y;
[lse, mu] = logsumexp_softmax(z);

ll = n'*z - N*lse;
lp = -0.5*(y-m)'*(Pdiag.*(y-m)) + 0.5*sum(log(Pdiag));

mu = max(mu,1e-12);
Hd2 = sum((H.^2).*mu,1)';
Qdiag = Pdiag + N*Hd2 + 1e-12;

Fv = ll + lp - 0.5*sum(log(Qdiag));

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

function y0 = y0_from(m, n, H)
n = double(n(:));
alpha = n + 1;
pi0 = alpha / sum(alpha);

z0 = log(max(pi0, realmin('double')));
z0 = z0 - mean(z0);
y_emp = H' * z0;

y0 = 0.7*y_emp + 0.3*double(m(:));
end