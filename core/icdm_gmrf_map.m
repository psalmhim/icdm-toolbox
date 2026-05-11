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

function [Ymap, KAP, Fvox] = icdm_gmrf_map(C, Ntot, M_ilr, P_ilr, H, Y0, dim3, idx, conn, ...
                                          lambda_gmrf, max_iter, tol_g, tol_s, jitter, ...
                                          pcg_tol, pcg_maxit)
% ICDM_GMRF_MAP  Joint MAP inference with Gaussian Markov random field prior.
%
%   [Ymap, KAP, Fvox] = icdm_gmrf_map(C, Ntot, M_ilr, P_ilr, H, Y0,
%       dim3, idx, conn, lambda_gmrf, max_iter, tol_g, tol_s, jitter,
%       pcg_tol, pcg_maxit)
%
%   Jointly optimises all voxels using a spatially coupled objective:
%     F(y) = sum_v [ell_v(y_v) + lp_v(y_v)] - (lambda/2) * y' (L x I) y
%   where ell_v is the multinomial log-likelihood in ILR (Eq. def-loglike),
%   lp_v is the Gaussian voxelwise prior, and L is the graph Laplacian
%   on the white-matter mask providing spatial regularisation.
%
%   Newton-CG with backtracking line search; PCG solves the (block-diagonal
%   + GMRF) linear system at each iteration.
%
%   INPUT
%     C           : [nI x K] streamline count matrix (masked voxels)
%     Ntot        : [nI x 1] total counts per voxel
%     M_ilr       : [nI x (K-1)] prior mean in ILR space
%     P_ilr       : [nI x (K-1)] diagonal prior precision per dimension
%     H           : [K x (K-1)] Helmert submatrix
%     Y0          : [nI x (K-1)] initial ILR coordinates
%     dim3        : [X Y Z] volume dimensions
%     idx         : linear indices of masked voxels
%     conn        : neighbourhood connectivity (6, 18, or 26)
%     lambda_gmrf : GMRF regularisation strength
%     max_iter    : maximum Newton iterations
%     tol_g, tol_s: gradient and step convergence tolerances
%     jitter      : diagonal jitter for numerical stability
%     pcg_tol     : PCG solver tolerance
%     pcg_maxit   : PCG maximum iterations
%
%   OUTPUT
%     Ymap : [nI x (K-1)] MAP ILR coordinates
%     KAP  : [nI x 1] posterior reliability kappa = tr(Q_v)/(K-1)
%     Fvox : [nI x 1] per-voxel Laplace free energy
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_SUBJECT_VB, HELMERT_SUBMATRIX

[nI, Km1] = size(Y0);
K = Km1 + 1;

% Build sparse Laplacian L on idx
L = build_laplacian_idx(dim3, idx, conn);          % nI x nI, SPD (semi), zero row-sum
% For numerical stability, add tiny diagonal
L = L + 1e-8*speye(nI);

% Vectorize y, m, precision
y = double(reshape(Y0, [], 1));                    % (nI*Km1) x 1
m = double(reshape(M_ilr, [], 1));
P = double(reshape(P_ilr, [], 1));                 % diag precisions per coord

% Convenience indexers
% Map between packed vector and voxel/coord blocks
pack  = @(M) reshape(M, [], 1);
unpack= @(v) reshape(v, [nI Km1]);

% Precompute Helmert once in both directions
H_T = H';

% Newton loop
Fvox = zeros(nI,1,'single');
for it=1:max_iter

    Y = unpack(y);                                  % nI x Km1
    % ---------- per-voxel softmax & stats ----------
    Z  = (H * Y')';                                 % nI x K
    MZ = max(Z, [], 2);
    A  = exp(Z - MZ);
    S  = sum(A,2);
    MU = A ./ S;                                    % nI x K
    MU = max(MU, 1e-12);                            % floors for stability
    MU = bsxfun(@rdivide, MU, sum(MU,2));

    % ---------- block-diagonal Hessians and gradients (likelihood+diag prior) ----------
    % For each voxel v:
    % g_v = H'*(n_v - N_v*mu_v) - diag(P_v)*(y_v - m_v)
    % H_v = diag(P_v) + N_v * H'*(Diag(mu_v) - mu_v*mu_v')*H
    g_vox = zeros(nI, Km1);
    Hdiag_blocks = cell(nI,1);  % store Km1xKm1 blocks (small), used in matvec & precond
    for v = 1:nI
        nv = double(C(v,:))'; Nv = double(Ntot(v));
        yv = Y(v,:)'; mv = double(M_ilr(v,:))';
        Pv = double(P_ilr(v,:))';
        muv = MU(v,:)';
        Sv  = diag(muv) - (muv*muv');
        gv  = H_T*(nv - Nv*muv) - diag(Pv)*(yv - mv);
        Hv  = diag(Pv) + Nv * (H_T * Sv * H) + jitter*eye(Km1);
        g_vox(v,:) = gv';
        Hdiag_blocks{v} = Hv;
    end

    % ---------- GMRF prior terms ----------
    % gradient adds: - lambda * ( (L ⊗ I) y )  ;   Hessian adds: + lambda * (L ⊗ I)
    g_local = pack(g_vox);
    % Apply (L ⊗ I) to y efficiently: for each coord d, multiply L*Y(:,d), then stack
    Ly = apply_L_kron_I(L, Y);                      % nI x Km1
    g_gmrf = -lambda_gmrf * pack(Ly);

    g_all = g_local + g_gmrf;

    % ---------- PCG solve: (H_blockdiag + lambda (L ⊗ I)) * step = g_all ----------
    Afun = @(x) apply_A(x, Hdiag_blocks, lambda_gmrf, L, nI, Km1);   % SPD operator
    Mfun = build_preconditioner(Hdiag_blocks, lambda_gmrf, L, nI, Km1); % simple SPD precond

    [step, flag] = pcg(Afun, g_all, pcg_tol, pcg_maxit, Mfun, [], zeros(size(g_all)));
    if flag~=0
        % fallback: single CG iteration with diagonal preconditioner
        warning('PCG did not fully converge (flag=%d). Proceeding with partial step.', flag);
    end

    % ---------- Backtracking line search on global objective ----------
    F_cur = global_objective(Y, C, Ntot, M_ilr, P_ilr, H, MU, lambda_gmrf, L);
    alpha = 1.0;
    for bt=1:10
        y_try = y + alpha*step;
        Yt    = unpack(y_try);
        Zt    = (H * Yt')';
        Mt    = max(Zt, [], 2);
        At    = exp(Zt - Mt);  St = sum(At,2);
        MUt   = At ./ St; MUt = max(MUt,1e-12); MUt = bsxfun(@rdivide, MUt, sum(MUt,2));
        F_try = global_objective(Yt, C, Ntot, M_ilr, P_ilr, H, MUt, lambda_gmrf, L);
        if F_try >= F_cur + 1e-12
            y = y_try; MU = MUt; break;
        end
        alpha = 0.5*alpha;
        if bt==10
            y = y; % reject
        end
    end

    % ---------- stopping checks ----------
    if norm(step,2) < tol_s*max(1.0, norm(y,2))
        % compute per-voxel precision summaries & per-voxel F (approx)
        [KAP, Fvox] = summarize_precision_and_F(Y, C, Ntot, M_ilr, P_ilr, H); % voxelwise (approx)
        Ymap = single(Y);
        return;
    end
    if norm(g_all,2) < tol_g*max(1.0, norm(y,2))
        [KAP, Fvox] = summarize_precision_and_F(Y, C, Ntot, M_ilr, P_ilr, H);
        Ymap = single(Y);
        return;
    end
end

% Reached max_iter; summarize anyway
Ymap = single(unpack(y));
[KAP, Fvox] = summarize_precision_and_F(unpack(y), C, Ntot, M_ilr, P_ilr, H);
end


function L = build_laplacian_idx(dim3, idx, conn)
% Sparse Laplacian on selected indices of a 3D grid (Neumann-like).
% conn ∈ {6,18,26}
Nx=dim3(1); Ny=dim3(2); Nz=dim3(3);
nI = numel(idx);

% Build a lookup from linear voxel index → compact [1..nI]
map = zeros(Nx*Ny*Nz,1,'uint32');
map(idx) = uint32(1:nI);

% neighbor offsets
offs = [];
if conn==6
    offs = [ -1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1 ];
elseif conn==18
    [X,Y,Z]=ndgrid(-1:1,-1:1,-1:1);
    offs = [X(:) Y(:) Z(:)]; offs = offs(any(abs(offs),2) & sum(abs(offs),2)<=2,:);
else % 26
    [X,Y,Z]=ndgrid(-1:1,-1:1,-1:1);
    offs = [X(:) Y(:) Z(:)]; offs = offs(any(abs(offs),2),:);
end

[I,J,V] = deal([],[],[]);
deg = zeros(nI,1);
[xi,yi,zi] = ind2sub([Nx Ny Nz], idx);
for t=1:size(offs,1)
    xn = xi + offs(t,1); yn = yi + offs(t,2); zn = zi + offs(t,3);
    valid = xn>=1 & xn<=Nx & yn>=1 & yn<=Ny & zn>=1 & zn<=Nz;
    nei_lin = sub2ind([Nx Ny Nz], xn(valid), yn(valid), zn(valid));
    src = map(idx(valid)); dst = map(nei_lin);
    keep = (dst>0); src = src(keep); dst = dst(keep);
    I = [I; double(src)]; J = [J; double(dst)]; V = [V; -ones(nnz(keep),1)];
    deg(src) = deg(src) + 1;
end

% Diagonal entries = degree
I = [I; (1:nI)']; J = [J; (1:nI)']; V = [V; deg];
L = sparse(I,J,V,nI,nI);
end


function Ly = apply_L_kron_I(L, Y)
% Y: nI x Km1  →  Ly: nI x Km1, each column L*Y(:,d)
Ly = L * Y;  % sparse * dense
end

function Ax = apply_A(x, Hdiag_blocks, lambda, L, nI, Km1)
% A = blockdiag(Hv) + lambda*(L ⊗ I)
% x, Ax are packed (nI*Km1 x 1)
X = reshape(x, [nI Km1]);

% blockdiag(Hv)*X  (small Km1xKm1 blocks)
BX = zeros(nI, Km1);
for v=1:nI
    BX(v,:) = (Hdiag_blocks{v} * X(v,:)')';
end

% lambda*(L ⊗ I)*X  → apply L to each column
LX = apply_L_kron_I(L, X);
AX = BX + lambda * LX;
Ax = AX(:);
end

function Mfun = build_preconditioner(Hdiag_blocks, lambda, L, nI, Km1)
% Simple block-Jacobi + Laplacian diagonal preconditioner:
% M ≈ blockdiag(Hv + lambda*deg(v) I)
deg = full(sum(L~=0,2));  % degree per voxel
D = cell(nI,1);
for v=1:nI
    D{v} = Hdiag_blocks{v} + lambda * deg(v) * eye(Km1);
    % Precompute Cholesky of each small Km1xKm1 block for fast solves
    [R,p] = chol((D{v}+D{v}')/2,'upper'); %#ok<NASGU>
    if p==0
        D{v} = R;         % store R (upper)
    else
        % fallback: add tiny jitter
        [R,~] = chol((D{v}+D{v}')/2 + 1e-6*eye(Km1),'upper');
        D{v} = R;
    end
end
Mfun = @(r) precond_apply(r, D, nI, Km1);
end

function z = precond_apply(r, D, nI, Km1)
R = reshape(r, [nI Km1]);
Z = zeros(nI, Km1);
for v=1:nI
    y = R(v,:)';
    % Solve (R'R) z = y  → two triangular solves
    z1 = D{v}' \ y;
    z2 = D{v}  \ z1;
    Z(v,:) = z2';
end
z = Z(:);
end

function F = global_objective(Y, C, Ntot, M_ilr, P_ilr, H, MU, lambda, L)
% Sum of voxelwise log-lik + log-prior  - (lambda/2)*y'(L⊗I)y
[nI, Km1] = size(Y); K = Km1+1;
Z = (H * Y')';               % nI x K
MZ = max(Z,[],2);
A  = exp(Z - MZ);
S  = sum(A,2);
lse= log(S) + MZ;

ll = 0; lp = 0;
for v=1:nI
    nv = double(C(v,:))';  Nv = double(Ntot(v));
    zv = Z(v,:)';          yv = Y(v,:)';  mv = double(M_ilr(v,:))';  Pv = double(P_ilr(v,:))';
    ll = ll + (nv' * zv - Nv * lse(v));
    lp = lp - 0.5*(yv-mv)'*(diag(Pv)*(yv-mv)) + 0.5*sum(log(Pv));
end
% GMRF quadratic
Qsp = 0;
for d=1:Km1
    qd = Y(:,d)' * (L * Y(:,d));
    Qsp = Qsp + qd;
end
F = ll + lp - 0.5*lambda * Qsp;
end

function [KAP, Fvox] = summarize_precision_and_F(Y, C, Ntot, M_ilr, P_ilr, H)
% Per-voxel summaries using block-diagonal Hessians only (standard voxelwise Laplace).
% This keeps backward compatibility with your existing ELBO bookkeeping.
[nI, Km1] = size(Y); K = Km1+1;
KAP  = zeros(nI,1,'single');
Fvox = zeros(nI,1,'single');
for v=1:nI
    yv = double(Y(v,:)'); mv = double(M_ilr(v,:))'; Pv = double(P_ilr(v,:))';
    zv = H*yv; [lse, muv] = logsumexp_softmax(zv);
    Nv = double(Ntot(v)); nv = double(C(v,:))';
    Sv = diag(muv) - (muv*muv');
    Qv = diag(Pv) + Nv*(H' * Sv * H) + 1e-10*eye(Km1);
    [ldQ,~] = logdet_pd(Qv);
    like_const = gammaln(Nv+1) - sum(gammaln(nv+1));
    ll = nv' * zv - Nv * lse;
    lp = -0.5*(yv-mv)'*(diag(Pv)*(yv-mv)) + 0.5*sum(log(Pv));
    Fv = ll + like_const + lp - 0.5*ldQ;   % voxelwise Laplace
    Fvox(v) = single(Fv);
    KAP(v)  = single(trace(Qv)/Km1);
end
end
