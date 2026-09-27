"""Spatially varying covariate-conditioned prior reference: simulation evidence.

An unconditional group prior assumes the same expected connectivity field for every
subject.  When a known characteristic is associated with systematic voxelwise
variation, that assumption leaves a spatially structured mismatch.  The predictive
field delta_v^(s) = (beta_v^pred)' x^(s) is meant to correct it, and because
beta_v^pred varies over voxels the correction is spatially heterogeneous -- the same
subject can need a positive shift at one voxel and a negative shift at another.

Generative model (ILR coordinates, K=68 targets):
    y_v^(s) = y_v^base + beta_v^true x_s + eps,   n_v^(s) ~ Multinomial(N_v, softmax(H y))
N_v spans a wide range so that r_data = tr(Q_data)/tr(Q_total) covers [0,1].

Held-out subjects are scored under two otherwise identical inference models -- same
Lambda_grp, same alpha_grp, same counts, same Newton MAP -- differing only in the
prior mean, with mu^grp and beta_hat^pred estimated ONLY from the training subjects:

    unconditioned   m_v^(s) = mu_v^grp
    conditioned     m_v^(s) = mu_v^grp + beta_hat_v^pred x_s

Saved outputs drive a four-panel supplementary figure:
  A  the spatial offset field ||delta_v^(s)|| for a younger and an older subject
  B  prior miscentering ||m_v - y_v,true|| for both priors, and its difference
  C  posterior composition error against the known subject-specific truth
  D  the gain against r_data, stratified by the size of the offset
plus one voxelwise inset (pi_grp, pi_cond, pi_true) at a large-offset, low-r_data voxel.
Run on server 89.
"""
import os, json
os.environ['OMP_NUM_THREADS'] = '16'; os.environ['OPENBLAS_NUM_THREADS'] = '16'
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K, K1 = 68, 67
V = 600
S_TRAIN, S_TEST = 50, 12
ALPHA_GRP = 0.25
N_REP = 4
NEWTON_IT = 15

def ilr_basis(K):
    Hm = np.zeros((K, K - 1))
    for j in range(K - 1):
        c = np.sqrt((j + 1) / (j + 2))
        Hm[:j + 1, j] = c / (j + 1)
        Hm[j + 1, j] = -c
    return Hm
H = ilr_basis(K)

def softmax_ilr(Y):
    Z = Y @ H.T
    Z -= Z.max(-1, keepdims=True)
    E = np.exp(Z)
    return E / E.sum(-1, keepdims=True)

def _curv(P, Nv):
    """Q_data = N * H'(diag(P) - PP')H, batched over voxels."""
    A = P[:, :, None] * H[None]                 # (V,K,K1)
    H2 = np.matmul(A.transpose(0, 2, 1), np.broadcast_to(H, (P.shape[0], K, K1)))
    PH = P @ H                                  # (V,K1)
    return Nv[:, None, None] * (H2 - PH[:, :, None] * PH[:, None, :])

def map_newton(cnt, Nv, m, lam):
    y = m.copy()
    for _ in range(NEWTON_IT):
        P = softmax_ilr(y)
        g = (cnt - Nv[:, None] * P) @ H - lam * (y - m)
        Q = _curv(P, Nv)
        Q[:, np.arange(K1), np.arange(K1)] += lam
        step = np.linalg.solve(Q, g[:, :, None])[:, :, 0]
        y = y + step
        if np.max(np.abs(step)) < 1e-8:
            break
    P = softmax_ilr(y)
    trd = np.trace(_curv(P, Nv), axis1=1, axis2=2)
    return y, P, trd / (trd + lam * K1)

def naive_ilr(c):
    comp = (c + 0.5) / (c.sum(-1, keepdims=True) + 0.5 * K)
    return np.log(comp) @ H

def one_rep(rep, bscale, keep=False):
    g = np.random.default_rng(7000 + 97 * rep + int(1000 * bscale))
    S = S_TRAIN + S_TEST
    ybase = g.normal(size=(V, K1)) * 1.5
    # ||beta_true|| varies over voxels; bscale sets its scale relative to the
    # between-subject residual norm SIGMA_E*sqrt(K1), so the covariate share of
    # between-subject variance is known and is reported alongside each row.
    mag = np.abs(g.normal(size=V)) * bscale
    btrue = g.normal(size=(V, K1))
    btrue *= (mag / np.maximum(np.linalg.norm(btrue, axis=1), 1e-12))[:, None]
    x = g.normal(size=S); x = (x - x.mean()) / x.std()
    Ytrue = ybase[None] + x[:, None, None] * btrue[None] + g.normal(size=(S, V, K1)) * 0.6
    Nv = np.exp(g.uniform(np.log(20), np.log(5000), size=V))
    CNT = np.empty((S, V, K))
    for s in range(S):
        P = softmax_ilr(Ytrue[s])
        CNT[s] = np.stack([g.multinomial(int(Nv[v]), P[v]) for v in range(V)])
    tr = np.arange(S_TRAIN); te = np.arange(S_TRAIN, S)
    Ytr = naive_ilr(CNT[tr])
    Xtr = np.column_stack([np.ones(S_TRAIN), x[tr]])
    B = np.linalg.pinv(Xtr) @ Ytr.reshape(S_TRAIN, V * K1)
    mu_grp = B[0].reshape(V, K1); bhat = B[1].reshape(V, K1)
    lam = ALPHA_GRP / max(np.var(Ytr.reshape(S_TRAIN, -1), axis=0).mean(), 1e-9)

    n = len(te)
    priU = np.zeros((n, V)); priC = np.zeros((n, V))
    posU = np.zeros((n, V)); posC = np.zeros((n, V))
    rd = np.zeros((n, V)); off = np.zeros((n, V))
    store = None
    for i, s in enumerate(te):
        dm = bhat * x[s]
        off[i] = np.linalg.norm(dm, axis=1)
        priU[i] = np.linalg.norm(mu_grp - Ytrue[s], axis=1)
        priC[i] = np.linalg.norm(mu_grp + dm - Ytrue[s], axis=1)
        Ptrue = softmax_ilr(Ytrue[s])
        _, Pu, r = map_newton(CNT[s], Nv, mu_grp, lam)
        _, Pc, _ = map_newton(CNT[s], Nv, mu_grp + dm, lam)
        posU[i] = np.abs(Pu - Ptrue).mean(1)
        posC[i] = np.abs(Pc - Ptrue).mean(1)
        rd[i] = r
        if keep and i == 0:
            store = dict(mu_grp=mu_grp, bhat=bhat, x_te=x[te], Ptrue=Ptrue,
                         Pu=Pu, Pc=Pc, Pgrp=softmax_ilr(mu_grp),
                         Pcond=softmax_ilr(mu_grp + dm), rd=r, off=off[i],
                         btrue_mag=np.linalg.norm(btrue, axis=1), Nv=Nv)
    return dict(priU=priU, priC=priC, posU=posU, posC=posC, rd=rd, off=off,
                btrue_mag=np.tile(np.linalg.norm(btrue, axis=1), (n, 1)), store=store)

res = {}
print(f'{"beta_scale":>10s} {"prior unc":>10s} {"prior cond":>11s} {"prior gain%":>12s} '
      f'{"post unc":>10s} {"post cond":>10s} {"post gain%":>11s} {"cond better":>12s}', flush=True)
for bs in (0.0, 0.5, 1.5, 3.0, 6.0):
    acc = {k: [] for k in ('priU', 'priC', 'posU', 'posC', 'rd', 'off', 'btrue_mag')}
    for rep in range(N_REP):
        o = one_rep(rep, bs, keep=(rep == 0 and abs(bs - 3.0) < 1e-9))
        for k in acc: acc[k].append(o[k])
        if o['store'] is not None:
            np.savez(f'{HERE}/predictive_prior_sim_panels.npz', **o['store'])
    A = {k: np.concatenate(v).ravel() for k, v in acc.items()}
    pg = 100 * (A['priU'].mean() - A['priC'].mean()) / A['priU'].mean()
    sg = 100 * (A['posU'].mean() - A['posC'].mean()) / A['posU'].mean()
    print(f'{bs:10.2f} {A["priU"].mean():10.4f} {A["priC"].mean():11.4f} {pg:12.3f} '
          f'{A["posU"].mean():10.5f} {A["posC"].mean():10.5f} {sg:11.3f} '
          f'{100*np.mean(A["posC"] < A["posU"]):11.1f}%', flush=True)
    bt = A['btrue_mag']
    r2 = float(np.mean(bt ** 2 / (bt ** 2 + (0.6 ** 2) * K1)))
    print(f'            implied covariate share of between-subject variance R^2 = {r2:.4f}', flush=True)
    ent = dict(beta_scale=bs, r2=r2, prior_unc=float(A['priU'].mean()), prior_cond=float(A['priC'].mean()),
               prior_gain_pct=float(pg), post_unc=float(A['posU'].mean()),
               post_cond=float(A['posC'].mean()), post_gain_pct=float(sg),
               frac_better=float(np.mean(A['posC'] < A['posU'])))
    if bs > 0:
        qr = np.quantile(A['rd'], [0.25, 0.5, 0.75]); qo = np.quantile(A['off'], [0.5])
        grid = []
        qb = np.quantile(A['btrue_mag'], [0.5])
        for om, oname in [(A['btrue_mag'] <= qb[0], 'small true offset'), (A['btrue_mag'] > qb[0], 'large true offset')]:
            row = []
            for rm in [A['rd'] <= qr[0], (A['rd'] > qr[0]) & (A['rd'] <= qr[1]),
                       (A['rd'] > qr[1]) & (A['rd'] <= qr[2]), A['rd'] > qr[2]]:
                m = om & rm
                row.append(float(100 * (A['posU'][m].mean() - A['posC'][m].mean()) / A['posU'][m].mean())
                           if m.sum() else float('nan'))
            grid.append(row); 
        ent['gain_grid'] = grid; ent['rdata_cuts'] = [float(v) for v in qr]
        ent['offset_median'] = float(qo[0])
        print('            posterior gain %, rows = TRUE offset (small,large), cols = r_data quartile:', flush=True)
        for row in grid:
            print('             ' + '  '.join(f'{v:7.3f}' for v in row), flush=True)
    res[f'bs{bs}'] = ent
json.dump(res, open(f'{HERE}/predictive_prior_sim.json', 'w'), indent=1)
print('\n=== DONE ===', flush=True)
