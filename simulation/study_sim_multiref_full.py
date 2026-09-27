#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-reference sensitivity analysis for iCDM simulation.

For each of N_REF HBN subjects used as the ground-truth reference anatomy,
runs the full S=20 simulation and reports key statistics.

Used for the ten-reference results (manuscript Table 2; Supplementary Fig. S2(j)-(l)).

Output:
  figures/multiref_results.csv     — per-reference statistics
  figures/multiref_summary.txt     — mean ± SD across references
"""

import numpy as np
import scipy.io as sio
from scipy.special import rel_entr
from scipy.stats import wilcoxon
import os, csv, multiprocessing as mp

try:
    import nibabel as nib
    HAS_NIB = True
except ImportError:
    HAS_NIB = False
    print("[WARN] nibabel not found — will skip NIfTI-based references")

# ============================================================
# CONFIG
# ============================================================
HBN_BASE   = os.environ.get("ICDM_HBN_BASE", "path/to/HBN/SC")   # folder with per-subject count maps (set ICDM_HBN_BASE)
NII_FNAME  = "WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii.gz"
CORONAL_Y  = 74

REF_SUBJECTS = [
    "NDARAT100AEQ",   # original reference
    "NDARAL828WXM",
    "NDARAM277WZT",
    "NDARAV610EY3",
    "NDARAX283MAK",
    "NDARBF998MBA",
    "NDARBN100LCD",
    "NDARBT436PMT",
    "NDARBU928LV0",
    "NDARCD182XT1",
]

S            = 20
SIGMA_SUBJ   = 0.15
GAMMA        = 0.35
ALPHA_DM     = 1.0
N_EB_ITER    = 1
SMOOTH_SIGMA = 1.0
LAMBDA_SCALE = 25.0
ALPHA_BLEND  = 0.5
TRIM_ALPHA   = 0.2
TAU_FLOOR    = 0.05
PREC_SCALE   = 1.0

OUTDIR = "figures/R1.1_multiref_full"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Helper functions
# ============================================================
def helmert_basis(K):
    # Orthonormal Helmert ILR basis (Egozcue): H^T H = I, H^T 1 = 0.
    # Matches MATLAB helmert_submatrix.m and pyicdm/helmert.py.
    H = np.zeros((K, K - 1))
    for j in range(K - 1):
        i = j + 1
        H[:i, j] = 1.0 / np.sqrt(i * (i + 1))
        H[i, j]  = -np.sqrt(i / (i + 1))
    assert np.abs(H.T @ H - np.eye(K - 1)).max() < 1e-10, "Helmert basis not orthonormal"
    return H

def clr_forward(pi):
    pi = np.clip(pi, 1e-12, 1)
    return np.log(pi) - np.log(pi).mean(axis=1, keepdims=True)

def clr_inverse(z):
    """Inverse CLR: exp then close to simplex.  (verbatim from study_sim_fig4.py)"""
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)

def ilr_forward(pi, H):
    return clr_forward(pi) @ H

def ilr_inverse_softmax(y, H):
    z = y @ H.T
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)

def _softmax_pi(y, H):
    """Vectorized softmax: (V,D) -> (V,K)."""
    z = y @ H.T
    z -= z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)

def icdm_map_voxelwise_vec(counts, H, prior_mean, prior_prec,
                            gamma_init=0.35, alpha_blend=0.5,
                            max_iter=30, tol=1e-6):
    """
    Vectorized iCDM MAP inference (all voxels processed as a batch).
    ~100× faster than the per-voxel Python loop for large V.
    """
    V, K = counts.shape
    D = K - 1
    Nv = counts.sum(axis=1).astype(float)

    if prior_prec.ndim == 1:
        lam = np.tile(prior_prec[:, None], (1, D))  # (V, D)
    else:
        lam = prior_prec.copy()

    # Initialization
    ct     = (counts + 1.0) ** gamma_init
    pi_ini = ct / ct.sum(axis=1, keepdims=True)
    y_ini  = ilr_forward(pi_ini, H)              # (V, D)
    y      = (1 - alpha_blend) * prior_mean + alpha_blend * y_ini

    valid = (Nv >= 20)
    active = valid.copy()

    for _ in range(max_iter):
        if not active.any():
            break
        idx = np.where(active)[0]
        y_a = y[idx];   m_a = prior_mean[idx]
        lam_a = lam[idx]; n_a = counts[idx]
        N_a = Nv[idx]

        # pi = softmax(H y)
        piv = _softmax_pi(y_a, H)               # (A, K)

        # Gradient: H^T (n - N pi) - lam*(y-m)  => vectorized as (n-Npi)@H
        grad = (n_a - N_a[:, None] * piv) @ H - lam_a * (y_a - m_a)  # (A, D)

        # Hessian: Q = diag(lam) + N * H^T S_pi H
        # H^T S_pi H = (sqrt(pi)*H)^T (sqrt(pi)*H) - (pi@H)(pi@H)^T
        sqH = np.sqrt(piv)[:, :, None] * H[None, :, :]        # (A, K, D)
        HpH = np.matmul(sqH.transpose(0, 2, 1), sqH)         # (A, D, D)  ~16× faster than einsum
        Hp  = piv @ H                                          # (A, D)
        HpH -= Hp[:, :, None] * Hp[:, None, :]

        Q = lam_a[:, :, None] * np.eye(D)[None, :, :] + N_a[:, None, None] * HpH

        # Batch solve
        dy = np.linalg.solve(Q, grad[:, :, None]).squeeze(-1)  # (A, D)
        y[idx] += dy

        norms = np.linalg.norm(dy, axis=1)
        active[idx[norms < tol]] = False

    # Set invalid voxels to prior mean (or init)
    y[~valid] = prior_mean[~valid]

    pi_hat = _softmax_pi(y, H)  # (V, K)

    # κ = tr(Q_final) / D  for valid voxels only
    kappa = np.zeros(V)
    if valid.any():
        vi = np.where(valid)[0]
        pf = pi_hat[vi]
        sqHf = np.sqrt(pf)[:, :, None] * H[None, :, :]
        HpHf = np.matmul(sqHf.transpose(0, 2, 1), sqHf)
        Hpf  = pf @ H
        HpHf -= Hpf[:, :, None] * Hpf[:, None, :]
        Qf   = lam[vi, :, None] * np.eye(D)[None, :, :] + Nv[vi, None, None] * HpHf
        kappa[vi] = np.einsum('aii->a', Qf) / D   # batch trace

    return pi_hat, y, kappa


def trimmed_mean(arr, axis=0, alpha=0.2):
    n = arr.shape[axis]; lo = int(np.floor(alpha * n)); hi = n - lo
    return np.mean(np.sort(arr, axis=axis)[lo:hi], axis=axis)


def mad_std(arr, axis=0):
    med = np.median(arr, axis=axis, keepdims=True)
    return 1.4826 * np.median(np.abs(arr - med), axis=axis)


def spatial_smooth_2d(field_flat, Hdim, Wdim, maskv, sigma=1.5):
    from scipy.ndimage import gaussian_filter
    out = field_flat.copy()
    wt  = maskv.reshape(Hdim, Wdim).astype(float)
    den_precomp = gaussian_filter(wt, sigma=sigma)
    for ch in range(field_flat.shape[1]):
        img = field_flat[:, ch].reshape(Hdim, Wdim)
        num = gaussian_filter(np.where(wt > 0, img, 0.0), sigma=sigma)
        out[:, ch] = (num / np.maximum(den_precomp, 1e-12)).ravel()
    return out


def js_divergence(p, q):
    m = 0.5 * (p + q)
    return 0.5 * (np.sum(rel_entr(p, m), axis=1) + np.sum(rel_entr(q, m), axis=1))

def per_subject_mae(pi_list, pi_ref, mask):
    return np.mean([np.abs(p[mask] - pi_ref[mask]).mean(axis=1) for p in pi_list], axis=0)

def per_subject_mae_paired(pi_list, ref_list, mask):
    return np.mean([np.abs(p[mask] - r[mask]).mean(axis=1)
                    for p, r in zip(pi_list, ref_list)], axis=0)

def per_subject_js_paired(pi_list, ref_list, mask):
    return np.mean([js_divergence(p[mask], r[mask])
                    for p, r in zip(pi_list, ref_list)], axis=0)

def per_subject_js(pi_list, pi_ref, mask):
    return np.mean([js_divergence(p[mask], pi_ref[mask]) for p in pi_list], axis=0)

# ============================================================
# Load reference anatomy
# ============================================================
def load_ref_from_nii(subj_id):
    path = os.path.join(HBN_BASE, subj_id, NII_FNAME)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    d = nib.load(path).get_fdata(dtype=np.float32)   # (X, Y, Z, K)
    return d[:, CORONAL_Y, :, :68].astype(np.float64) # (128, 96, 68)

# ============================================================
# Core simulation (one reference anatomy)
# ============================================================
def run_simulation(args):
    """Worker function: args = (subj_id, C_mat, seed)."""
    subj_id, C_mat, seed = args
    import numpy as np

    K = 68
    C = C_mat[:, :, :K].astype(float)
    Hdim, Wdim = C.shape[:2]
    V = Hdim * Wdim

    maskv = (C.sum(axis=2) >= 20).ravel()
    Ct    = (C + 1.0) ** GAMMA
    pi_true = Ct / Ct.sum(axis=2, keepdims=True)
    pi_true = pi_true.reshape(V, K)

    H_ilr = helmert_basis(K)
    D     = K - 1
    y_true = ilr_forward(pi_true, H_ilr)

    np.random.seed(seed)
    Nv_real = C.reshape(V, K).sum(axis=1).astype(int)

    counts_all = []
    pi_subj_true = []
    for s in range(S):
        y_s  = y_true + SIGMA_SUBJ * np.random.randn(V, D)
        pi_s = ilr_inverse_softmax(y_s, H_ilr)
        pi_subj_true.append(pi_s)
        n_s  = np.zeros((V, K), dtype=float)
        for v in range(V):
            if maskv[v] and Nv_real[v] >= 20:
                n_s[v] = np.random.multinomial(Nv_real[v], pi_s[v])
        counts_all.append(n_s)

    # Baselines
    pi_naive_all, pi_dm_all, pi_clr_all = [], [], []
    for s in range(S):
        c = counts_all[s]; Ns = c.sum(axis=1, keepdims=True)
        pi_nv = c / np.maximum(Ns, 1)
        pi_naive_all.append(pi_nv)
        pi_dm_all.append((c + ALPHA_DM) / (Ns + K * ALPHA_DM))
        # CLR baseline: sigma=1.5 as in study_sim_fig4.py, NOT SMOOTH_SIGMA
        clr_c = clr_forward(pi_nv)
        clr_s = spatial_smooth_2d(clr_c, Hdim, Wdim, maskv, sigma=1.5)
        pi_clr_all.append(clr_inverse(clr_s))

    # iCDM-single: per-subject MAP against a DM-derived prior, no cross-subject pooling
    pi_single_all = []
    for s in range(S):
        c = counts_all[s]
        Ns_arr = c.sum(axis=1).astype(float)
        med_N_s = float(np.median(Ns_arr[maskv]))
        prior_mean_s = ilr_forward(pi_dm_all[s], H_ilr)
        prior_prec_s = np.where(maskv, LAMBDA_SCALE * med_N_s / np.maximum(Ns_arr, 1), 1.0)
        pi_s1, _, _ = icdm_map_voxelwise_vec(
            c, H_ilr, prior_mean_s, prior_prec_s,
            gamma_init=GAMMA, alpha_blend=ALPHA_BLEND, max_iter=30, tol=1e-6)
        pi_single_all.append(pi_s1)

    # iCDM-EB
    Ns_arr_eb = np.stack([counts_all[s].sum(axis=1) for s in range(S)])
    med_N     = float(np.median(Ns_arr_eb[:, maskv]))
    y_hat_all  = [ilr_forward(pi_dm_all[s], H_ilr) for s in range(S)]
    kappa_all  = [np.ones(V) for _ in range(S)]

    for eb_iter in range(N_EB_ITER):
        y_stack = np.stack(y_hat_all)
        mu_grp  = trimmed_mean(y_stack, axis=0, alpha=TRIM_ALPHA)
        tau_grp = np.maximum(mad_std(y_stack, axis=0), TAU_FLOOR)
        mu_grp  = spatial_smooth_2d(mu_grp, Hdim, Wdim, maskv, sigma=SMOOTH_SIGMA)
        # Group-prior precision: inverse between-subject dispersion scaled by one global constant
        # (the simulation counterpart of the population-transfer strength alpha_grp).
        lam_grp = 1.0 / (tau_grp ** 2)
        Lam_grp = LAMBDA_SCALE * med_N * lam_grp / np.maximum(lam_grp.mean(), 1e-6)

        y_hat_all, kappa_all, pi_prop_all = [], [], []
        for s in range(S):
            pi_s, y_s, kappa_s = icdm_map_voxelwise_vec(
                counts_all[s], H_ilr, mu_grp, Lam_grp,
                gamma_init=GAMMA, alpha_blend=ALPHA_BLEND,
                max_iter=30, tol=1e-6)
            y_hat_all.append(y_s); kappa_all.append(kappa_s); pi_prop_all.append(pi_s)

    mae_naive = per_subject_mae(pi_naive_all, pi_true, maskv)
    mae_dm    = per_subject_mae(pi_dm_all,    pi_true, maskv)
    mae_eb    = per_subject_mae(pi_prop_all,  pi_true, maskv)
    mae_clr    = per_subject_mae(pi_clr_all,    pi_true, maskv)
    mae_single = per_subject_mae(pi_single_all, pi_true, maskv)

    # group-only: the EB population mean applied unchanged to every subject
    pi_group     = ilr_inverse_softmax(mu_grp, H_ilr)
    pi_group_all = [pi_group for _ in range(S)]
    mae_group    = per_subject_mae(pi_group_all, pi_true, maskv)

    _EST = {"naive": pi_naive_all, "dm": pi_dm_all, "clr": pi_clr_all,
            "single": pi_single_all, "group": pi_group_all, "eb": pi_prop_all}

    _js, _maes, _jss, _rho = {}, {}, {}, {}
    from scipy.stats import spearmanr as _spr
    for _n, _lst in _EST.items():
        _js[_n]   = per_subject_js(_lst, pi_true, maskv).mean()
        _maes[_n] = per_subject_mae_paired(_lst, pi_subj_true, maskv).mean()
        _jss[_n]  = per_subject_js_paired(_lst, pi_subj_true, maskv).mean()
        _m = np.mean(np.asarray(_lst), axis=0)
        _rho[_n] = float(np.nanmean([_spr(pi_true[maskv, _k], _m[maskv, _k])[0]
                                     for _k in range(pi_true.shape[1])]))

    pct_vs_naive = (mae_naive.mean() - mae_eb.mean()) / mae_naive.mean() * 100
    pct_vs_dm    = (mae_dm.mean()    - mae_eb.mean()) / mae_dm.mean()    * 100
    diff         = mae_dm - mae_eb
    cohen_d      = diff.mean() / diff.std()
    win_rate     = (mae_eb < mae_dm).mean() * 100
    _, p_w       = wilcoxon(mae_dm, mae_eb, alternative='greater')

    return {
        "subject":      subj_id,
        "n_vox":        int(maskv.sum()),
        "mae_naive":    mae_naive.mean(),
        "mae_dm":       mae_dm.mean(),
        "mae_eb":       mae_eb.mean(),
        "mae_clr":      mae_clr.mean(),
        "mae_single":   mae_single.mean(),
        "mae_group":    mae_group.mean(),
        "js_naive":      _js["naive"],
        "js_dm":         _js["dm"],
        "js_clr":        _js["clr"],
        "js_single":     _js["single"],
        "js_group":      _js["group"],
        "js_eb":         _js["eb"],
        "maeS_naive":    _maes["naive"],
        "maeS_dm":       _maes["dm"],
        "maeS_clr":      _maes["clr"],
        "maeS_single":   _maes["single"],
        "maeS_group":    _maes["group"],
        "maeS_eb":       _maes["eb"],
        "jsS_naive":     _jss["naive"],
        "jsS_dm":        _jss["dm"],
        "jsS_clr":       _jss["clr"],
        "jsS_single":    _jss["single"],
        "jsS_group":     _jss["group"],
        "jsS_eb":        _jss["eb"],
        "rho_naive":     _rho["naive"],
        "rho_dm":        _rho["dm"],
        "rho_clr":       _rho["clr"],
        "rho_single":    _rho["single"],
        "rho_group":     _rho["group"],
        "rho_eb":        _rho["eb"],
        "pct_vs_naive": pct_vs_naive,
        "pct_vs_dm":    pct_vs_dm,
        "cohen_d":      cohen_d,
        "win_rate":     win_rate,
        "p_wilcox":     p_w,
    }

# ============================================================
# Main
# ============================================================
def main():
    import time

    jobs = []
    for i, subj in enumerate(REF_SUBJECTS):
        try:
            C_mat = load_ref_from_nii(subj) if HAS_NIB else None
            if C_mat is None:
                raise RuntimeError("nibabel unavailable")
        except Exception as e:
            if subj == "NDARAT100AEQ":
                mat   = sio.loadmat("icdm84.mat")
                C_mat = mat["C"][:, :, :68].astype(float)
                print(f"  [fallback to icdm84.mat for {subj}]")
            else:
                print(f"  [SKIP {subj}]: {e}")
                continue
        jobs.append((subj, C_mat, 42 + i))

    print(f"Running {len(jobs)} references (serial, ETA ~{len(jobs)*4} min)")
    t0 = time.time()
    results = []
    for j in jobs:
        print(f"  Starting {j[0]}...", flush=True)
        r = run_simulation(j)
        results.append(r)
        print(f"  {j[0]}: vs_DM={r['pct_vs_dm']:.1f}%, d={r['cohen_d']:.3f}", flush=True)
    t1 = time.time()

    results = [r for r in results if r is not None]
    print(f"Done in {t1-t0:.0f}s ({len(results)} references)")

    # Save CSV
    _M = ["naive","dm","clr","single","group","eb"]
    fields = (["subject","n_vox","mae_naive","mae_dm","mae_clr","mae_single",
               "mae_group","mae_eb"]
              + [f"js_{n}" for n in _M] + [f"maeS_{n}" for n in _M]
              + [f"jsS_{n}" for n in _M] + [f"rho_{n}" for n in _M]
              + ["pct_vs_naive","pct_vs_dm","cohen_d","win_rate","p_wilcox"])
    with open(f"{OUTDIR}/multiref_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results: w.writerow({k: r[k] for k in fields})

    pct_dm    = np.array([r["pct_vs_dm"]    for r in results])
    pct_naive = np.array([r["pct_vs_naive"] for r in results])
    d_vals    = np.array([r["cohen_d"]      for r in results])
    win_vals  = np.array([r["win_rate"]     for r in results])

    summary = (
        f"\n{'='*60}\n"
        f"MULTI-REFERENCE SENSITIVITY ({len(results)} subjects)\n"
        f"{'='*60}\n"
        f"  MAE vs DM:    {pct_dm.mean():.1f}% ± {pct_dm.std():.1f}%  [{pct_dm.min():.1f}–{pct_dm.max():.1f}%]\n"
        f"  MAE vs Naive: {pct_naive.mean():.1f}% ± {pct_naive.std():.1f}%  [{pct_naive.min():.1f}–{pct_naive.max():.1f}%]\n"
        f"  Cohen's d:    {d_vals.mean():.2f} ± {d_vals.std():.2f}  [{d_vals.min():.2f}–{d_vals.max():.2f}]\n"
        f"  Win rate:     {win_vals.mean():.1f}% ± {win_vals.std():.1f}%\n"
        f"{'='*60}\n"
    )
    print(summary)
    for r in results:
        print(f"  {r['subject']}: vs_DM={r['pct_vs_dm']:.1f}%, d={r['cohen_d']:.3f}, win={r['win_rate']:.1f}%")

    with open(f"{OUTDIR}/multiref_summary.txt", "w") as f:
        f.write(summary)
        for r in results:
            f.write(f"  {r['subject']}: vs_DM={r['pct_vs_dm']:.1f}%, "
                    f"d={r['cohen_d']:.3f}, win={r['win_rate']:.1f}%\n")
    print(f"\n[SAVED] {OUTDIR}/multiref_results.csv, multiref_summary.txt")


if __name__ == "__main__":
    main()
