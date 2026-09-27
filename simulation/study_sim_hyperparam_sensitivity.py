#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-at-a-time hyperparameter sensitivity analysis for iCDM.

Hyperparameters examined (Supplementary Table S2):
  trimming fraction (TRIM_ALPHA)
  group-prior spatial smoothing (SMOOTH_SIGMA, proxy for lambda_sigma)
  prior precision scale (LAMBDA_SCALE, proxy for regularization strength)

Each parameter is swept independently while all others are held at default.
Defaults: TRIM_ALPHA=0.20, SMOOTH_SIGMA=1.0, LAMBDA_SCALE=25.0

gamma=0.35, N_EB_ITER=1, ALPHA_BLEND=0.5 are fixed throughout.

Outputs:
  figures/R3.5_hyperparam_sensitivity/hyperparam_results.csv
  figures/R3.5_hyperparam_sensitivity/hyperparam_summary.txt
"""

import numpy as np
import scipy.io as sio
from scipy.stats import wilcoxon
import os, csv, time

from study_sim_multiref import (
    helmert_basis, ilr_forward, icdm_map_voxelwise_vec,
    trimmed_mean, mad_std, spatial_smooth_2d,
    ALPHA_DM, TAU_FLOOR,
)

# ============================================================
# CONFIG — defaults
# ============================================================
MAT_PATH     = "icdm84.mat"
OUTDIR       = "figures/R3.5_hyperparam_sensitivity"
os.makedirs(OUTDIR, exist_ok=True)

GAMMA        = 0.35
N_EB_ITER    = 1
ALPHA_BLEND  = 0.5
S            = 20
SIGMA_SUBJ   = 0.15
SEEDS        = [42, 43, 44]

# Defaults
DEF_TRIM_ALPHA   = 0.20
DEF_SMOOTH_SIGMA = 1.0
DEF_LAMBDA_SCALE = 25.0

# Sweep grids
TRIM_VALUES   = [0.00, 0.05, 0.10, 0.20]   # 0 = no trimming
SMOOTH_VALUES = [0.0,  0.5,  1.0,  2.0]    # 0 = no smoothing
LAMBDA_VALUES = [5.0,  10.0, 25.0, 50.0, 100.0]

# ============================================================
# Load reference anatomy (same as gamma sensitivity)
# ============================================================
mat   = sio.loadmat(MAT_PATH)
C_ref = mat["C"][:, :, :68].astype(float)
Hdim, Wdim, K = C_ref.shape
V  = Hdim * Wdim
D  = K - 1
H  = helmert_basis(K)

maskv   = (C_ref.sum(axis=2) >= 20).ravel()
Nv_real = C_ref.reshape(V, K).sum(axis=1).astype(int)

# Ground-truth pi (gamma=0.35 tempering, same as reported analysis)
Ct      = (C_ref + 1.0) ** GAMMA
pi_true = (Ct / Ct.sum(axis=2, keepdims=True)).reshape(V, K)
y_true  = ilr_forward(pi_true, H)

valid_vox = np.where(maskv & (Nv_real >= 20))[0]

print(f"Reference: {Hdim}x{Wdim}, K={K}, mask={maskv.sum()} voxels")
print(f"Fixed: gamma={GAMMA}, N_EB_ITER={N_EB_ITER}, S={S}")
print()


def _softmax_pi(y):
    z = y @ H.T
    z -= z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


# ============================================================
# Core simulation for one (trim, smooth, lam, seed) condition
# ============================================================
def run_one(trim_alpha, smooth_sigma, lambda_scale, seed):
    rng = np.random.default_rng(seed)

    # Generate S subjects
    counts_all = []
    for _ in range(S):
        y_s  = y_true + SIGMA_SUBJ * rng.standard_normal((V, D))
        pi_s = _softmax_pi(y_s)
        n_s  = np.zeros((V, K), dtype=float)
        for v in valid_vox:
            if Nv_real[v] >= 1:
                n_s[v] = rng.multinomial(Nv_real[v], pi_s[v])
        counts_all.append(n_s)

    # DM baseline
    pi_dm_all = []
    for c in counts_all:
        Ns = c.sum(axis=1, keepdims=True)
        pi_dm_all.append((c + ALPHA_DM) / (np.maximum(Ns, 1) + K * ALPHA_DM))

    # Lambda scaling (same as multiref)
    Ns_arr = np.stack([c.sum(axis=1) for c in counts_all])
    med_N  = float(np.median(Ns_arr[:, maskv]))

    # iCDM-EB
    y_hat_all   = [ilr_forward(p, H) for p in pi_dm_all]
    pi_prop_all = list(pi_dm_all)
    kappa_all   = [np.zeros(V) for _ in range(S)]

    for _iter in range(N_EB_ITER):
        y_stack = np.stack(y_hat_all)
        mu_grp  = trimmed_mean(y_stack, axis=0, alpha=trim_alpha)
        tau_grp = np.maximum(mad_std(y_stack, axis=0), TAU_FLOOR)
        mu_grp  = spatial_smooth_2d(mu_grp, Hdim, Wdim, maskv, sigma=smooth_sigma)
        lam_grp = 1.0 / (tau_grp ** 2)
        lam_sc  = (lambda_scale * med_N * lam_grp
                   / np.maximum(lam_grp.mean(), 1e-6))

        y_hat_all   = []
        kappa_all   = []
        pi_prop_all = []
        for c in counts_all:
            pi_s, y_s, kap_s = icdm_map_voxelwise_vec(
                c, H, mu_grp, lam_sc,
                gamma_init=GAMMA, alpha_blend=ALPHA_BLEND,
                max_iter=30, tol=1e-6)
            y_hat_all.append(y_s)
            kappa_all.append(kap_s)
            pi_prop_all.append(pi_s)

    def vox_mae(pi_list):
        return np.mean(
            [np.abs(p[maskv] - pi_true[maskv]).mean(axis=1) for p in pi_list],
            axis=0)

    mae_dm_vox = vox_mae(pi_dm_all)
    mae_eb_vox = vox_mae(pi_prop_all)

    pct_vs_dm = (mae_dm_vox - mae_eb_vox).mean() / mae_dm_vox.mean() * 100
    diff      = mae_dm_vox - mae_eb_vox
    cohen_d   = diff.mean() / (diff.std() + 1e-12)
    _, p_wil  = wilcoxon(mae_dm_vox, mae_eb_vox, alternative='greater')

    return dict(
        trim_alpha=trim_alpha, smooth_sigma=smooth_sigma,
        lambda_scale=lambda_scale, seed=seed,
        mae_dm=float(mae_dm_vox.mean()),
        mae_eb=float(mae_eb_vox.mean()),
        pct_vs_dm=float(pct_vs_dm),
        cohen_d=float(cohen_d),
        p_wilcoxon=float(p_wil),
    )


# ============================================================
# Sweep 1: TRIM_ALPHA (smooth and lambda fixed at default)
# ============================================================
print("=" * 60)
print("Sweep 1: TRIM_ALPHA  (smooth=%.1f, lambda=%.0f)" %
      (DEF_SMOOTH_SIGMA, DEF_LAMBDA_SCALE))
print("=" * 60)
records_trim = []
for trim in TRIM_VALUES:
    for seed in SEEDS:
        t0 = time.time()
        r  = run_one(trim, DEF_SMOOTH_SIGMA, DEF_LAMBDA_SCALE, seed)
        dt = time.time() - t0
        records_trim.append(r)
        marker = " <- default" if trim == DEF_TRIM_ALPHA else ""
        print(f"  trim={trim:.2f}  seed={seed}  vs_DM={r['pct_vs_dm']:+.1f}%"
              f"  d={r['cohen_d']:.3f}  ({dt:.0f}s){marker}", flush=True)

# ============================================================
# Sweep 2: SMOOTH_SIGMA (trim and lambda fixed at default)
# ============================================================
print()
print("=" * 60)
print("Sweep 2: SMOOTH_SIGMA  (trim=%.2f, lambda=%.0f)" %
      (DEF_TRIM_ALPHA, DEF_LAMBDA_SCALE))
print("=" * 60)
records_smooth = []
for sigma in SMOOTH_VALUES:
    for seed in SEEDS:
        t0 = time.time()
        r  = run_one(DEF_TRIM_ALPHA, sigma, DEF_LAMBDA_SCALE, seed)
        dt = time.time() - t0
        records_smooth.append(r)
        marker = " <- default" if sigma == DEF_SMOOTH_SIGMA else ""
        print(f"  sigma={sigma:.1f}  seed={seed}  vs_DM={r['pct_vs_dm']:+.1f}%"
              f"  d={r['cohen_d']:.3f}  ({dt:.0f}s){marker}", flush=True)

# ============================================================
# Sweep 3: LAMBDA_SCALE (trim and smooth fixed at default)
# ============================================================
print()
print("=" * 60)
print("Sweep 3: LAMBDA_SCALE  (trim=%.2f, smooth=%.1f)" %
      (DEF_TRIM_ALPHA, DEF_SMOOTH_SIGMA))
print("=" * 60)
records_lam = []
for lam in LAMBDA_VALUES:
    for seed in SEEDS:
        t0 = time.time()
        r  = run_one(DEF_TRIM_ALPHA, DEF_SMOOTH_SIGMA, lam, seed)
        dt = time.time() - t0
        records_lam.append(r)
        marker = " <- default" if lam == DEF_LAMBDA_SCALE else ""
        print(f"  lambda={lam:.0f}  seed={seed}  vs_DM={r['pct_vs_dm']:+.1f}%"
              f"  d={r['cohen_d']:.3f}  ({dt:.0f}s){marker}", flush=True)

# ============================================================
# Summarise and save
# ============================================================
all_records = records_trim + records_smooth + records_lam
csv_path = f"{OUTDIR}/hyperparam_results.csv"
fields   = ["trim_alpha", "smooth_sigma", "lambda_scale", "seed",
            "mae_dm", "mae_eb", "pct_vs_dm", "cohen_d", "p_wilcoxon"]
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(all_records)

def summarise(records, key, values, label):
    print(f"\n--- {label} ---")
    print(f"{'Value':>10}  {'MAE-red vs DM':>14}  {'Cohen d':>8}")
    for val in values:
        rows = [r for r in records if r[key] == val]
        m = np.mean([r['pct_vs_dm'] for r in rows])
        s = np.std( [r['pct_vs_dm'] for r in rows])
        d = np.mean([r['cohen_d']   for r in rows])
        marker = " *" if val == {
            "trim_alpha": DEF_TRIM_ALPHA,
            "smooth_sigma": DEF_SMOOTH_SIGMA,
            "lambda_scale": DEF_LAMBDA_SCALE,
        }[key] else ""
        print(f"{val:>10.2f}  {m:>+8.1f}+/-{s:.1f}%  {d:>8.3f}{marker}")

txt_path = f"{OUTDIR}/hyperparam_summary.txt"
import sys
with open(txt_path, "w") as f:
    old = sys.stdout; sys.stdout = f
    summarise(records_trim,   "trim_alpha",   TRIM_VALUES,   "TRIM_ALPHA")
    summarise(records_smooth, "smooth_sigma", SMOOTH_VALUES, "SMOOTH_SIGMA")
    summarise(records_lam,    "lambda_scale", LAMBDA_VALUES, "LAMBDA_SCALE")
    sys.stdout = old

summarise(records_trim,   "trim_alpha",   TRIM_VALUES,   "TRIM_ALPHA")
summarise(records_smooth, "smooth_sigma", SMOOTH_VALUES, "SMOOTH_SIGMA")
summarise(records_lam,    "lambda_scale", LAMBDA_VALUES, "LAMBDA_SCALE")

print(f"\n[SAVED]\n  {csv_path}\n  {txt_path}")
