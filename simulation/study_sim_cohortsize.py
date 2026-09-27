#!/usr/bin/env python3
"""
Cohort-size sensitivity analysis for iCDM-EB.
Question: does the empirical-Bayes advantage depend on cohort size S? (Supplementary Fig. S2(a)-(c))

For each S in [5, 10, 20, 50, 100] and each of N_SEEDS random seeds:
  - Generate S virtual subjects from icdm84.mat ground-truth
  - Run iCDM-EB (empirical Bayes) and baseline estimators
  - Record MAE reductions and Cohen's d vs DM

Outputs:
  figures/cohortsize_results.csv    per-condition results
  figures/cohortsize_summary.txt    human-readable summary
  figures/cohortsize_figure.png     publication figure
"""

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import rel_entr
from scipy.stats import wilcoxon
from scipy.ndimage import gaussian_filter
import os, time

# ============================================================
# CONFIG
# ============================================================
MAT_PATH   = "icdm84.mat"
OUTDIR     = "figures/R1.2_cohortsize"
os.makedirs(OUTDIR, exist_ok=True)

S_VALUES   = [5, 10, 20, 50, 100]
N_SEEDS    = 3                # repeat each S with 3 different random seeds
SEEDS      = [42, 43, 44]

ALPHA_DM     = 1.0
N_EB_ITER    = 1
SMOOTH_SIGMA = 1.0
LAMBDA_SCALE = 25.0
ALPHA_BLEND  = 0.5
TRIM_ALPHA   = 0.2
TAU_FLOOR    = 0.05
SIGMA_SUBJ   = 0.15
GAMMA        = 0.35

# ============================================================
# Basis helpers  (identical to study_sim_multiref.py)
# ============================================================
def helmert_basis(K):
    # Orthonormal Helmert ILR basis: H^T H = I, H^T 1 = 0.
    H = np.zeros((K, K - 1))
    for j in range(K - 1):
        i = j + 1
        H[:i, j] = 1.0 / np.sqrt(i * (i + 1))
        H[i, j]  = -np.sqrt(i / (i + 1))
    assert np.abs(H.T @ H - np.eye(K - 1)).max() < 1e-10, "Helmert basis not orthonormal"
    return H

def _softmax_pi(y, H):
    z = y @ H.T
    z -= z.max(axis=1, keepdims=True)
    e  = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)

def ilr_forward(pi, H):
    pi  = np.clip(pi, 1e-12, 1.0)
    clr = np.log(pi) - np.log(pi).mean(axis=1, keepdims=True)
    return clr @ H

def trimmed_mean(arr, axis=0, alpha=0.2):
    n  = arr.shape[axis]
    lo = int(np.floor(alpha * n))
    hi = n - lo
    return np.mean(np.sort(arr, axis=axis).take(range(lo, hi), axis=axis), axis=axis)

def mad_std(arr, axis=0):
    med = np.median(arr, axis=axis, keepdims=True)
    return 1.4826 * np.median(np.abs(arr - med), axis=axis)

def spatial_smooth_2d(field, Hd, Wd, maskv, sigma=1.0):
    out = field.copy()
    weight = maskv.reshape(Hd, Wd).astype(float)
    for ch in range(field.shape[1]):
        img = field[:, ch].reshape(Hd, Wd)
        num = gaussian_filter(np.where(weight > 0, img, 0.0), sigma=sigma)
        den = gaussian_filter(weight, sigma=sigma)
        out[:, ch] = (num / np.maximum(den, 1e-12)).ravel()
    return out

def icdm_map_voxelwise_vec(counts, H, prior_mean, prior_prec,
                            gamma_init=0.35, alpha_blend=0.5,
                            max_iter=30, tol=1e-6):
    """Vectorized Newton-Laplace MAP for all voxels simultaneously."""
    V, K  = counts.shape
    D     = K - 1
    Nv    = counts.sum(axis=1)
    lam   = prior_prec if prior_prec.ndim == 2 else np.tile(prior_prec[:, None], (1, D))

    ct     = (counts + 1.0) ** gamma_init
    pi_init = ct / ct.sum(axis=1, keepdims=True)
    y_init  = ilr_forward(pi_init, H)
    y       = (1 - alpha_blend) * prior_mean + alpha_blend * y_init

    active = Nv >= 20
    pi_hat = pi_init.copy()
    y_hat  = y_init.copy()
    kappa  = np.zeros(V)

    for _ in range(max_iter):
        idx  = np.where(active)[0]
        if idx.size == 0:
            break
        n_a   = counts[idx]
        N_a   = Nv[idx]
        y_a   = y[idx]
        m_a   = prior_mean[idx]
        lam_a = lam[idx]

        piv  = _softmax_pi(y_a, H)
        grad = (n_a - N_a[:, None] * piv) @ H - lam_a * (y_a - m_a)

        sqH  = np.sqrt(piv)[:, :, None] * H[None, :, :]
        HpH  = np.matmul(sqH.transpose(0, 2, 1), sqH)
        Hp   = piv @ H
        HpH -= Hp[:, :, None] * Hp[:, None, :]

        Q    = lam_a[:, :, None] * np.eye(D)[None] + N_a[:, None, None] * HpH
        dy   = np.linalg.solve(Q, grad[:, :, None]).squeeze(-1)
        y[idx] += dy
        active[idx[np.linalg.norm(dy, axis=1) < tol]] = False

    valid = Nv >= 20
    if valid.any():
        vi   = np.where(valid)[0]
        pf   = _softmax_pi(y[vi], H)
        sqHf = np.sqrt(pf)[:, :, None] * H[None, :, :]
        HpHf = np.matmul(sqHf.transpose(0, 2, 1), sqHf)
        Hpf  = pf @ H
        HpHf -= Hpf[:, :, None] * Hpf[:, None, :]
        Qf   = lam[vi, :, None] * np.eye(D)[None] + Nv[vi, None, None] * HpHf
        kappa[vi]   = np.trace(Qf, axis1=1, axis2=2) / D
        pi_hat[vi]  = pf
        y_hat[vi]   = y[vi]

    return pi_hat, y_hat, kappa

# ============================================================
# Load reference data
# ============================================================
mat   = sio.loadmat(MAT_PATH)
C_ref = mat["C"][:, :, :68].astype(float)
Hdim, Wdim, K = C_ref.shape
V    = Hdim * Wdim
D    = K - 1
H    = helmert_basis(K)

maskv    = (C_ref.sum(axis=2) >= 20).ravel()
Nv_real  = C_ref.reshape(V, K).sum(axis=1).astype(int)

# Ground-truth π from gamma-tempered counts
Ct      = (C_ref + 1.0) ** GAMMA
pi_true = (Ct / Ct.sum(axis=2, keepdims=True)).reshape(V, K)
y_true  = ilr_forward(pi_true, H)

print(f"Reference: V={V}, mask={maskv.sum()} voxels, K={K}")
print(f"Running S in {S_VALUES}, {N_SEEDS} seeds each")
print()

# ============================================================
# Core simulation for one (S, seed)
# ============================================================
def run_one(S, seed):
    rng = np.random.default_rng(seed)

    # --- Generate S subjects ---
    counts_all = []
    for _ in range(S):
        y_s   = y_true + SIGMA_SUBJ * rng.standard_normal((V, D))
        z_s   = y_s @ H.T
        z_s  -= z_s.max(axis=1, keepdims=True)
        e_s   = np.exp(z_s)
        pi_s  = e_s / e_s.sum(axis=1, keepdims=True)
        n_s   = np.zeros((V, K), dtype=float)
        for v in np.where(maskv & (Nv_real >= 20))[0]:
            n_s[v] = rng.multinomial(Nv_real[v], pi_s[v])
        counts_all.append(n_s)

    # --- Naive & DM baselines ---
    pi_naive_all, pi_dm_all = [], []
    for c in counts_all:
        Ns = c.sum(axis=1, keepdims=True)
        pi_naive_all.append(c / np.maximum(Ns, 1))
        pi_dm_all.append((c + ALPHA_DM) / (np.maximum(Ns, 1) + K * ALPHA_DM))

    # --- iCDM-EB ---
    y_hat_all   = [ilr_forward(p, H) for p in pi_dm_all]
    kappa_all   = [np.ones(V)] * S
    pi_prop_all = None

    for _ in range(N_EB_ITER):
        y_stack = np.stack(y_hat_all)
        mu_grp  = trimmed_mean(y_stack, axis=0, alpha=TRIM_ALPHA)
        tau_grp = np.maximum(mad_std(y_stack, axis=0), TAU_FLOOR)
        mu_grp  = spatial_smooth_2d(mu_grp, Hdim, Wdim, maskv, sigma=SMOOTH_SIGMA)
        lam_grp = 1.0 / (tau_grp ** 2)

        Ns_arr = np.stack([c.sum(axis=1) for c in counts_all])
        med_N  = np.median(Ns_arr[:, maskv])
        lam_sc = LAMBDA_SCALE * med_N * lam_grp / np.maximum(lam_grp.mean(), 1e-6)

        y_hat_all   = []
        kappa_all   = []
        pi_prop_all = []
        for c in counts_all:
            pi_s, y_s, kap_s = icdm_map_voxelwise_vec(
                c, H, mu_grp, lam_sc,
                gamma_init=GAMMA, alpha_blend=ALPHA_BLEND)
            y_hat_all.append(y_s)
            kappa_all.append(kap_s)
            pi_prop_all.append(pi_s)

    # Per-voxel MAE averaged over S subjects  (shape: n_mask_voxels,)
    # This gives stable Cohen's d regardless of S (n~3600, not S)
    def vox_mae(pi_list):
        return np.mean([np.abs(p[maskv] - pi_true[maskv]).mean(axis=1)
                        for p in pi_list], axis=0)

    mae_naive_vox = vox_mae(pi_naive_all)
    mae_dm_vox    = vox_mae(pi_dm_all)
    mae_eb_vox    = vox_mae(pi_prop_all)

    pct_vs_naive = (mae_naive_vox - mae_eb_vox).mean() / mae_naive_vox.mean() * 100
    pct_vs_dm    = (mae_dm_vox    - mae_eb_vox).mean() / mae_dm_vox.mean()    * 100

    diff    = mae_dm_vox - mae_eb_vox
    cohen_d = diff.mean() / (diff.std() + 1e-12)

    _, p_wil = wilcoxon(mae_dm_vox, mae_eb_vox, alternative='greater')
    win_rate  = np.mean(mae_eb_vox < mae_dm_vox) * 100

    return dict(S=S, seed=seed,
                mae_naive=mae_naive_vox.mean(), mae_dm=mae_dm_vox.mean(),
                mae_eb=mae_eb_vox.mean(),
                pct_vs_naive=pct_vs_naive, pct_vs_dm=pct_vs_dm,
                cohen_d=cohen_d, win_rate=win_rate, p_wilcoxon=p_wil)

# ============================================================
# Main loop
# ============================================================
records = []
t_start = time.time()

for S in S_VALUES:
    for seed in SEEDS:
        t0 = time.time()
        r  = run_one(S, seed)
        dt = time.time() - t0
        records.append(r)
        print(f"  S={S:3d}  seed={seed}: vs_DM={r['pct_vs_dm']:.1f}%  "
              f"d={r['cohen_d']:.3f}  win={r['win_rate']:.1f}%  ({dt:.0f}s)", flush=True)

print(f"\nTotal: {time.time()-t_start:.0f}s")

# ============================================================
# Aggregate by S
# ============================================================
import csv

csv_path = f"{OUTDIR}/cohortsize_results.csv"
fields   = list(records[0].keys())
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(records)

S_arr            = np.array(S_VALUES)
pct_dm_mean      = np.zeros(len(S_VALUES))
pct_dm_sd        = np.zeros(len(S_VALUES))
pct_naive_mean   = np.zeros(len(S_VALUES))
pct_naive_sd     = np.zeros(len(S_VALUES))
cohend_mean      = np.zeros(len(S_VALUES))
cohend_sd        = np.zeros(len(S_VALUES))
win_mean         = np.zeros(len(S_VALUES))
win_sd           = np.zeros(len(S_VALUES))

for i, S in enumerate(S_VALUES):
    rows = [r for r in records if r["S"] == S]
    pct_dm_mean[i]    = np.mean([r["pct_vs_dm"]    for r in rows])
    pct_dm_sd[i]      = np.std( [r["pct_vs_dm"]    for r in rows])
    pct_naive_mean[i] = np.mean([r["pct_vs_naive"]  for r in rows])
    pct_naive_sd[i]   = np.std( [r["pct_vs_naive"]  for r in rows])
    cohend_mean[i]    = np.mean([r["cohen_d"]       for r in rows])
    cohend_sd[i]      = np.std( [r["cohen_d"]       for r in rows])
    win_mean[i]       = np.mean([r["win_rate"]      for r in rows])
    win_sd[i]         = np.std( [r["win_rate"]      for r in rows])

# ============================================================
# Summary text
# ============================================================
summary_path = f"{OUTDIR}/cohortsize_summary.txt"
sep = "=" * 60
with open(summary_path, "w") as f:
    f.write(sep + "\n")
    f.write("COHORT-SIZE SENSITIVITY (iCDM-EB vs DM)\n")
    f.write(f"{N_SEEDS} seeds per S level\n")
    f.write(sep + "\n\n")
    f.write(f"{'S':>6}  {'vs DM (%)':>14}  {'vs Naive (%)':>14}  {'Cohen d':>10}  {'Win (%)':>10}\n")
    f.write("-" * 60 + "\n")
    for i, S in enumerate(S_VALUES):
        f.write(f"{S:6d}  "
                f"{pct_dm_mean[i]:6.1f}±{pct_dm_sd[i]:.1f}  "
                f"{pct_naive_mean[i]:8.1f}±{pct_naive_sd[i]:.1f}  "
                f"{cohend_mean[i]:8.3f}±{cohend_sd[i]:.3f}  "
                f"{win_mean[i]:7.1f}±{win_sd[i]:.1f}\n")
    f.write("\n" + sep + "\n")

with open(summary_path) as f:
    print(f.read())

# ============================================================
# Figure: 3-panel publication figure
# ============================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), dpi=150)
fig.subplots_adjust(wspace=0.38)

x_log = np.log10(S_arr)
xt    = S_arr
xl    = [str(s) for s in S_arr]

# Mark S=20 reference line
ref_idx = S_VALUES.index(20)

# Panel A: MAE reduction vs DM
ax = axes[0]
ax.errorbar(S_arr, pct_dm_mean, yerr=pct_dm_sd,
            fmt="-o", color="#1f77b4", linewidth=1.8, markersize=6,
            capsize=4, capthick=1.2, label="vs DM")
ax.axvline(20, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(pct_dm_mean[ref_idx], color="gray", linestyle=":", linewidth=0.7, alpha=0.5)
ax.set_xscale("log")
ax.set_xticks(S_arr); ax.set_xticklabels(xl)
ax.set_xlabel("Cohort size S", fontsize=10)
ax.set_ylabel("MAE reduction vs DM (%)", fontsize=10)
ax.set_title("A  EB gain over DM", fontweight="bold", fontsize=10, loc="left")
ax.set_ylim(bottom=0)
ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)

# Panel B: MAE reduction vs Naive
ax = axes[1]
ax.errorbar(S_arr, pct_naive_mean, yerr=pct_naive_sd,
            fmt="-s", color="#ff7f0e", linewidth=1.8, markersize=6,
            capsize=4, capthick=1.2, label="vs Naive")
ax.axvline(20, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(pct_naive_mean[ref_idx], color="gray", linestyle=":", linewidth=0.7, alpha=0.5)
ax.set_xscale("log")
ax.set_xticks(S_arr); ax.set_xticklabels(xl)
ax.set_xlabel("Cohort size S", fontsize=10)
ax.set_ylabel("MAE reduction vs Naive (%)", fontsize=10)
ax.set_title("B  EB gain over Naive", fontweight="bold", fontsize=10, loc="left")
ax.set_ylim(bottom=0)
ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)

# Panel C: Cohen's d
ax = axes[2]
ax.errorbar(S_arr, cohend_mean, yerr=cohend_sd,
            fmt="-^", color="#2ca02c", linewidth=1.8, markersize=6,
            capsize=4, capthick=1.2, label="Cohen's d")
ax.axvline(20, color="gray", linestyle="--", linewidth=0.8, alpha=0.7,
           label="S=20 (manuscript)")
ax.axhline(cohend_mean[ref_idx], color="gray", linestyle=":", linewidth=0.7, alpha=0.5)
ax.set_xscale("log")
ax.set_xticks(S_arr); ax.set_xticklabels(xl)
ax.set_xlabel("Cohort size S", fontsize=10)
ax.set_ylabel("Cohen's d (EB vs DM)", fontsize=10)
ax.set_title("C  Effect size", fontweight="bold", fontsize=10, loc="left")
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, frameon=False)
ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)

# Annotate S=20 value on each panel
for ax_i, yvals, fmt in zip(axes,
                              [pct_dm_mean, pct_naive_mean, cohend_mean],
                              ["{:.0f}%", "{:.0f}%", "{:.2f}"]):
    y20 = yvals[ref_idx]
    ax_i.annotate(fmt.format(y20),
                  xy=(20, y20), xytext=(6, 4), textcoords="offset points",
                  fontsize=8, color="gray")

fig.suptitle(
    f"iCDM-EB cohort-size sensitivity  (N_seeds={N_SEEDS}, ref=NDARAT100AEQ)",
    fontsize=9, y=1.02, color="dimgray")

fig_path = f"{OUTDIR}/cohortsize_figure.png"
fig.savefig(fig_path, bbox_inches="tight", dpi=150)
fig.savefig(fig_path.replace(".png", ".pdf"), bbox_inches="tight")
plt.close()

print(f"\n[SAVED]")
print(f"  {csv_path}")
print(f"  {summary_path}")
print(f"  {fig_path}")
