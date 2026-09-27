# iCDM Toolbox

MATLAB implementation of **individualized connection distribution mapping (iCDM)**, the method described in

> Park H-J, Kim E, Park J, Lee J, Lee D, Eo J, Lee D, Jeong S-O. *Voxelwise White-Matter Compositional Connectivity Mapping from Diffusion MRI via Reliability-Aware Hierarchical Inference.* Medical Image Analysis (under review).

iCDM treats the voxelwise streamline--target counts produced by a tractography pipeline as finite evidence about a latent, subject-specific connectivity composition. Inference is performed in isometric log-ratio (ILR) coordinates by combining each subject's counts with a spatially varying empirical-Bayes population prior, and every estimate is reported together with the amount and source of its support.

## Method in brief

- **Likelihood.** Multinomial working likelihood for the relative streamline--target incidences at each voxel, written in ILR coordinates (Helmert basis).
- **Population prior.** Estimated once from prior-independent, count-derived subject representations warped to template space: trimmed mean, MAD-based between-subject dispersion, and isotropic precision
  `Lambda_grp = alpha_grp * lambda_grp * I`, with `lambda_grp = 1 / max(mean(tau^2), eps)` (capped at 1e4) and spatial smoothing of the group parameters. The prior is held fixed during subject-level inference; no posterior quantity is fed back into it.
- **Subject-level inference.** Voxelwise MAP estimation by Newton--conjugate-gradient optimization and a Laplace approximation of the posterior, in native diffusion space.
- **Reported indices.**
  - `kappa_data` -- prior-independent data information (likelihood curvature at the relative streamline frequencies);
  - `kappa_post` -- posterior certainty (mean posterior curvature at the MAP);
  - `r_data = tr(Q_data) / tr(Q_post)` -- fraction of posterior precision supplied by the subject's own counts. With the isotropic prior, `r_data = 1 - alpha_grp * lambda_grp / kappa_post`.
- **Subject characteristics.** Two modes: (i) *association testing* -- the characteristic under test is excluded from the prior and estimated downstream on the inferred ILR coordinates (ordinary least squares, permutation-based TFCE / cluster-extent FWE, Freedman--Lane for nuisance covariates); (ii) *individualized reference* -- an independently learned or cross-fitted coefficient field `beta_pred` shifts the prior mean (optional; off by default).

## Requirements

- MATLAB R2020a or later (Parallel Computing Toolbox optional, for `parfor`)
- SPM (SPM12 or SPM25) for NIfTI I/O and deformation fields
- Tractography and parcellation are upstream of iCDM (e.g. MRtrix3, FreeSurfer) and are not part of this toolbox
- Simulations: Python 3 with `numpy`, `scipy`, `matplotlib`

## Repository layout

```
core/         method (MATLAB)
examples/     run_icdm_example.m -- end-to-end driver with the settings used in the paper
simulation/   simulation scripts (Python)
```

### `core/`

| file | role |
|---|---|
| `icdm_compose_subject.m` | voxelize a count map into streamline--target incidence counts and prepare the subject record |
| `helmert_submatrix.m` | orthonormal Helmert basis of the ILR transform |
| `icdm_evaluate_group_mask.m`, `icdm_build_subject_warps.m` | template-space analysis domain and native <-> template warp tables |
| `icdm_build_warp_to_mni.m`, `icdm_build_warp_to_native.m`, `icdm_warp_to_mni.m`, `icdm_warp_to_native.m`, `icdm_warp_4d.m`, `icdm_apply_warp.m`, `icdm_trilinear_index_weight.m` | spatial transport of ILR fields and prior parameters |
| `icdm_init_group_prior.m` | initialize the population-prior structure |
| `icdm_population_eb.m` | driver: population prior from count-derived representations, then subject-level inference under the fixed prior |
| `icdm_update_group_prior.m` | robust group mean and dispersion, group precision, spatial smoothing |
| `icdm_run_subject_vb.m`, `icdm_subject_vb.m` | subject-level Newton--CG MAP, Laplace curvature, `kappa_data`, `kappa_post` |
| `icdm_prepare_design_vector.m`, `icdm_design_vector.m` | design assembly for downstream association |
| `icdm_estimate_beta.m`, `icdm_beta_estimator.m`, `icdm_fl_oneperm.m`, `icdm_tfce.m` | downstream association on inferred ILR coordinates with permutation inference |
| `getfield_default.m` | option defaults |

## Usage

```matlab
addpath(genpath('/path/to/icdm-toolbox'));
addpath('/path/to/spm');
covariates.names  = {'age','sex'};
covariates.values = [age(:) sex(:)];
run_icdm_example(icdm_files, names, covariates, '/path/to/template.nii', 'out');
```

`icdm_files` are per-subject 4-D images of streamline--target counts (voxel x target) in native diffusion space. The example uses the settings of the reported analyses: `alpha_grp = 0.25`, trimming fraction 0.2, group-prior smoothing `eta = 0.5` (6-connected), variance floor 1e-6, precision cap 1e4, Newton--CG with at most 20 iterations. The population-transfer strength `alpha_grp` should be calibrated for each new dataset with criteria that do not use the effect under test (e.g. cross-validated predictive likelihood).

Per-subject outputs (`<outdir>/iter_002/<name>_vb.mat`) include the MAP ILR coordinates (`y_ilr_mni`), the count-derived coordinates (`y_data_mni`), `kappa_data_mni`, `kappa_mni` (posterior certainty) and the validity mask `valid_mni`.

## Simulations (`simulation/`)

| script | purpose |
|---|---|
| `study_sim_multiref_full.py` (with helper `study_sim_multiref.py`) | compositional recovery across ten reference anatomies (manuscript Table 2) |
| `study_sim_cohortsize.py` | sensitivity to cohort size |
| `study_sim_hyperparam_sensitivity.py` | one-at-a-time sensitivity to trimming, group-prior smoothing and prior scale |
| `predictive_prior_sim.py` | covariate-conditioned (individualized) reference; fully synthetic |

The reference-based scripts construct the simulation ground truth from voxelwise streamline--target counts of Healthy Brain Network participants. These derived data are not redistributed here: `icdm84.mat` holds `C`, an `H x W x 84` count array from one coronal slice, and the multi-reference script reads per-participant count maps from the folder given by the environment variable `ICDM_HBN_BASE`. HBN data are available from the Healthy Brain Network under its data-use terms.

## License

BSD 3-Clause (see `LICENSE`).
