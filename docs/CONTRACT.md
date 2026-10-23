# Input/output contract

This document specifies what iCDM reads, how streamlines become counts, and what it writes. Every rule in Section 2 is asserted by `examples/minimal_synthetic/test_counting_rules.m`, and `examples/minimal_synthetic/run_minimal_example.m` checks a full run against `expected_outputs.json`.

## 1. Inputs

| Input | Format | Convention |
|---|---|---|
| Tractogram | MRtrix `.tck` (vertices in world mm). TrackVis `.trk` is read as vertices in voxel units of the label image (not the TrackVis voxmm convention), without per-point scalars or properties; the reported analyses and the examples use `.tck`. | One file per subject (per visit). |
| Label image | Integer NIfTI | Label value `k` is target channel `k`; `0` is unlabelled. It must be in the same world space as the tractogram, and its grid and affine define the count grid. |
| Deformation | SPM deformation `y_*.nii` (`[X Y Z 1 3]`, world mm) | Defined on the template grid. Each template voxel stores the subject-space world coordinate it samples. Expected at `<counts folder>/dartel/y_fs_t1.nii`. |
| Template | NIfTI | Defines the analysis grid (template space). |

## 2. Counting: `core/icdm_tracks_to_counts.m`

1. **Voxel assignment.** Each vertex is mapped to the label grid with `inv(atlas.mat) * [x; y; z; 1]` and rounded to the nearest voxel (single precision). Vertices outside the grid are dropped.
2. **Endpoint assignment.** If an end vertex lies in an unlabelled voxel, the streamline is extended along its terminal direction by one and then two voxels. Labelled voxels reached in this way are appended; the unlabelled gap voxel is not.
3. **Along-path assignment.** The streamline's target set is the set of distinct labels over all of its voxels. It therefore includes its endpoint regions and any labelled structure it passes through (for example a subcortical nucleus).
4. **Repeated encounters.** Each target is counted once per streamline, however many of its voxels the streamline visits. Consecutive repeats of the same voxel collapse to one visit. A voxel re-entered after the streamline has left it is visited again and counted again.
5. **Unassigned streamlines.** A streamline whose target set is empty contributes nothing. A streamline with one assigned end contributes to that target only.
6. **Weighting.** Every visited voxel receives `+1` in each channel of the streamline's target set. Streamline weights (SIFT/SIFT2) are not applied; the reported analyses used unfiltered tractograms.
7. **Storage.** Counts are `uint16` (maximum 65,535 per voxel and channel). The highest single-channel count in the reported HBN data was below 35,000.

The resulting `n_{v,k}` counts streamline–target incidences. `N_v = sum_k n_{v,k}` exceeds the number of distinct streamlines whenever streamlines have more than one target.

## 3. Target channels

`atlas/icdm_remap_dk84.m` converts a FreeSurfer `aparc+aseg` image to the 84-target label image. Channel order (`atlas/target_lookup_84.tsv`):

| Channels | Content | FreeSurfer labels |
|---|---|---|
| 1–34 | left Desikan–Killiany cortex | 1001–1003, 1005–1035 (1004 corpus callosum excluded) |
| 35–68 | right Desikan–Killiany cortex | 2001–2003, 2005–2035 |
| 69–84 | subcortical and cerebellar | 10, 11, 12, 13, 17, 18, 26, 49, 50, 51, 52, 53, 54, 58, 8, 47 |

This is **not** the MRtrix `fs_default` order. The 68-target analyses (HBN, HCP) model channels 1–68 (`atlas/target_lookup_68.tsv`). The tumour analyses model all 84 channels. Channels outside the modelled set still enter the native validity threshold (Section 5).

## 4. Spatial conventions

- The count image inherits the label image's grid and affine. The reported analyses used FreeSurfer labels resampled (nearest neighbour) to a 1.5 mm grid aligned with the diffusion data.
- Subject-to-template mapping uses the SPM deformation (Section 1). `icdm_compose_subject` builds trilinear native→template and template→native operators from it. For a chain of transforms (for example subject → age-specific template → common space), compose the fields first; see `docs/NEW_POPULATION.md`.
- Composition (`H`, Helmert basis `helmert_submatrix(K)`, `K x (K-1)`): ILR coordinates are `y = H' log(pi)` and `pi = softmax(H y)`.

## 5. Validity masks

| Mask | Definition | Where |
|---|---|---|
| Native mask | `sum over all channels in the count image > 5` (`opts.thresh`) | `icdm_subj_<id>.mat`: `idx_native` |
| Group analysis domain | template voxels covered by the native masks of at least 30% of subjects | `group_mask.mat`, `grp.idx_mni` |
| Subject template mask | warp support of the native mask `>= 0.5` | `<id>_vb.mat`: `valid_mni` |
| Common-valid set | `valid_mni` true in every subject | used for voxelwise group statistics |

## 6. Outputs

Per subject, `<outdir>/iter_002/<id>_vb.mat` (template space, `grp.idx_mni` order). Fields are NaN outside `valid_mni`.

| Field | Size | Meaning |
|---|---|---|
| `y_ilr_mni` | `N x (K-1)` | individualized estimate: MAP ILR coordinates |
| `y_data_mni` | `N x (K-1)` | count-derived (naive) composition `ilr((n + 1/2) / (N + K/2))`, prior-independent |
| `kappa_data_mni` | `N x 1` | `kappa_data`: prior-independent data information, `tr(N_v H' S(pi_bar) H) / (K-1)` at the relative frequencies `pi_bar = n / N_v` |
| `kappa_mni` | `N x 1` | `kappa_post`: posterior certainty, mean posterior curvature at the MAP |
| `valid_mni` | `N x 1` | subject template mask |

`kappa_data_mni` and `kappa_mni` are set to 1 (`kappa_base`) at valid voxels where the warped value is not positive.

Population prior, `<outdir>/group_icdm_iter_002.mat`: `mu_ilr_mni` (prior mean), `tau2` (between-subject dispersion), `kappa_mni` (`lambda_grp`; the prior precision applied is `alpha_grp * lambda_grp * I`), `idx_mni`, `dim_mni`, `H`.

`core/icdm_posterior_support.m` returns, per subject:

| Output | Meaning |
|---|---|
| `N_mni` | `N_v` warped to template space with the subject's operator |
| `kappa_data`, `kappa_post` | as above |
| `r_data` | posterior data fraction `tr(Q_data) / tr(Q_post)` at the MAP: `N c(pi_hat) / (N c(pi_hat) + alpha_grp * lambda_grp)`, with `c(pi) = mean diag(H' (diag(pi) - pi pi') H)`; `0` where `N_v = 0` |
| `valid` | `valid_mni` |

`examples/minimal_synthetic/run_minimal_example.m` writes these as `outputs/sub-XX_{kappa_data,kappa_post,r_data,valid}.nii` and summarizes them in `outputs/summary.json`. The reference values are in `expected_outputs.json`.

## 7. Settings of the reported analyses

| Setting | Value |
|---|---|
| Population-transfer strength `alpha_grp` | 0.25 (HBN; calibrate for new data with `icdm_select_alpha`) |
| Trimming fraction | 0.20 |
| Group-prior smoothing `eta` | 0.5, 6-connected |
| Variance floor / precision cap | 1e-6 / 1e4 |
| Newton–CG | at most 20 iterations, 20 CG iterations, gradient and step tolerance 1e-6 |
| Native threshold / group coverage | > 5 incidences / 30% of subjects |
