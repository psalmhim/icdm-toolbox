# Minimal executable example

Tractogram + atlas → counts → population prior → individualized estimates, on synthetic inputs (about one minute; MATLAB with SPM12 or SPM25 on the path).

```matlab
addpath('/path/to/spm');
addpath('/path/to/icdm-toolbox/examples/minimal_synthetic');
test_counting_rules('/tmp/icdm_rules');          % executable specification of the counting rules
run_minimal_example('/tmp/icdm_minimal');        % full pipeline, checked against expected_outputs.json
```

| File | Role |
|---|---|
| `test_counting_rules.m` | Seven hand-built streamlines on a 12³ atlas; asserts the counts produced by every rule in `docs/CONTRACT.md` (endpoint snapping, along-path subcortical targets, repeated encounters, unassigned streamlines, out-of-grid vertices) |
| `make_synthetic_cohort.m` | 8 subjects: 30 × 30 × 20 atlas (2 mm; six cortical targets, one subcortical structure), 3,000-streamline `.tck` files, identity deformations |
| `run_minimal_example.m` | counts (`icdm_tracks_to_counts`) → validity masks and warps (`icdm_compose_subject`) → prior and MAP/Laplace inference (`icdm_population_eb`, settings of the paper, six model targets) → `N_v`, `kappa_data`, `kappa_post`, `r_data` (`icdm_posterior_support`) |
| `check_outputs.m`, `expected_outputs.json` | reference values: valid-voxel counts and `N_v` totals exactly, median `kappa_data`, `kappa_post`, `r_data` to 1e-4 |
| `write_tck.m`, `make_identity_deformation.m` | helpers |

Outputs are written under the working folder: `inputs/` (synthetic data and `sub-XX/counts.nii`), `icdm/` (prior and per-subject estimates), and `outputs/` (`sub-XX_{kappa_data,kappa_post,r_data,valid}.nii`, `summary.json`).
