# Longitudinal infant example (UNC-4D templates)

A configurable workflow for applying iCDM to infant data with age-specific templates, run here on a small synthetic dataset (5 participants × 3 visits, about 3, 6 and 12 months). It demonstrates configuration and data organisation; it is not a validation of iCDM in infant data.

```matlab
addpath('/path/to/spm');
addpath('/path/to/icdm-toolbox/examples/infant_unc4d');
TT = run_unc4d_example('/tmp/icdm_unc4d');                  % synthetic data + unc4d_config.json
TT = run_unc4d_example('/my/study', 'my_config.json');      % your data, same layout
```

`unc4d_config.json` sets the template age bins and files, the deformations, the common analysis space, the reference column, the prior mode, and the calibration. Replace the template paths with your local UNC-4D files. The steps, from template selection to `trajectory_table.tsv`, are described in `docs/NEW_POPULATION.md`.
