# iCDM Toolbox

**Individualized Connection Distribution Mapping: A Hierarchical Bayesian Framework for Voxelwise Compositional Connectivity Inference from Diffusion MRI**

## Overview

iCDM is a MATLAB toolbox for estimating individualized voxelwise compositional connectivity distributions from diffusion MRI tractography. The framework models streamline counts as multinomial observations on the simplex and performs inference in isometric log-ratio (ILR) coordinates using Newton-Laplace optimization within an iterative empirical Bayes loop.

### Key features

- **ILR representation**: Maps compositional connectivity from the simplex to unconstrained Euclidean space via Helmert-based ILR transform
- **Newton-Laplace MAP inference**: Voxelwise posterior estimation with multinomial likelihood and Gaussian prior in ILR space
- **Iterative empirical Bayes**: Alternates subject-level inference (E-step) and robust group aggregation (M-step) with trimmed means and MAD-based dispersion
- **Posterior reliability**: Curvature-derived certainty index kappa from the Laplace approximation
- **Covariate modeling**: Reliability-weighted regression of ILR coordinates on subject-level covariates (e.g., age)
- **Spatial regularization**: Optional GMRF prior and neighbourhood smoothing of group-level parameters
- **DARTEL warping**: Native-to-MNI and MNI-to-native spatial normalization via SPM12

## Requirements

- **MATLAB** R2020a or later
- **SPM12** (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) - required for NIfTI I/O, DARTEL warping, and spatial normalization
- **Parallel Computing Toolbox** (optional, for parfor acceleration)
- **FreeSurfer** (optional, for `icdm_fs_mgz_to_vep` parcellation conversion)

## Installation

1. Download or clone this repository
2. Add the `icdmcore/` directory to your MATLAB path:
   ```matlab
   addpath('/path/to/icdmcore');
   ```
3. Ensure SPM12 is on your MATLAB path:
   ```matlab
   addpath('/path/to/spm12');
   ```

## Quick start

### Core pipeline

```matlab
% 1. Initialise
K = 68;  % number of cortical targets (Desikan-Killiany)
opts = struct('gamma', 0.35, 'max_iter', 3, ...);
grp = icdm_init_group_prior(K, outdir, opts);

% 2. Run iterative empirical Bayes
icdm_population_eb(subjects, K, outdir, opts, grp);

% 3. Examine results
icdm_check_group(outdir);
```

### Key functions

| Function | Description |
|----------|-------------|
| `helmert_submatrix` | Orthonormal Helmert basis for ILR transform |
| `ilr_inverse` | Inverse ILR (softmax mapping to simplex) |
| `icdm_subject_vb` | Subject-level Newton-Laplace MAP inference |
| `icdm_population_eb` | Iterative empirical Bayes main loop |
| `icdm_update_group_prior` | Robust group aggregation (M-step) |
| `icdm_estimate_beta` | Covariate regression in ILR space |
| `icdm_gmrf_map` | GMRF-regularised joint MAP inference |
| `spatial_smooth_prior` | Spatial smoothing of group prior |

## File organization

```
icdmcore/
  helmert_submatrix.m      % ILR basis construction
  ilr_inverse.m            % Inverse ILR transform
  icdm_subject_vb.m        % Subject-level inference (E-step)
  icdm_population_eb.m     % Iterative EB main loop
  icdm_update_group_prior.m% Group aggregation (M-step)
  icdm_init_group_prior.m  % Prior initialisation
  icdm_estimate_beta*.m    % Covariate regression
  icdm_gmrf_map.m          % GMRF-regularised MAP
  spatial_smooth_prior.m   % Spatial smoothing
  get_prior_native.m       % Covariate-informed prior
  icdm_warp_to_mni.m       % Native -> MNI warping
  icdm_warp_to_native.m    % MNI -> native warping
  icdm_*_plot_*.m          % Visualization functions
  icdm_kl_js.m             % KL/JS divergence metrics
  study_icdm_train_hbn*.m  % Example study scripts
  LICENSE                  % BSD 3-Clause License
  README.md                % This file
```

## Citation

If you use this software, please cite:

```bibtex
@article{park2026icdm,
  title={Individualized Connection Distribution Mapping: A Hierarchical
         Bayesian Framework for Voxelwise Compositional Connectivity
         Inference from Diffusion MRI},
  author={Park, Hae-Jeong and Kim, Euisun and Park, Jiyoung and
          Eo, Jinseok and Lee, Dongha},
  journal={NeuroImage},
  year={2026}
}
```

## License

This software is released under the BSD 3-Clause License. See [LICENSE](LICENSE) for details.

## Acknowledgments

This software uses the following third-party tools:

- **SPM12** (Wellcome Centre for Human Neuroimaging, UCL) - NIfTI I/O, DARTEL spatial normalization. SPM12 is distributed under the GNU GPL v2. https://www.fil.ion.ucl.ac.uk/spm/
- **FreeSurfer** (Laboratory for Computational Neuroimaging, MGH) - Cortical parcellation. FreeSurfer is freely available for non-commercial use. https://surfer.nmr.mgh.harvard.edu/
- **Healthy Brain Network (HBN)** dataset (Child Mind Institute) - Pediatric diffusion MRI data used for development and validation. https://healthybrainnetwork.org/

## Contact

Hae-Jeong Park, Ph.D.
Department of Nuclear Medicine, Yonsei University College of Medicine
Email: parkhj@yonsei.ac.kr
