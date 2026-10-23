# Applying iCDM to a new population or template

This guide follows `examples/infant_unc4d/run_unc4d_example.m`, a configurable longitudinal workflow for infant data with age-specific UNC-4D templates (`unc4d_config.json`). The example runs on small synthetic inputs. It demonstrates configuration and data organisation; it is not a validation of iCDM in infant data.

UNC-4D refers to the UNC/UMN 4D infant brain volumetric atlas built from the Baby Connectome Project cohort:

> Chen L, Wu Z, Hu D, Wang Y, Zhao F, Zhong T, Lin W, Wang L, Li G. A 4D infant brain volumetric atlas based on the UNC/UMN Baby Connectome Project (BCP) cohort. *NeuroImage*. 2022;253:119097. https://doi.org/10.1016/j.neuroimage.2022.119097

The atlas is distributed at https://www.nitrc.org/projects/uncbcp_4d_atlas/. See that page for the current release (later versions extend the age range beyond the published one) and its terms of use.

## 1. Upstream data

- Use one tractography protocol for every subject and visit (algorithm, streamline budget, seeding, termination). iCDM models the counts this protocol produces; it does not correct differences between protocols.
- Provide a target label image per visit in the tractogram's world space, with the channel order of `atlas/target_lookup_68.tsv` or `_84.tsv` (Section 3 of `docs/CONTRACT.md`). Infant parcellations must be mapped to the same channel list at every age. The UNC/UMN 4D atlas provides the Desikan cortical parcellation, mapped consistently across ages, together with subcortical labels, which can serve as the basis for the Desikan–Killiany target channels.

## 2. Template selection

Assign each visit the age-appropriate template, i.e. the one whose age bin contains the visit (`templates[].age_months` in the configuration). With UNC-4D, list the month-specific templates used by your registration pipeline and their age bins. A visit outside every bin stops the run instead of being matched silently.

## 3. Mapping to a common analysis space

Group analysis needs a single space. Choose one common template (e.g. one UNC-4D month) and provide, for every age-specific template, the deformation from the common space to it. For each visit the example composes

- `y_to_<template>.nii`: subject → age-specific template (defined on the age-specific template grid), and
- `y_common_to_<template>.nii`: common space → age-specific template (defined on the common grid)

into a single field on the common grid with SPM (`spm_deformations`, `comp = {subject field, common field}`; SPM evaluates the later field first). The result is saved where `icdm_compose_subject` expects it (`<counts folder>/dartel/y_fs_t1.nii`). Deformations from other software must first be converted to SPM `y_*.nii` format.

## 4. Reference sample and prior construction

- Mark reference visits in `sessions.tsv` (`reference = 1`). Typically these are participants who meet the study's inclusion and quality-control criteria and represent the population the prior should describe. The prior is built once from their count-derived representations and then held fixed. Visits outside the reference sample are estimated under it without contributing to it.
- **Developmental trajectories (association testing):** use one covariate-agnostic prior pooled across ages (`prior_mode: "pooled"`), balanced across the age range. Age must not enter the prior when age is the effect under test.
- **Individualized reference (prediction):** an age-specific prior (one per age bin) or a covariate-conditioned prior mean can be used. It encodes age and is therefore unsuitable for testing age effects.
- The reference sample should be large enough for stable trimmed means and dispersions. The paper's cohort-size simulation characterises this dependence.

## 5. Parameter calibration

| Parameter | Recommendation |
|---|---|
| `alpha_grp` | Calibrate on the reference sample with `core/icdm_select_alpha.m` (cross-validated held-out predictive likelihood; the in-sample marginal likelihood is also reported). Neither criterion uses covariates. For repeated visits, fold by participant (`opts.groups`). HBN: 0.26 (cross-validated) and 0.38 (marginal), with 0.25 adopted. |
| Trimming fraction, smoothing `eta`, variance floor, precision cap | Keep the defaults (0.20, 0.5, 1e-6, 1e4); they are fixed constants of the method, not tuned per dataset. |
| Native threshold (> 5 incidences), group coverage (30%) | Data-domain thresholds. Sparser tractograms (smaller brains, lower streamline budgets) may need different values; fix them before any analysis of the effect of interest. |
| ILR basis, channel order | Determined by the target list; not tuned. |

## 6. Inference and outputs

Every visit is estimated under the fixed prior with the calibrated `alpha_grp`. Outputs are organised by participant and visit:

```
derivatives/icdm/
  calibration.json                  alpha_grp, criterion, reference size, fixed settings
  prior/                            population prior (group_icdm_iter_002.mat) and its reference fits
  final/                            <sub>_<ses>_vb.mat for every visit under the fixed prior and calibrated alpha_grp
  sub-XX/ses-YY/
    counts.nii                      streamline-target incidence counts
    dartel/y_fs_t1.nii              composed deformation to the common space
    icdm_vb.mat                     y_ilr_mni, y_data_mni, kappa_data_mni, kappa_mni, valid_mni
    kappa_data.nii, kappa_post.nii, r_data.nii, valid.nii
  trajectory_table.tsv              one row per visit: participant, visit, age, template, reference flag,
                                    alpha_grp, valid voxels, median kappa_data / kappa_post / r_data, estimate file
```

`trajectory_table.tsv` indexes the per-visit estimates for longitudinal models, for example mixed-effects models of the ILR coordinates with a participant random effect.
