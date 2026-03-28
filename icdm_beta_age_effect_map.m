% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the copyright holder nor the names of its
%      contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES.
%
% REQUIRES: SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

function E = icdm_beta_age_effect_map(Beta, age_index)
% ICDM_BETA_AGE_EFFECT_MAP  Compute voxelwise age-effect magnitude from beta fields.
%
%   E = icdm_beta_age_effect_map(Beta, age_index)
%
%   Computes the L2-norm of the age-related beta coefficients across all
%   ILR dimensions at every MNI voxel, yielding a scalar summary of how
%   strongly connectivity changes with age at each location.
%
%   When age_index is a scalar, the corresponding row of Beta is extracted
%   and the effect magnitude is sqrt(sum(beta_age.^2, 2)).  When
%   age_index is a vector (e.g., multiple spline basis columns), the L2
%   norm is taken jointly across all basis functions and ILR dimensions.
%
%   Inputs
%     Beta      : [P x Nmni x (K-1)] voxelwise regression coefficients.
%     age_index : scalar or vector of row indices identifying age
%                 regressors in the design matrix.
%
%   Output
%     E : [Nmni x 1] age-effect magnitude per voxel.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_ESTIMATE_BETA, ICDM_ESTIMATE_BETA_BAYES
if isscalar(age_index)
    beta_age = squeeze(Beta(age_index,:,:));   % [Nmni × K1]
    E = sqrt(sum(beta_age.^2, 2));                   % effect magnitude
else
    E = icdm_beta_age_effect_map_all(Beta, age_index);
end
end

function E = icdm_beta_age_effect_map_all(Beta, age_basis_indices)

% Beta: [P × Nmni × (K-1)]
% age_basis_indices: e.g. [2 3 4] for 3 spline basis
%
% Output:
%   E(v) = L2 norm of all spline basis effects across ILR dimensions

beta_all = squeeze(Beta(age_basis_indices,:,:));  % [#basis × Nmni × K1]

E = sqrt( sum(beta_all.^2, [1 3]) ).';  % [Nmni × 1]

end
