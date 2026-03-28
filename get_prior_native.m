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

function pred_mni = get_prior_native(grp, subj, opts, idx_native)
% GET_PRIOR_NATIVE  Compute covariate-predicted prior mean in MNI space.
%
%   pred_mni = get_prior_native(grp, subj, opts, idx_native)
%
%   Computes the covariate-informed prior mean for a given subject by
%   multiplying the group-level regression coefficients Beta with the
%   subject's covariate vector (manuscript Eq. subject_prior):
%
%     m_v^(s) = (1-zeta) * mu_grp + zeta * Beta' * x^(s)
%
%   This function returns the Beta' * x^(s) component in MNI index space.
%   The blending with mu_grp is done in icdm_subject_vb.
%
%   INPUT
%     grp        : group prior struct with .Beta [P x Nmni x K1]
%     subj       : subject struct (used for design vector construction)
%     opts       : options struct (design specification)
%     idx_native : (optional) native-space voxel indices
%
%   OUTPUT
%     pred_mni : [Nmni x K1] predicted ILR prior from covariates
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_SUBJECT_VB, ICDM_ESTIMATE_BETA, ICDM_DESIGN_VECTOR
if nargin < 4
    idx_native = [];
end
if isempty(idx_native) && isfield(subj,'idx_native') && ~isempty(subj.idx_native)
    idx_native = subj.idx_native;
end
pred_mni = [];
% --------------------------------------------------------------
% 0) Check Beta field
% --------------------------------------------------------------
if ~isfield(grp,'Beta') || isempty(grp.Beta)
    return;
end

grp.Beta = grp.Beta;           % [P × Nmni × K1]
idx_mni    = grp.idx_mni(:);
dim_mni    = double(grp.dim_mni);
Nmni       = numel(idx_mni);
K1         = size(grp.Beta,3);

% --------------------------------------------------------------
% 1) build subject covariate vector
% --------------------------------------------------------------
try
    x = icdm_design_vector(subj, opts);   % [1 × P]
catch
    return;
end

if any(~isfinite(x))
    return;
end

% --------------------------------------------------------------
% 2) predict β in MNI index space
% --------------------------------------------------------------
pred_mni = zeros(Nmni, K1, 'single');
for d = 1:K1
    Bd = double(squeeze(grp.Beta(:,:,d)));      % [P × Nmni]
    pred_mni(:,d) = single( (x * Bd).' );
end

end
