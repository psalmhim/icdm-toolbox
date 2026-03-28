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

function VBc = icdm_predict_individual_cdm(subj, VBs, grp, opts)
% ICDM_PREDICT_INDIVIDUAL_CDM  Predict individual connectivity via group-individual fusion.
%
%   VBc = icdm_predict_individual_cdm(subj, VBs, grp, opts)
%
%   Produces voxelwise compositional connectivity maps for a single
%   subject by combining group-level covariate-matched predictions with
%   the subject's own posterior estimate.  The three outputs are:
%
%     pi_grp   : predicted connectivity from the group beta-field model,
%                Y_pred(v,:) = Beta(:,v,:)' * x, then ILR-inverse.
%     pi_ind   : subject posterior connectivity from VB inference.
%     pi_fused : kappa-weighted fusion of pi_grp and pi_ind,
%                w_ind = kappa/(kappa+tau), w_grp = tau/(kappa+tau).
%
%   This implements the individual-level prediction step of the iCDM
%   framework (manuscript Section on individual prediction).
%
%   Inputs
%     subj : subject struct (with property fields for design vector).
%     VBs  : subject VB output struct with fields:
%              .y_ilr_mni  [Nmni x (K-1)]  posterior ILR in MNI.
%              .kappa_mni  [Nmni x 1]      posterior precision.
%     grp  : group-level struct with fields:
%              .Beta  [P x Nmni x (K-1)]   regression coefficients.
%              .H     [K x (K-1)]          Helmert sub-matrix.
%     opts : struct with design_spec, idx_mni, fusion_tau.
%
%   Output
%     VBc : struct; VBc.pred contains pi_grp, pi_ind, pi_fused, Y_pred.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_ESTIMATE_BETA, ICDM_SUBJECT_VB, ICDM_DESIGN_VECTOR

fprintf('=== Predicting Individual Brain Connectivity ===\n');

Nmni = numel(opts.idx_mni);
H    = grp.H;
K    = size(H,1);
K1   = K - 1;

%% --------------------------------------------------------------
% 1. Make design vector for this subject (supports spline)
%% --------------------------------------------------------------
x = icdm_design_vector(subj, opts);      % row vector [1 × P]
x = double(x(:));                        % column vector [P × 1]
P = numel(x);

fprintf('  Design vector length = %d (supports spline)\n', P);

%% --------------------------------------------------------------
% 2. Compute voxelwise predicted ILR from group Beta
%      Y_pred(v,:) = (Beta(:,v,:))' * x
%% --------------------------------------------------------------
fprintf('  Computing voxelwise ILR prediction...\n');

Beta = grp.Beta;   % [P × Nmni × K1]
if size(Beta,1) ~= P
    error('grp.Beta first dimension (%d) does not match design length (%d).', size(Beta,1), P);
end

Y_pred = zeros(Nmni, K1, 'single');

for v = 1:Nmni
    Bv = squeeze(Beta(:,v,:));    % [P × K1]
    yv = Bv.' * x;                % [K1 × 1]
    Y_pred(v,:) = single(yv);
end

%% --------------------------------------------------------------
% 3. Convert group prediction ILR ? probabilities
%% --------------------------------------------------------------
fprintf('  Converting group prediction ILR ? probability...\n');

PI_grp = zeros(Nmni, K, 'single');

for v = 1:Nmni
    y = double(Y_pred(v,:)).';
    z = H * y;
    ez = exp(z);
    PI_grp(v,:) = single(ez ./ sum(ez));
end

%% --------------------------------------------------------------
% 4. Convert individual's posterior ILR ? probabilities
%% --------------------------------------------------------------
fprintf('  Converting individual posterior ILR ? probability...\n');

Y_ind = VBs.y_ilr_mni;                   % [Nmni × K1]
PI_ind = zeros(Nmni, K, 'single');

for v = 1:Nmni
    y = double(Y_ind(v,:)).';
    z = H * y;
    ez = exp(z);
    PI_ind(v,:) = single(ez ./ sum(ez));
end

%% --------------------------------------------------------------
% 5. ?-based voxelwise fusion
%       w_ind = ? / (?+?)
%       w_grp = ? / (?+?)
%   supports spline automatically
%% --------------------------------------------------------------
tau = getfield_default(opts, 'fusion_tau', 1.0);

kappa = VBs.kappa_mni(:);          % [Nmni × 1]
kappa = max(kappa, eps('single'));
if tau/kappa<0.1
    tau=median(kappa); 
    fprintf('tau=%f is determined by median kappa\n',tau);
end
w_ind = kappa ./ (kappa + tau);
w_grp = tau   ./ (kappa + tau);

PI_fused = PI_ind .* w_ind + PI_grp .* w_grp;

%% --------------------------------------------------------------
% 6. Store result
%% --------------------------------------------------------------
VBc.pred.pi_grp   = PI_grp;
VBc.pred.pi_ind   = PI_ind;
VBc.pred.pi_fused = PI_fused;
VBc.pred.Y_pred   = Y_pred;

fprintf('=== Prediction Completed ===\n');

end
