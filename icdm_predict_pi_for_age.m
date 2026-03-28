% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI", NeuroImage (2026).
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

function PI4 = icdm_predict_pi_for_age(age0, Beta, H, opts, sex0)
% ICDM_PREDICT_PI_FOR_AGE  Predict a voxelwise compositional map for a given age.
%
%   PI4 = icdm_predict_pi_for_age(age0, Beta, H, opts, sex0)
%
%   Constructs a design vector for the specified age (and optionally sex),
%   applies the group beta-field model at every MNI voxel to predict the
%   ILR coordinates, and converts to the probability simplex via the
%   Helmert ILR-inverse transform.  Supports both simple linear and
%   spline-based age regressors as specified in opts.design_spec.
%
%   Inputs
%     age0 : target age (scalar, in original units).
%     Beta : [P x Nmni x K1] group-level regression coefficients.
%     H    : [K x K1] Helmert sub-matrix for ILR inverse.
%     opts : struct with design_spec, design.mu, design.sd, idx_mni,
%            dim_mni.
%     sex0 : sex covariate value (default 0).
%
%   Output
%     PI4 : [X x Y x Z x K] single 4-D volume of predicted compositional
%           probabilities in MNI space.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PREDICT_INDIVIDUAL_CDM, ICDM_ESTIMATE_BETA, ILR_INVERSE

if nargin < 5
    sex0 = 0; % default = 0
end

% ------------------------------------------------------------
% 0) Basic sizes (from Beta, not from spec)
% ------------------------------------------------------------
P     = size(Beta, 1);   % #predictors
Nmni  = size(Beta, 2);   % #voxels
K1    = size(Beta, 3);   % ILR dim (K-1)
K     = size(H,1);       % #targets, H: [K x K1]

% Safety check (optional)
if size(H,2) ~= K1
    error('icdm_predict_pi_for_age: H second dim (%d) ~= K1 from Beta (%d).', size(H,2), K1);
end

% ------------------------------------------------------------
% 1) Build design vector x(age0, sex0) of length P
% ------------------------------------------------------------
spec = opts.design_spec;

x = zeros(P,1);

for j = 1:length(spec)

    entry = spec{j};

    % Case A: simple name ('const','age','sex',...)
    if ischar(entry) || isstring(entry)
        name = char(entry);

        switch lower(name)
            case 'const'
                x(j) = 1;

            case 'age'
                mu = opts.design.mu.age;
                sd = opts.design.sd.age;
                if sd == 0, sd = 1; end
                x(j) = (age0 - mu)/sd;

            case 'sex'
                x(j) = sex0;

            otherwise
                x(j) = 0;
        end

    % Case B: struct specification
    %   e.g., struct('var','age','type','spline','K',3,'idx',k)
    elseif isstruct(entry)

        base  = '';
        btype = '';
        Kb    = 1;

        if isfield(entry,'var'),  base  = entry.var;  end
        if isfield(entry,'type'), btype = entry.type; end
        if isfield(entry,'K'),    Kb    = entry.K;    end

        if strcmpi(base,'age') && strcmpi(btype,'spline')
            mu = opts.design.mu.age;
            sd = opts.design.sd.age;
            if sd == 0, sd = 1; end
            a  = (age0 - mu)/sd;

            b = icdm_make_spline_basis(a, Kb);   % 1×K

            if isfield(entry,'idx')
                idx = entry.idx;
                if idx >= 1 && idx <= numel(b)
                    x(j) = b(idx);
                else
                    x(j) = 0;
                end
            else
                x(j:j+Kb-1) = b(1:Kb);
            end
        else
            x(j) = 0;
        end

    else
        x(j) = 0;
    end
end

x = double(x(:));   % P×1

% ------------------------------------------------------------
% 2) Predict ILR for every voxel
%    Beta: [P × Nmni × K1]
%    For each voxel v:
%       Bv = [P × K1]
%       y  = Bv.' * x   → [K1 × 1]
% ------------------------------------------------------------
Yhat = zeros(Nmni, K1, 'single');   % ILR predictions per voxel

for v = 1:Nmni
    Bv = squeeze(Beta(:, v, :));   % [P × K1]
    % ensure 2D
    if isvector(Bv)
        Bv = reshape(Bv, [P, K1]);
    end
    yv = x'*Bv;                 % [K1 × 1]
    Yhat(v,:) = single(yv(:));   % row 1×K1
end

% ------------------------------------------------------------
% 3) Convert ILR → simplex using Helmert matrix H (K×K1)
% ------------------------------------------------------------
dim = opts.dim_mni;
Xn = dim(1); Yn = dim(2); Zn = dim(3);

PI4 = zeros(Xn,Yn,Zn,K,'single');

for i = 1:Nmni
    y = double(Yhat(i,:)).';   % [K1 × 1]

    % ILR inverse: z = H * y  (log-ratio→log)
    z = H * y;                 % [K × 1]

    % stable softmax
    z = z - max(z);
    a = exp(z);
    s = sum(a);
    if s == 0
        pi = ones(K,1) / K;
    else
        pi = a / s;
    end

    lin = opts.idx_mni(i);
    [xi, yi, zi] = ind2sub([Xn Yn Zn], lin);
    PI4(xi,yi,zi,:) = single(pi);
end

end


% ================================================================
% Natural cubic spline basis (simple truncated-power basis)
% ================================================================
function b = icdm_make_spline_basis(x, K)
% Domain assumed z-scored: x in [-2,2] approx

knots = linspace(-2,2,K);
b = zeros(1,K);

for i = 1:K
    b(i) = max(0, x - knots(i)).^3;
end
end
