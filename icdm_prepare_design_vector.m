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

function opts = icdm_prepare_design_vector(subjects, opts)
% ICDM_PREPARE_DESIGN_VECTOR  Pre-compute design metadata (z-scoring, types, splines).
%
%   opts = icdm_prepare_design_vector(subjects, opts)
%
%   Inspects the design specification in opts.design_spec and prepares the
%   normalisation parameters needed by icdm_design_vector.  For each
%   continuous variable the sample mean and standard deviation are
%   computed across all subjects and stored in opts.design.mu / sd.
%   Variable types (continuous / categorical) are inferred automatically
%   or from naming conventions.  Spline entries store knot metadata.
%
%   Example design_spec:
%     { 'const', 'age', struct('var','age','type','spline','K',3), 'sex' }
%
%   Inputs
%     subjects : [S x 1] struct array with property fields.
%     opts     : struct containing opts.design_spec.
%
%   Output
%     opts : updated struct with opts.design.mu, opts.design.sd,
%            opts.design.type, opts.design.spline.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_DESIGN_VECTOR, BUILD_AGE_SPLINE_BASIS

spec = opts.design_spec;
S = numel(subjects);

opts.design.mu   = struct();
opts.design.sd   = struct();
opts.design.type = struct();
opts.design.spline = struct();

for i = 1:numel(spec)

    entry = spec{i};

    % ============================================================
    % Case 1: STRING (normal variable)
    % ============================================================
    if ischar(entry) || isstring(entry)

        name = char(entry);

        if strcmpi(name, 'const')
            opts.design.type.const = 'const';
            continue;
        end

        % Infer variable type
        lowername = lower(name);
        if contains(lowername,'age')
            vtype = 'continuous';
        elseif contains(lowername,'sex') || contains(lowername,'gender')
            vtype = 'categorical';
        else
            % auto detect
            vec = arrayfun(@(s) s.property.(name), subjects);
            if numel(unique(vec)) <= 5
                vtype = 'categorical';
            else
                vtype = 'continuous';
            end
        end

        % Store type
        opts.design.type.(name) = vtype;

        % Z-scoring
        if strcmp(vtype,'continuous')
            vec = arrayfun(@(s) s.property.(name), subjects);
            m = mean(vec);
            d = std(vec); if d==0, d=1; end
            opts.design.mu.(name) = m;
            opts.design.sd.(name) = d;
        else
            opts.design.mu.(name) = 0;
            opts.design.sd.(name) = 1;
        end

        continue;
    end


    % ============================================================
    % Case 2: STRUCT  → spline or GAM term
    %   struct('var','age','type','spline','K',3)
    % ============================================================
    if isstruct(entry)

        if ~isfield(entry,'var') || ~isfield(entry,'type')
            error('Spline spec must contain fields "var" and "type".');
        end

        basevar = entry.var;
        vtype   = entry.type;

        % e.g., "age_spline3"
        fname = sprintf('%s_%s_%d', basevar, vtype, entry.K);
        fname = matlab.lang.makeValidName(fname);

        if strcmp(vtype,'spline')

            opts.design.type.(fname) = 'spline';
            opts.design.spline.(fname).var = basevar;
            opts.design.spline.(fname).K   = entry.K;

            % base variable must be continuous
            vec = arrayfun(@(s) s.property.(basevar), subjects);
            m = mean(vec); d = std(vec); if d==0, d=1; end
            opts.design.mu.(basevar) = m;
            opts.design.sd.(basevar) = d;

        else
            error('Unsupported design struct type: %s', vtype);
        end

        continue;
    end

    error('Invalid design_spec entry: must be string or struct.');
end

end
