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

function x = icdm_design_vector(subj, opts)
% ICDM_DESIGN_VECTOR  Construct a subject's design-matrix row from the design spec.
%
%   x = icdm_design_vector(subj, opts)
%
%   Evaluates the design specification (opts.design_spec) for a single
%   subject and returns a numeric row vector.  Supported entry types:
%     - 'const'       : inserts 1 (intercept).
%     - character name: retrieves subj.property.<name>, z-scores
%                       continuous variables using opts.design.mu / sd.
%     - struct with fields var/K : evaluates a natural cubic spline basis
%                       for the named variable with K knots.
%
%   Inputs
%     subj : subject struct with subj.property.<name> fields.
%     opts : options struct containing design_spec, design.type,
%            design.mu, and design.sd (prepared by
%            icdm_prepare_design_vector).
%
%   Output
%     x : [1 x P] numeric design-matrix row for this subject.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_PREPARE_DESIGN_VECTOR, ICDM_ESTIMATE_BETA

spec = opts.design_spec;
x_list = {};

for i = 1:numel(spec)
    entry = spec{i};

    % ------------------ constant --------------------
    if ischar(entry) && strcmpi(entry,'const')
        x_list{end+1} = 1;
        continue;
    end

    % ------------------ normal variable --------------------
    if ischar(entry)
        name = entry;
        val = subj.property.(name);

        if strcmp(opts.design.type.(name),'continuous')
            val = (val - opts.design.mu.(name)) / opts.design.sd.(name);
        end

        x_list{end+1} = val;
        continue;
    end

    % ------------------ spline struct --------------------
    if isstruct(entry)
        basevar = entry.var;
        Kspl    = entry.K;

        fname = sprintf('%s_spline_%d', basevar, Kspl);
        fname = matlab.lang.makeValidName(fname);

        % z-scored base value
        v = subj.property.(basevar);
        v = (v - opts.design.mu.(basevar)) / opts.design.sd.(basevar);

        % build natural cubic spline basis
        bx = icdm_make_spline_basis(v, Kspl);

        x_list = [x_list, num2cell(bx)]; 
        continue;
    end
end

x = cell2mat(x_list);

end

function b = icdm_make_spline_basis(x, K)
% Natural cubic spline basis with K knots on [-2,2] (z-scored domain)
knots = linspace(-2,2,K);
b = zeros(1,K);
for i = 1:K
    b(i) = max(0, x - knots(i)).^3;
end
end
