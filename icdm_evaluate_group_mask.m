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

function gopts = icdm_evaluate_group_mask(subjects, min_coverage_frac, outfile)
% ICDM_EVALUATE_GROUP_MASK  Compute the group-level MNI analysis mask from subject coverage.
%
%   gopts = icdm_evaluate_group_mask(subjects, min_coverage_frac, outfile)
%
%   For each subject the native WM mask is warped forward to MNI space and
%   a coverage count is accumulated.  Voxels where at least
%   min_coverage_frac of all subjects contribute data are retained in the
%   group analysis mask.  The result is cached to outfile so that
%   subsequent calls skip the expensive warping step.
%
%   Inputs
%     subjects          : struct array from icdm_compose_subject.
%     min_coverage_frac : minimum fraction of subjects required
%                         (e.g. 0.3 = 30 %).  Default 0.3.
%     outfile           : path to save/load the group mask MAT.
%
%   Output
%     gopts : struct with fields
%               .dim_mni  [1 x 3] MNI grid dimensions.
%               .idx_mni  uint32 linear indices of included voxels.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_COMPOSE_SUBJECT, ICDM_POPULATION_EB
if nargin<3, outfile = ''; end
if nargin<2, min_coverage_frac = 0.3; end
S = numel(subjects);
if S == 0
    error('No subjects given.');
end

if exist(outfile,'file')
    fprintf('[GroupMask] Loading existing group mask file: %s\n', outfile);
    L = load(outfile, 'dim_mni','idx_mni');
    gopts.dim_mni = L.dim_mni;
    gopts.idx_mni = L.idx_mni;
    return;
end

fprintf('[GroupMask] Evaluating group MNI coverage over %d subjects...\n', S);
% Add each subject's MNI mask index to coverage count
cov_count=[];
for s=1:S
    if rem(s,10)==1
        fprintf('  Processing subject %d/%d: %s\n', s, S, subjects(s).id);
    end
    info = load(subjects(s).datafile);
    M  = zeros(info.dim_native,"single");
    M(info.idx_native) = 1;
    % native -> MNI (nearest neighbor)
    Yout = icdm_warp_4d(M, info.V, subjects(s).dartel_flow, subjects(s).template, +1, 0);
    dim_mni = size(Yout);
    if numel(dim_mni) < 3
        dim_mni(3) = 1;
    end
    idx_mask = find(Yout > 0.5);
    if s==1
        Nvox_mni = prod(dim_mni);
        cov_count = zeros(Nvox_mni,1,'uint16');
    end
    cov_count(idx_mask) = cov_count(idx_mask) + 1;
end

idx_mni = uint32(find(cov_count / S >= min_coverage_frac));

gopts.dim_mni      = dim_mni;
gopts.idx_mni      = uint32(idx_mni);

fprintf('[GroupMask] dim=[%d %d %d], coverage voxels=%d (%.2f%%)\n', ...
    dim_mni(1),dim_mni(2),dim_mni(3), ...
    numel(idx_mni), 100*numel(idx_mni)/Nvox_mni);

if ~isempty(outfile)
    save(outfile,'idx_mni','dim_mni', '-v7.3');
    fprintf('[GroupMask] Saved group mask info to %s\n', outfile);
end
end

