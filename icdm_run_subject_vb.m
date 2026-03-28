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

function OUTS = icdm_run_subject_vb(subjects, K, grp, opts, iter_dir)
% ICDM_RUN_SUBJECT_VB  Execute subject-level variational Bayes for all subjects.
%
%   OUTS = icdm_run_subject_vb(subjects, K, grp, opts, iter_dir)
%
%   Iterates over the subject array and runs icdm_subject_vb for each
%   subject, caching the result as <sid>_vb.mat in iter_dir.  On
%   subsequent calls cached results are loaded directly, making the
%   function safe for interrupted runs.  This constitutes the E-step of
%   the iCDM population EB algorithm.
%
%   Inputs
%     subjects : [S x 1] struct array from icdm_compose_subject.
%     K        : number of compositional targets.
%     grp      : group-level prior struct (from icdm_update_group_prior).
%     opts     : VB and design options struct.
%     iter_dir : directory where per-subject VB cache files are stored.
%
%   Output
%     OUTS : {S x 1} cell array of subject VB output structs.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_SUBJECT_VB, ICDM_POPULATION_EB
S = numel(subjects);
OUTS = cell(S,1);

for s = 1:S
    subj = subjects(s);
    sid  = subj.id;
    fprintf(' [VB] Subject %d/%d: %s\n', s, S, sid);
    cache_file = fullfile(iter_dir, [sid '_vb.mat']);
    % Try loading cache
    if ~exist(cache_file,'file')
        try
            tic
            OUT = icdm_subject_vb(subj, K, grp, opts);
            elapsed_time = toc;
            save(cache_file, '-struct', 'OUT', '-v7.3');
            fprintf('  Computed VB: %s (%.2f seconds)\n', sid, elapsed_time);
        catch
            fprintf('ERROR!!! in %s\n',sid);
        end
    end
end

for s = 1:S
    subj = subjects(s);
    sid  = subj.id;
    fprintf(' [VB] Subject %d/%d: %s\n', s, S, sid);
    cache_file = fullfile(iter_dir, [sid '_vb.mat']);
    need_recalc = true;

    % Try loading cache
    if exist(cache_file,'file')
        try
            OUTS{s} = load(cache_file);
            fprintf('  Loaded cached VB: %s\n', sid);
            need_recalc = false;
        catch
            warning('  Cache load failed. Recomputing: %s', sid);
        end
    end

    % Compute VB
    if need_recalc
        try
            tic
            OUT = icdm_subject_vb(subj, K, grp, opts);
            elapsed_time = toc;
            save(cache_file, '-struct', 'OUT', '-v7.3');
            OUTS{s} = OUT;
            fprintf('  Computed VB: %s (%.2f seconds)\n', sid, elapsed_time);
        catch
            fprintf('ERROR!!! in %s\n',sid);
        end
    end
end
end
