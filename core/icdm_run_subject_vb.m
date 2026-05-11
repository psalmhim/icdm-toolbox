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
