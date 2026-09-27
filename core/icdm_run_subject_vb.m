function OUTS = icdm_run_subject_vb(subjects, K, grp, opts, iter_dir)
% Run subject-level VB for each subject (with cache).
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
