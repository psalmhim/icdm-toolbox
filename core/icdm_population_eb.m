function icdm_population_eb(subjects, K, outdir, opts, grp)
% ICDM_POPULATION_EB Covariate-free empirical-Bayes population fitting.
%
% Iteration 1 estimates the population prior from the count-derived subject
% representations; iteration 2 estimates every subject under that fixed prior.
% No covariate enters. Outputs: group_icdm_iter_XXX.mat (prior), iter_XXX/<id>_vb.mat
% (subject estimates) and subject_vb_iter_<final>.mat. Association with covariates is
% analysed downstream on the saved estimates.

fprintf('\n=== ICDM Population EB (covariate-free prior) ===\n');
if ~exist(outdir,'dir'), mkdir(outdir); end
S = numel(subjects);
if S < 3, error('At least three subjects are required.'); end
if opts.max_iter < 2
    error(['opts.max_iter must be at least 2: iteration 1 estimates the ' ...
        'data-derived population prior and iteration 2 applies that prior.']);
end

log_file = fullfile(outdir,'vb_log.txt');
diary(log_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>
opts.idx_mni = grp.idx_mni;
opts.dim_mni = grp.dim_mni;
stale_fields = intersect({'Beta','stats','beta_info','inference_config'},fieldnames(grp));
if ~isempty(stale_fields), grp = rmfield(grp,stale_fields); end
grp.model_version = 'covariate_free_eb_v2';

for it = 1:opts.max_iter
    fprintf('\n[EB] Iteration %d/%d\n',it,opts.max_iter);
    fniter = fullfile(outdir,sprintf('group_icdm_iter_%03d.mat',it));
    vbfile = fullfile(outdir,sprintf('subject_vb_iter_%03d.mat',it));
    is_final = (it == opts.max_iter);

    if exist(fniter,'file')
        saved = load(fniter);
        if group_file_complete(saved) && (~is_final || exist(vbfile,'file'))
            fprintf('[EB] Complete iteration %d found; loading it.\n',it);
            grp = saved;
            continue;
        end
    end

    iter_dir = fullfile(outdir,sprintf('iter_%03d',it));
    if ~exist(iter_dir,'dir'), mkdir(iter_dir); end

    opts.w_beta = 0;
    fprintf('[EB] Running subject-wise MAP/Laplace for %d subjects...\n',S);
    t0 = tic;
    VBs = icdm_run_subject_vb(subjects,K,grp,opts,iter_dir);
    fprintf('[EB] Subject estimation finished in %.2f seconds.\n',toc(t0));

    old_mu = grp.mu_ilr_mni;
    old_kappa = grp.kappa_mni;
    fprintf('[EB] Updating covariate-free population prior...\n');
    t0 = tic;
    [grp.mu_ilr_mni,grp.kappa_mni,grp.tau2] = ...
        icdm_update_group_prior(VBs,opts,it);
    fprintf('[EB] Population update finished in %.2f seconds.\n',toc(t0));

    [dmu,dkappa] = parameter_change(old_mu,old_kappa, ...
        grp.mu_ilr_mni,grp.kappa_mni);
    grp.eb_change = struct('mu_relative',dmu,'kappa_relative',dkappa);
    fprintf('[EB] Relative change: mu %.6g | kappa %.6g\n',dmu,dkappa);

    save(fniter,'-struct','grp','-v7.3');
    if is_final
        save(vbfile,'VBs','-v7.3');
    end
    fprintf('[EB] Saved EB iteration to %s\n',fniter);

    fprintf('[EB] Iteration %d complete.\n',it);
end
end


function tf = group_file_complete(g)
required = {'mu_ilr_mni','kappa_mni','tau2','idx_mni','dim_mni','H'};
tf = isfield(g,'model_version') && ...
    strcmp(g.model_version,'covariate_free_eb_v2') && ...
    all(cellfun(@(f) isfield(g,f) && ~isempty(g.(f)),required));
end


function [dmu,dkappa] = parameter_change(mu0,k0,mu1,k1)
if isempty(mu0) || ~isequal(size(mu0),size(mu1))
    dmu = Inf;
else
    good = isfinite(mu0) & isfinite(mu1);
    dmu = norm(double(mu1(good)-mu0(good))) / max(norm(double(mu0(good))),eps);
end
if isempty(k0) || ~isequal(size(k0),size(k1))
    dkappa = Inf;
else
    good = isfinite(k0) & isfinite(k1);
    dkappa = norm(double(k1(good)-k0(good))) / max(norm(double(k0(good))),eps);
end
end
