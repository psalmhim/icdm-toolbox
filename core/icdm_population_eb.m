function icdm_population_eb(subjects, K, outdir, opts, grp)
% ICDM_POPULATION_EB Covariate-free empirical-Bayes population fitting.
%
% Subject estimates and the population prior are iterated without using age,
% sex, or any other covariate. Covariate effects and permutation inference
% are evaluated once, downstream, after the final EB estimates are frozen.

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
n_perm = getfield_default(opts,'n_perm',2000);
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
    resume_inference = false;

    if exist(fniter,'file')
        saved = load(fniter);
        if group_file_complete(saved)
            if ~is_final
                fprintf('[EB] Complete iteration %d found; loading it.\n',it);
                grp = saved;
                continue;
            elseif stats_complete(saved,n_perm)
                fprintf('[EB] Complete final inference found; loading it.\n');
                grp = saved;
                print_stats(grp.stats);
                continue;
            elseif exist(vbfile,'file')
                fprintf('[EB] Final EB fit found; resuming downstream inference.\n');
                tmp = load(vbfile,'VBs');
                if ~isfield(tmp,'VBs')
                    error('Resume file %s does not contain VBs.',vbfile);
                end
                VBs = tmp.VBs;
                grp = saved;
                resume_inference = true;
            else
                fprintf('[EB] Incomplete final inference and no subject checkpoint; rerunning final iteration.\n');
            end
        end
    end

    iter_dir = fullfile(outdir,sprintf('iter_%03d',it));
    if ~exist(iter_dir,'dir'), mkdir(iter_dir); end

    if ~resume_inference
        opts.w_beta = 0; % legacy option explicitly disabled
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
    end

    if is_final
        fprintf('[EB] Estimating downstream covariate effects and sign-flip inference...\n');
        opts_inf = opts;
        opts_inf.w_beta = 0;
        % Default is OLS (beta is estimated downstream); the kappa-weighted 'wls-ridge'
        % estimator can still be requested explicitly.
        opts_inf.estimator = getfield_default(opts,'estimator','ols');
        opts_inf.inference = 'signflip';
        opts_inf.n_perm = n_perm;
        opts_inf.batch_size = getfield_default(opts,'batch_size',256);
        opts_inf.perm_seed = getfield_default(opts,'perm_seed',42);
        opts_inf.do_cluster = getfield_default(opts,'do_cluster',true);
        opts_inf.cluster_alpha = getfield_default(opts,'cluster_alpha',0.05);
        opts_inf.cluster_form_q = getfield_default(opts,'cluster_form_q',0.99);
        opts_inf.cluster_min_size = getfield_default(opts,'cluster_min_size',20);
        opts_inf.fdr_q = getfield_default(opts,'fdr_q',0.05);
        opts_inf.perm_mode = getfield_default(opts,'perm_mode','safe');
        opts_inf.usePar = getfield_default(opts,'usePar',0);
        opts_inf.use_pca = getfield_default(opts,'use_pca',true);
        opts_inf.pca_mode = getfield_default(opts,'pca_mode','auto');
        opts_inf.pca_var_ratio = getfield_default(opts,'pca_var_ratio',0.85);
        opts_inf.pca_max_rank = getfield_default(opts,'pca_max_rank',15);
        opts_inf.pca_min_rank = getfield_default(opts,'pca_min_rank',5);
        opts_inf.grp_kappa_mni = grp.kappa_mni(:);
        % Population precision now reflects inverse between-subject variance,
        % so the legacy count-curvature threshold of 50 is not meaningful.
        opts_inf.grp_kappa_thresh = getfield_default(opts,'grp_kappa_thresh',0);

        t0 = tic;
        [grp.Beta,grp.stats,grp.beta_info] = ...
            icdm_estimate_beta(subjects,VBs,K,opts_inf);
        fprintf('[EB] Downstream inference finished in %.2f seconds.\n',toc(t0));
        grp.inference_config = struct('n_perm',opts_inf.n_perm, ...
            'perm_seed',opts_inf.perm_seed,'estimator',opts_inf.estimator);
        save(fniter,'-struct','grp','-v7.3');
        print_stats(grp.stats);
    end
    fprintf('[EB] Iteration %d complete.\n',it);
end
end


function tf = group_file_complete(g)
required = {'mu_ilr_mni','kappa_mni','tau2','idx_mni','dim_mni','H'};
tf = isfield(g,'model_version') && ...
    strcmp(g.model_version,'covariate_free_eb_v2') && ...
    all(cellfun(@(f) isfield(g,f) && ~isempty(g.(f)),required));
end


function tf = stats_complete(g,n_perm)
tf = isfield(g,'stats') && isstruct(g.stats) && ...
    isfield(g.stats,'T_null_max') && ...
    numel(g.stats.T_null_max) >= n_perm && ...
    all(isfinite(g.stats.T_null_max(:)));
if tf && isfield(g,'inference_config') && isfield(g.inference_config,'n_perm')
    tf = (g.inference_config.n_perm == n_perm);
end
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


function print_stats(stats)
fields = {'T_obs','T_null_max','p_fwer','p_uncorr','q_map'};
labels = {'T_obs','T_null','p_fwer','p_unc','q'};
for j = 1:numel(fields)
    if isfield(stats,fields{j})
        x = double(stats.(fields{j})(:));
        x = x(isfinite(x));
        if ~isempty(x)
            fprintf('%-7s min %.6g | mean %.6g | max %.6g\n', ...
                labels{j},min(x),mean(x),max(x));
        end
    end
end
end
