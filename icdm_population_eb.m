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

function icdm_population_eb(subjects, K, outdir, opts, grp)
% ICDM_POPULATION_EB  Iterative empirical Bayes loop for hierarchical iCDM.
%
%   icdm_population_eb(subjects, K, outdir, opts, grp)
%
%   Implements the alternating E-M scheme described in the manuscript
%   (Section: Iterative empirical Bayes inference):
%     E-step: per-subject Newton-Laplace MAP inference (icdm_subject_vb)
%     M-step: robust group aggregation of ILR posteriors
%             (icdm_update_group_prior), followed by covariate regression
%             (icdm_estimate_beta)
%
%   Each iteration refines group-level priors (mu_grp, kappa_grp, tau2,
%   Beta) and propagates them back to subject-level inference.
%
%   INPUT
%     subjects : [S x 1] struct array with per-subject data paths
%     K        : number of compositional components
%     outdir   : output directory for iteration results
%     opts     : options struct (.max_iter, .w_beta, .vb, .precision, ...)
%     grp      : initial group prior (from icdm_init_group_prior)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_SUBJECT_VB, ICDM_UPDATE_GROUP_PRIOR, ICDM_ESTIMATE_BETA

fprintf('\n=== ICDM Population EB (FINAL) ===\n');

S = numel(subjects);
orig_w_beta = opts.w_beta;

log_file = fullfile(outdir, 'vb_log.txt');
diary(log_file);

for it = 1:opts.max_iter
    fprintf('\n[EB] Iteration %d/%d\n', it, opts.max_iter);
    fniter=fullfile(outdir,sprintf('group_icdm_iter_%03d.mat', it));
    if exist(fniter,'file')
        fprintf('[EB] Iteration %d already done. Skipping...\n', it);
        g = load(fniter);
        req = {'mu_ilr_mni','kappa_mni','tau2','Beta','idx_mni','dim_mni','H'};
        for f = req
            if ~isfield(g,f{1})
                error('Saved grp for iter %d is missing field %s', it, f{1});
            end
        end
        grp = g;
        continue;
    end
    iter_dir = fullfile(outdir, sprintf('iter_%03d',it));
    if ~exist(iter_dir,'dir'), mkdir(iter_dir); end

    if it == 1
        opts.w_beta = 0;
    else
        opts.w_beta = orig_w_beta;
    end

    %% SUBJECT VB
    fprintf('[EB] Running subject wise VB of %d subjects...\n', S);
    tic
    VBs = icdm_run_subject_vb(subjects, K, grp, opts, iter_dir);
    elapsed_time=toc;
    fprintf('[EB] Subject wise VB done taking %g seconds.\n', elapsed_time);

    %% GROUP UPDATE
    fprintf('[EB] Updating group parameters...\n');
    tic
    [grp.mu_ilr_mni, grp.kappa_mni, grp.tau2] = icdm_update_group_prior(VBs,opts,it);
    elapsed_time=toc;
    fprintf('[EB] Group parameters updated taking %g seconds.\n', elapsed_time);

    %% BETA UPDATE
    fprintf('[EB] Updating group beta...\n');
    tic
    grp.Beta = icdm_estimate_beta(subjects, VBs, K, opts);
    elapsed_time=toc;
    if opts.verbose > 0
        fprintf('[EB] Group beta updated (iter %03d) taking %g seconds.\n', it, elapsed_time);
    end

    %% SAVE
    fprintf('[EB] Saving iteration results...\n');
    save(fniter, '-struct','grp','-v7.3');

    fprintf('[EB] Iteration %d DONE.\n', it);
end

end
