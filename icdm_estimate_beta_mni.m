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

function beta_idx_mat = icdm_estimate_beta_mni(group_path, subjects, istat_path, K, opts, save_nii)
% ICDM_ESTIMATE_BETA_MNI  Voxelwise weighted ridge regression for ILR beta in MNI space.
%
%   beta_idx_mat = icdm_estimate_beta_mni(group_path, subjects,
%                                          istat_path, K, opts, save_nii)
%
%   Estimates per-voxel regression coefficients (Beta) in MNI space using
%   weighted ridge regression.  Subject-level ILR means are loaded from
%   *_mni_mu.mat files, kappa weights from *_mni_kappa.mat or native
%   subj_stats.mat (warped forward if needed).  The ridge regularization
%   and kappa-tempering produce stable estimates even with incomplete
%   coverage.
%
%   Inputs
%     group_path  : folder containing *_mni_mu.mat, mni_mask_fix.mat.
%     subjects    : struct array (fields icdm_4d, dartel_flow, template).
%     istat_path  : folder with subj_***/subj_stats.mat (native kappa).
%     K           : number of compositional targets (labels).
%     opts        : struct with optional fields:
%         .design_matrix      [S x P] explicit design override.
%         .design_spec        cell, e.g. {'const','age','age2'}.
%         .ridge_lambda       regularization (default 1e-3).
%         .beta_weight_alpha  kappa tempering exponent (default 0.5).
%         .beta_coverage_gamma coverage tempering (default 0.25).
%     save_nii    : logical; write NIfTI volumes of Beta (default false).
%
%   Output
%     beta_idx_mat : path to saved MAT file containing
%                    Beta [P x Nmni x (K-1)], idx_mni, dim_mni.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_ESTIMATE_BETA, ICDM_ESTIMATE_BETA_BAYES
%
% -------------------------------------------------------------------------

if nargin < 6 || isempty(save_nii), save_nii = false; end
if nargin < 5, opts = struct; end

beta_idx_mat = fullfile(group_path, 'beta_ilr_mni.mat');
if exist(beta_idx_mat,'file')==2
    if save_nii
        beta_nii = fullfile(group_path, 'beta_ilr_mni.nii');
        if ~exist(beta_nii,'file')
            Lb   = load(beta_idx_mat,'Beta','idx_mni','dim_mni');
            Vtpl = pick_template(subjects);
            materialize_beta_nii_from_idx(Vtpl,Lb.Beta, Lb.idx_mni,beta_nii);
        end
    end
    return;
end

% -------------------- Load mask / geometry --------------------
Lmask = load(fullfile(group_path,'mni_mask_fix.mat'),'idx_mni','dim_mni');
idx_mni = Lmask.idx_mni; dim_mni = double(Lmask.dim_mni);
nI = numel(idx_mni);
Km1 = K - 1;
S   = numel(subjects);

% -------------------- Design matrix --------------------
if isfield(opts,'design_matrix') && ~isempty(opts.design_matrix)
    X = opts.design_matrix;             % S×P
else
    spec = getfield_default(opts,'design_spec', {'const','age','age2'});%,'sex','age*sex'});
    X = zeros(S, numel(spec));
    for s = 1:S
        X(s,:) = make_design_vector_flex(subjects(s), spec);
    end
end
[Schk,P] = size(X);
assert(Schk==S, 'Design matrix must be S×P with S subjects.');

Xd = double(X);

% -------------------- Collect ILR Y (per subject, in MNI) --------------------
H = helmert_submatrix(K);
Yall = zeros(S, nI, Km1, 'single');

mats = dir(fullfile(group_path,'*_mni_mu.mat'));
assert(numel(mats)==S, 'Expected %d *_mni_mu.mat under %s', S, group_path);

for s = 1:S
    L  = load(fullfile(mats(s).folder, mats(s).name), 'mu');
    MU = reshape(renorm_mu(L.mu), [], K);
    clr = log(max(MU, realmin('single'))) - mean(log(max(MU, realmin('single'))),2);
    Y   = clr * H;  % [nI × Km1]
    Yall(s,:,:) = single(Y);
end

% -------------------- Build κ weights in MNI --------------------
Wkappa = ones(S, nI, 'single');
for s = 1:S
    subj=subjects(s);
    [~,subjname]=fileparts(fileparts(subj.icdm_4d));

    kap_mni_mat  = fullfile(group_path, sprintf('%s_mni_kappa.mat', subjname));

    subdir = fullfile(istat_path, subjname);
    stats_path = fullfile(subdir, 'subj_stats.mat');
    if exist(kap_mni_mat,'file')==2
        Ls = load(kap_mni_mat,'kappa_mni');
        if isfield(Ls,'kappa_mni') && ~isempty(Ls.kappa_mni)
            v = reshape(max(single(Ls.kappa_mni(:,:,:,1)),0), [], 1);
            Wkappa(s,:) = v;
            continue;
        end
    else
        if exist(stats_path,'file')==2
            Ls = load(stats_path,'dim','idx','kappa');
            if isfield(Ls,'kappa') && ~isempty(Ls.kappa)
                % scatter native κ -> 3D, warp to MNI, sample at idx_mni
                kap3d = zeros(prod(Ls.dim),1,'single'); kap3d(Ls.idx) = single(Ls.kappa);
                kap3d = reshape(kap3d, Ls.dim);
                kap4d = reshape(kap3d, [Ls.dim 1]);

                Vref = spm_vol(subjects(s).icdm_4d); Vref=Vref(1);
                kap_mni4d = generic_warp_4d(kap4d, Vref, ...
                            getfield_default(subjects(s),'dartel_flow',''), ...
                            getfield_default(subjects(s),'template',''), +1, 1);
                v = reshape(max(single(kap_mni4d(:,:,:,1)),0), [], 1);
                Wkappa(s,:) = v(idx_mni);
            end
        end
    end
end

% Optional: temper by coverage (raises weights where more subjects cover)
if ~isempty(coverage_mni_last)
    cov_all = max(min(reshape(single(opts.coverage_mni), [],1),1),0);
    gamma_cov = getfield_default(opts,'beta_coverage_gamma', 0.25);
    if gamma_cov ~= 0
        Wkappa = bsxfun(@times, Wkappa, (cov_all'.^gamma_cov));
    end
end

% Temper and normalize per voxel
alpha_w = getfield_default(opts,'beta_weight_alpha', 0.5);
Wkappa = max(Wkappa, eps('single')).^alpha_w;
medW   = median(Wkappa, 1);
medW   = max(medW, eps('single'));
Wkappa = bsxfun(@rdivide, Wkappa, medW);

% -------------------- Weighted ridge per voxel & ILR dim --------------------
rl = getfield_default(opts,'ridge_lambda', 1e-3);
Ireg = eye(P);
Beta = zeros(P, nI, Km1, 'single');

dq = parallel.pool.DataQueue;
progress=0;Niter=nI;h = waitbar(0,'Starting...');
afterEach(dq, @updateProgress);
tic
for d = 1:Km1
    Yd = Yall(:,:,d);   % [S × nI]
    progress=0;Niter=nI;%h = waitbar(0,'Starting...');
    parfor i = 1:nI
        send(dq, 1);
        y = double(Yd(:,i));
        w = double(Wkappa(:,i));
        msk = isfinite(y) & isfinite(w) & (w > 0);
        if nnz(msk) < P
            % fallbacks for underdetermined voxel
            if nnz(msk) >= P
                Xm = Xd(msk,:); ym = y(msk);
                B  = (Xm.'*Xm + rl*Ireg) \ (Xm.'*ym);
            else
                B  = zeros(P,1);
            end
            Beta(:,i,d) = single(B); continue;
        end
        Xm = Xd(msk,:); wm = w(msk); ym = y(msk);
        sw  = sqrt(wm);
        Xsw = bsxfun(@times, Xm, sw);
        ysw = sw .* ym;
        B   = (Xsw.'*Xsw + rl*Ireg) \ (Xsw.'*ysw);
        Beta(:,i,d) = single(B);
    end
end

% -------------------- Save compact MAT; optional NIfTI --------------------
Vtpl = pick_template(subjects);
save(beta_idx_mat, 'Beta','idx_mni','dim_mni','K','P','Vtpl','-v7.3');
toc;

if save_nii
    beta_nii = fullfile(group_path,'beta_ilr_mni.nii');
    materialize_beta_nii_from_idx(Beta, idx_mni, dim_mni, Vtpl, beta_nii);
end


function updateProgress(~)
    progress = progress + 1;
    if rem(progress,round(Niter/100))~=1, return; end
    waitbar(progress/nI,h, ...
            sprintf('Progress: %d / %d', progress, Niter));
end

end % ===== end main =====


% =================== Helpers (self-contained) ===================

function H = helmert_submatrix(K)
H = zeros(K, K-1);
for i=1:(K-1)
    H(1:i, i) =  1 / sqrt(i*(i+1));
    H(i+1, i) = -i / sqrt(i*(i+1));
end
end


function mu = renorm_mu(mu_in)
mu = max(mu_in, 0);
sz = size(mu); K = sz(2);
sums = sum(mu,2);
zero_vox = (sums <= 0);
mu(zero_vox, :) = 1/K;
sums(~zero_vox) = 1 ./ sums(~zero_vox);
mu(~zero_vox,:) = mu(~zero_vox,:) .* sums(~zero_vox);
end

function mu4d = renorm_mu4d(mu4d_in)
mu = max(mu4d_in, 0);
sz = size(mu); K = sz(4);
mu = reshape(mu, [], K);
sums = sum(mu,2);
zero_vox = (sums <= 0);
mu(zero_vox, :) = 1/K;
sums(~zero_vox) = 1 ./ sums(~zero_vox);
mu(~zero_vox,:) = mu(~zero_vox,:) .* sums(~zero_vox);
mu4d = reshape(mu, sz);
end

function Vtpl = pick_template(subjects)
iTpl = find(arrayfun(@(s) isfield(s,'template') && ~isempty(s.template), subjects), 1);
assert(~isempty(iTpl), 'No subject.template found for MNI header.');
Vtpl = spm_vol(subjects(iTpl).template); Vtpl = Vtpl(1);
end

function materialize_beta_nii_from_idx(Vtpl, Beta, idx_mni,out_path)
[P, ~, Km1] = size(Beta);
[p1,f1,e1]=fileparts(out_path);
for p=1:P
    out_path1 = fullfile(p1,sprintf('%s_beta%d.nii',f1,p));
    Vout = repmat(Vtpl,Km1, 1);
    for f=1:Km1
        Vout(f).fname   = out_path1;
        Vout(f).n       = [f 1];
        Vout(f).dt      = [spm_type('float32') 0];
        Vout(f).pinfo   = [1;0;0];
        Vout(f).descrip = sprintf('β (ILR) p=%d, d=%d', p, f);
        vol = zeros(Vout(f).dim(1:3),'single');
        vol(idx_mni) = single(Beta(p,:,f));
        spm_write_vol(Vout(f), reshape(vol,Vout(f).dim(1:3)));
    end
end
end


function val = getfield_default(S, field, defaultVal)
if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
    val = S.(field);
else
    val = defaultVal;
end
end

function x = make_design_vector_flex(subj, spec)
if nargin<2 || isempty(spec)
    spec = {'const','age','age2'};%,'sex','age*sex'};
end
P = numel(spec); x = zeros(1,P);
age = getfield_default(subj,'age',0);
sex = mapsex(getfield_default(subj,'sex',0));
for t=1:P
    tok = lower(spec{t});
    switch tok
        case 'const',    v = 1;
        case 'age',      v = age;
        case 'age2',     v = age.^2;
        case 'sex',      v = sex;
        case {'age*sex','sex*age'}, v = age*sex;
        otherwise
            if isfield(subj, spec{t}), v = scalar(subj.(spec{t}));
            else, v = 0;
            end
    end
    x(t) = scalar(v);
end
end

function s = mapsex(v)
if ischar(v) || isstring(v)
    vs = lower(string(v));
    if any(vs == ["m","male","1"]), s = 1; return; end
    if any(vs == ["f","female","0","2"]), s = 0; return; end
    d = str2double(vs); if isnan(d), v = 0; else, v = d; end
end
v = scalar(v); if v==2, v=0; end
s = double(v >= 0.5);
end

function v = scalar(v)
if isempty(v), v=0; return; end
if iscell(v), v=v{1}; end
if ischar(v) || isstring(v)
    d = str2double(string(v)); if isnan(d), v=0; else, v=d; end
end
v = double(v); if ~isscalar(v), v=v(1); end
if ~isfinite(v), v=0; end
end
