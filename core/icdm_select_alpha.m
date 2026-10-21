function A = icdm_select_alpha(vb_files, grp, opts)
% ICDM_SELECT_ALPHA  Age-independent calibration of the population-transfer strength alpha_grp.
%
%   A = icdm_select_alpha(vb_files, grp, opts)
%
% Uses only the prior-independent count-derived representations (y_data_mni, kappa_data_mni) of the
% reference sample; no covariate (age, group, ...) enters, so the effect later tested cannot steer alpha.
%   data term   : y_data ~ N(m, kappa_data^-1 I)
%   prior       : m ~ N(mu_grp, (alpha * lambda_grp)^-1 I)
%   marginal    : y_data ~ N(mu_grp, (kappa_data^-1 + (alpha*lambda_grp)^-1) I)
% Criteria
%   A.alpha_ml : in-sample type-II marginal likelihood
%   A.alpha_cv : K-fold cross-validated held-out predictive likelihood, with the fold prior
%                (trimmed mean, lambda = 1/max(mean variance, 1e-6)) re-estimated without the held-out subjects
%
% Inputs
%   vb_files : cell array of <id>_vb.mat files (icdm_subject_vb outputs) for the reference sample
%   grp      : population prior (struct or file) with idx_mni, mu_ilr_mni, kappa_mni
%   opts     : .grid (default logspace(-2,1.6,80)), .nfold (5), .trim (0.20), .seed (7)
%              .groups (optional, one label per file, e.g. participant ID): folds are assigned per group,
%              so repeated visits of one participant never fall on both sides of a fold
if nargin < 3, opts = struct; end
grid  = getfield_default(opts,'grid',logspace(-2,1.6,80));
nfold = getfield_default(opts,'nfold',5);
trim  = getfield_default(opts,'trim',0.20);
seed  = getfield_default(opts,'seed',7);
if ischar(grp), grp = load(grp,'idx_mni','mu_ilr_mni','kappa_mni'); end
S = numel(vb_files); Nv = numel(grp.idx_mni); K1 = size(grp.mu_ilr_mni,2);

valid = true(Nv,1);
for s = 1:S, d = load(vb_files{s},'valid_mni'); valid = valid & logical(d.valid_mni(:)); end
vi = find(valid); nv = numel(vi);
mu = double(grp.mu_ilr_mni(vi,:)); kg = double(grp.kappa_mni(vi));
YD = zeros(S,nv,K1,'single'); KD = zeros(S,nv,'single');
for s = 1:S
    d = load(vb_files{s},'y_data_mni','kappa_data_mni');
    yy = single(d.y_data_mni); YD(s,:,:) = reshape(yy(vi,:),1,nv,K1);
    KD(s,:) = single(d.kappa_data_mni(vi))';
end
KDv = max(double(KD),1e-9);

R2 = zeros(S,nv);
for s = 1:S, dif = double(squeeze(YD(s,:,:))) - mu; R2(s,:) = sum(dif.^2,2)'; end
kgR = repmat(kg',S,1);
nA = numel(grid); mll = zeros(1,nA);
for i = 1:nA
    sig2 = 1./KDv + 1./(grid(i)*kgR);
    mll(i) = sum(sum(-0.5*K1*log(2*pi*sig2) - 0.5*R2./sig2));
end

rng(seed);
if isfield(opts,'groups') && ~isempty(opts.groups)
    [~,~,gid] = unique(opts.groups(:)); gfold = mod(randperm(max(gid)),nfold) + 1; fold = gfold(gid)';
else
    fold = mod(randperm(S),nfold) + 1;
end
cvll = zeros(1,nA);
for k = 1:nfold
    tr = find(fold ~= k); ho = find(fold == k); ntr = numel(tr); nt = floor(trim*ntr);
    mu_tr = zeros(nv,K1); tau_tr = zeros(nv,K1);
    for dd = 1:K1
        col = sort(double(YD(tr,:,dd)),1); col = col(nt+1:ntr-nt,:);
        mu_tr(:,dd) = mean(col,1)'; tau_tr(:,dd) = var(col,0,1)';
    end
    kg_tr = 1./max(mean(tau_tr,2),1e-6);
    for j = ho(:)'
        dif = double(squeeze(YD(j,:,:))) - mu_tr; r2 = sum(dif.^2,2)'; kd = max(double(KD(j,:)),1e-9);
        for i = 1:nA
            sig2 = 1./kd + 1./(grid(i)*kg_tr');
            cvll(i) = cvll(i) + sum(-0.5*K1*log(2*pi*sig2) - 0.5*r2./sig2);
        end
    end
end

A.grid = grid; A.mll = mll; A.cvll = cvll; A.nv = nv;
A.alpha_ml = refine(grid, mll);
A.alpha_cv = refine(grid, cvll);
end

function a = refine(grid, f)
% grid argmax with parabolic refinement on log(alpha)
[~,i] = max(f); a = grid(i);
if i > 1 && i < numel(grid)
    la = log(grid(i-1:i+1)); y = f(i-1:i+1); den = y(1) - 2*y(2) + y(3);
    if den ~= 0, a = exp(la(2) + 0.5*(y(1)-y(3))/den*(la(3)-la(2))); end
end
end
