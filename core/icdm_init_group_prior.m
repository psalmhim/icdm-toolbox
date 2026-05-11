function grp = icdm_init_group_prior(K, outdir, opts)
% ICDM_INIT_GROUP_PRIOR  Initialise group-level prior for iterative EB.
%
%   grp = icdm_init_group_prior(K, outdir, opts)
%
%   Creates the initial group prior structure for the first EB iteration.
%   Prior mean mu is initialised to zero (uninformative), prior precision
%   kappa to a baseline value, and covariate effects Beta to zero.
%
%   OUTPUT fields:
%     grp.mu_ilr_mni : [Nmni x (K-1)] initial group mean (zeros)
%     grp.kappa_mni  : [Nmni x 1]     initial precision (= kappa_base)
%     grp.Beta       : [P x Nmni x (K-1)] initial covariate effects (zeros)
%     grp.H          : [K x (K-1)]    Helmert submatrix
%     grp.idx_mni    : MNI-space voxel indices
%     grp.dim_mni    : MNI volume dimensions
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_POPULATION_EB, HELMERT_SUBMATRIX
Nmni = numel(opts.idx_mni);
K1   = K - 1;
P    = numel(opts.design_spec) - 1;      % slopes only (age)

grp.idx_mni= opts.idx_mni;
grp.dim_mni= opts.dim_mni;
grp.idx_regions = opts.idx_regions;
grp.mu_ilr_mni = zeros(Nmni, K1, 'single');
grp.kappa_mni  = opts.precision.kappa_base * ones(Nmni, 1, 'single');
grp.Beta       = zeros(P, Nmni, K1, 'single');   
grp.H          = helmert_submatrix(K);

save(fullfile(outdir,'group_prior.mat'),'-struct','grp','-v7.3');
end