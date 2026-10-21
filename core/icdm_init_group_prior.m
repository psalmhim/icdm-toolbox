function grp = icdm_init_group_prior(K, outdir, opts)
% =========================================================================
% icdm_init_group_prior.m
%
% 인덱스 기반 초기 group prior 생성
%
%   grp.mu_ilr_mni : [Nmni x (K-1)] 초기 0
%   grp.kappa_mni  : [Nmni x 1]      = kappa_base
%   grp.H          : [K x (K-1)] Helmert matrix
%
% (disk에 group_prior.mat으로도 저장)
% =========================================================================
Nmni = numel(opts.idx_mni);
K1   = K - 1;

grp.idx_mni= opts.idx_mni;
grp.dim_mni= opts.dim_mni;
grp.idx_regions = opts.idx_regions;
grp.mu_ilr_mni = zeros(Nmni, K1, 'single');
grp.kappa_mni  = opts.precision.kappa_base * ones(Nmni, 1, 'single');
grp.H          = helmert_submatrix(K);

save(fullfile(outdir,'group_prior.mat'),'-struct','grp','-v7.3');
end