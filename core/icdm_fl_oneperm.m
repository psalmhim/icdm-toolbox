function [Tmax, Tperm_all, maxCluster] = icdm_fl_oneperm( ...
    X, Yhat_red, R_red, Kap_alpha, w_cov_all, age_idx, lambda, dof, ...
    batch_size, doClusterNull, dim_mni, idx_mni, cluster_form_q, idx)
% ICDM_FL_ONEPERM  Standalone (parfor-safe) copy of the Freedman-Lane one-permutation
% statistic used in icdm_estimate_beta. Identical computation; exists as a file
% function so parfor workers can resolve it (a local subfunction cannot be called
% reliably from a parallel pool). The serial inference path still uses the in-file
% subfunction, so this copy does not change the standard result.

[S,Nmni,K1] = size(Yhat_red); %#ok<ASGLU>
P = size(X,2);

if nargin < 14 || isempty(idx)
    idx = randperm(S);
end

Tmax = 0;
maxCluster = 0;
Tperm_all = zeros(Nmni,1,'single');

nbatch = ceil(Nmni / batch_size);
for b = 1:nbatch
    v_start = (b-1)*batch_size + 1;
    v_end   = min(b*batch_size, Nmni);
    for v = v_start:v_end
        Yh = double(squeeze(Yhat_red(:,v,:)));   % [S×K1]
        Rv = double(squeeze(R_red(:,v,:)));      % [S×K1]

        % Freedman–Lane: permute reduced residuals
        Yperm = Yh + Rv(idx,:);

        % voxel weights
        w  = double(Kap_alpha(:,v) * w_cov_all(v));
        sw = sqrt(w);

        % full model refit
        Xsw = X .* sw;
        Ysw = Yperm .* sw;
        XtX = (Xsw'*Xsw) + lambda*eye(P);
        Bp  = XtX \ (Xsw'*Ysw);

        % studentized T
        R = Yperm - X*Bp;
        WRSS = sum((R.^2).*w,1);
        sigma2 = WRSS ./ dof;

        XtX_inv = XtX \ eye(P);
        var_age = sigma2(:) .* XtX_inv(age_idx,age_idx);

        beta_age = Bp(age_idx,:).';
        Tperm = sum(beta_age.^2 ./ var_age);

        Tperm_all(v) = single(Tperm);
        if Tperm > Tmax, Tmax = Tperm; end
    end
end

if doClusterNull
    thr = quantile(double(Tperm_all), cluster_form_q);
    Tvol = nan(dim_mni,'single');
    Tvol(idx_mni) = Tperm_all;
    mask = Tvol > thr;
    CC = bwconncomp(mask,26);
    if CC.NumObjects > 0
        maxCluster = max(cellfun(@numel, CC.PixelIdxList));
    else
        maxCluster = 0;
    end
end
end
