function Beta = icdm_beta_estimator(X, Yall, Kap, coverage, opts)

[S, Nmni, K1] = size(Yall);
P = size(X,2);

Beta = zeros(P, Nmni, K1, 'single');

alpha  = getfield_default(opts,'beta_alpha',0.5);
gamma  = getfield_default(opts,'beta_gamma',0.2);
lambda = getfield_default(opts,'ridge_lambda',0);

Kap_alpha = max(Kap, eps).^alpha;
w_cov_all = max(coverage, eps).^gamma;

for v = 1:Nmni

    Yv = reshape(Yall(:,v,:), [S K1]);

    switch lower(opts.estimator)

        case 'ols'
            XtX = X'*X;
            B = XtX \ (X'*Yv);

        case 'wls'
            w  = Kap_alpha(:,v) * w_cov_all(v);
            sw = sqrt(w);
            Xsw = X .* sw;
            Ysw = Yv .* sw;
            B = (Xsw'*Xsw) \ (Xsw'*Ysw);

        case 'wls-ridge'
            w  = Kap_alpha(:,v) * w_cov_all(v);
            sw = sqrt(w);
            Xsw = X .* sw;
            Ysw = Yv .* sw;
            XtX = Xsw'*Xsw + lambda*eye(P);
            B = XtX \ (Xsw'*Ysw);

        otherwise
            error('Unknown estimator.');
    end

    Beta(:,v,:) = single(B);
end

end