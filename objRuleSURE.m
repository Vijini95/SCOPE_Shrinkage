function sureVal = objRuleSURE(lambda, k, wt_data, detailIdx, sigma2, distName, distParams)
% objScopeSURE computes the SURE criterion:
%
% Rhat(lambda,k)
% =
% sum_i { (y_i A_i - y_i)^2 - sigma^2 }
% +
% 2 sigma^2 sum_i { A_i + k lambda |y_i| |F*(lambda y_i)|^(k-1) f*(lambda y_i) }
%
% where
% A_i = |F*(lambda y_i)|^k,
% F*(u) = 2F(u)-1,
% f*(u) = 2f(u).

    if ~isfinite(lambda) || ~isfinite(k) || lambda <= 0 || k <= 0
        sureVal = inf;
        return;
    end

    y = wt_data(detailIdx);
    y = y(:);

    u = lambda * y;

    [Fstar, fstar] = centeredCDFandDensity(u, distName, distParams);

    absF = abs(Fstar);

    % Numerical protection near Fstar = 0
    eps0 = 1e-12;
    absFsafe = max(absF, eps0);

    A = absFsafe.^k;

    shrinkTerm = (y .* A - y).^2 - sigma2;

    divergenceTerm = A + k * lambda .* abs(y) .* absFsafe.^(k-1) .* fstar;

    sureTerms = shrinkTerm + 2 * sigma2 * divergenceTerm;

    sureVal = sum(sureTerms);

    if ~isfinite(sureVal)
        sureVal = inf;
    end
end
