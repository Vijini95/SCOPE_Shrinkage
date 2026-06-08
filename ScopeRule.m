function R = ScopeRule(x, distName, distParams, lambda, k)
% myRule  Compute Rule(x) = x * |F*(lambda*x)|^k
% where F*(p) = 2*F(p) - 1 and F is a CDF from myCDF.
%
% Usage:
%   R = myRule(x, 'normal',   [0 1],   2.0, 3);
%   R = myRule(x, 'cauchy',   [0 1],   1.5, 2);
%   R = myRule(x, 'logistic', [0 1],   0.8, 1);
%
% Inputs:
%   x          : scalar or vector
%   distName   : distribution name string (as used in myCDF)
%   distParams : parameter vector for that distribution
%   lambda     : scalar parameter
%   k          : scalar parameter (typically >= 0)
%
% Output:
%   R          : Rule(x) values (same size as x)

% ---- basic checks ----
if nargin < 5
    error('Usage: R = myRule(x, distName, distParams, lambda, k)');
end
if ~isscalar(lambda) || ~isfinite(lambda)
    error('lambda must be a finite scalar.');
end
if ~isscalar(k) || ~isfinite(k)
    error('k must be a finite scalar.');
end

% ---- compute p = lambda*x ----
p = lambda .* x;

% ---- compute F(p) using myCDF ----
Fval = myCDF(p, distName, distParams);

% ---- compute F*(p) = 2F(p) - 1 ----
Fstar = 2 .* Fval - 1;

% ---- Rule(x) = x * |Fstar|^k ----
R = x .* (abs(Fstar) .^ k);

end
