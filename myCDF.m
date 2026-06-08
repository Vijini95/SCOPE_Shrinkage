function F = myCDF(x, distName, params)
% myCDF  Compute CDF for a selected distribution by name.
%
% Usage examples:
%   F = myCDF(x, 'normal',      [mu sigma])
%   F = myCDF(x, 'cauchy',      [x0 gamma])
%   F = myCDF(x, 'logistic',    [mu s])
%   F = myCDF(x, 'exponential', [lambda])
%
% Notes on parameters:
%   - normal:      [mu, sigma]          sigma > 0
%   - cauchy:      [x0, gamma]          gamma > 0   (location, scale)
%   - logistic:    [mu, s]              s > 0       (location, scale)
%   - exponential: [lambda]             lambda > 0  (rate)
%   - uniform:     [a, b]               a < b
%   - gamma:       [k, theta]           k>0, theta>0 (shape, scale)
%   - beta:        [a, b]               a>0, b>0
%   - lognormal:   [mu, sigma]          sigma>0 (log-space params)
%   - poisson:     [lambda]             lambda>=0
%   - binomial:    [n, p]               n>=0 integer, 0<=p<=1
%   - t:           [nu]                 nu > 0
%   - chisquare:   [nu]                 nu > 0

if nargin < 3
    params = [];
end

% Dictionary of CDFs
CDF = struct();

% Continuous
CDF.normal      = @(x,p) normcdf(x, p(1), p(2));
CDF.exponential = @(x,p) expcdf(x, 1/p(1));          % p(1)=lambda (rate)
CDF.uniform     = @(x,p) unifcdf(x, p(1), p(2));
CDF.gamma       = @(x,p) gamcdf(x, p(1), p(2));
CDF.beta        = @(x,p) betacdf(x, p(1), p(2));
CDF.lognormal   = @(x,p) logncdf(x, p(1), p(2));
CDF.t           = @(x,p) tcdf(x, p(1));
CDF.chisquare   = @(x,p) chi2cdf(x, p(1));

% New: Cauchy and Logistic (manual formulas; no toolbox dependency)
CDF.cauchy      = @(x,p) (1/pi) .* atan((x - p(1))./p(2)) + 0.5;   % [x0 gamma]
CDF.logistic    = @(x,p) 1 ./ (1 + exp(-(x - p(1))./p(2)));        % [mu s]

% Discrete
CDF.poisson     = @(x,p) poisscdf(x, p(1));
CDF.binomial    = @(x,p) binocdf(x, p(1), p(2));

% Evaluate
distName = lower(strrep(distName,'-','')); % allows "chi-square" and "chisquare"

if ~isfield(CDF, distName)
    error('Distribution "%s" not found.', distName);
end

% Optional: basic parameter check (lightweight)
validateParams(distName, params);

F = CDF.(distName)(x, params);

end

