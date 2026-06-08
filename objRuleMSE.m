function mseAvg = objRuleMSE(lambda, k, yTrue, Nrep, coarsest, filt, distName, distParams)
% objRuleMSE Average MSE over Nrep noise realizations for given lambda,k.

% Hard safety (keeps optimizer stable)
if ~isfinite(lambda) || ~isfinite(k) || lambda <= 0 || k <= 0
    mseAvg = inf;
    return;
end

n = numel(yTrue);
mse = zeros(Nrep,1);

for r = 1:Nrep
    noise  = randn(1,n);      % sigma = 1
    yNoisy = yTrue + noise;

    wt_data = dwtr(yNoisy, coarsest, filt);
    wt_hat  = ScopeRule(wt_data, distName, distParams, lambda, k);
    y_hat   = idwtr(wt_hat, coarsest, filt);

    mse(r) = mean((y_hat - yTrue).^2);
end

mseAvg = mean(mse);
end
