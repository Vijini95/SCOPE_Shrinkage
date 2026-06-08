close all force; clear; clc

addpath '/Users/hvimalajeewa2/Documents/TAMU_Documents/TAMU/DemosNew'
addpath '/Users/hvimalajeewa2/Documents/TAMU_Documents/TAMU/WaveletSrinkage/Codes /codesMatLab5.2'
addpath './MlFunctions/'

%% Fixed setup
randn('seed',5);

J = 10; 
n = 2^J;
t = linspace(0,1,n);

L = 4;
SNR = 5;

sigma  = 1;
sigma2 = sigma^2;

%% Select signal and wavelet filter
for sType = 4
   switch sType 
       case 1
           fun='Wave';      
           wtype='Symmlet'; 
           filtersize = 8;
       case 2
           fun='Blocks';    
           wtype='Haar';    
           filtersize = 2;
       case 3
           fun='HeaviSine'; 
           wtype='Symmlet'; 
           filtersize = 8;
       case 4
           fun='Doppler';   
           wtype='Symmlet'; 
           filtersize = 8;
       case 5
           fun='Bumps';     
           wtype='Symmlet'; 
           filtersize = 4;
   end
end   

filt = MakeONFilter(wtype, filtersize);

coarsest = J - L;

%% True signal only for simulation check
yTrue0 = MakeSignal(fun, n);
yTrue  = sqrt(SNR)/std(yTrue0) * yTrue0;

%% Generate one noisy observation
noise  = sigma * randn(1,n);
yNoisy = yTrue + noise;

%% Wavelet transform
wt_data = dwtr(yNoisy, coarsest, filt);

% Wavelab-style indexing:
% scaling coefficients: 1:2^coarsest
% detail coefficients:  2^coarsest+1:n
nScaling = 2^coarsest;
detailIdx = (nScaling+1):n;

%% Distribution for F
distName   = 'normal';
distParams = [0 1];

%% SURE grid search
lambdaGrid = linspace(0.01, 2.0, 40);
kGrid      = linspace(0.10, 20.0, 40);

SUREsurf = zeros(numel(lambdaGrid), numel(kGrid));

bestSURE = inf;
bestLam  = NaN;
bestK    = NaN;

fprintf('SURE grid search: %d x %d = %d evaluations...\n', ...
    numel(lambdaGrid), numel(kGrid), numel(lambdaGrid)*numel(kGrid));

for i = 1:numel(lambdaGrid)

    lam = lambdaGrid(i);

    for j = 1:numel(kGrid)

        kk = kGrid(j);

        sureVal = objScopeSURE(lam, kk, wt_data(detailIdx), ...
                               sigma2, distName, distParams);

        SUREsurf(i,j) = sureVal;

        if sureVal < bestSURE
            bestSURE = sureVal;
            bestLam  = lam;
            bestK    = kk;
        end

    end
end

fprintf('Best grid SURE: lambda = %.6f, k = %.6f, SURE = %.6g\n', ...
    bestLam, bestK, bestSURE);

%% Local refinement using fminsearch
u0 = log([bestLam, bestK]);

opts = optimset('Display','iter', ...
                'MaxFunEvals', 150, ...
                'MaxIter', 150, ...
                'TolX',1e-5, ...
                'TolFun',1e-6);

uHat = fminsearch(@(u) objScopeSURE(exp(u(1)), exp(u(2)), ...
                                    wt_data(detailIdx), ...
                                    sigma2, distName, distParams), ...
                  u0, opts);

lambdaHat = exp(uHat(1));
kHat      = exp(uHat(2));

finalSURE = objScopeSURE(lambdaHat, kHat, wt_data(detailIdx), ...
                         sigma2, distName, distParams);

fprintf('\nOptimized by SURE: lambda = %.6f, k = %.6f, SURE = %.6g\n', ...
    lambdaHat, kHat, finalSURE);

%% Apply SURE-tuned SCOPE shrinkage
wt_hat = wt_data;

wt_hat(detailIdx) = ScopeRuleSUREVersion(wt_data(detailIdx), ...
                                         distName, distParams, ...
                                         lambdaHat, kHat);

% Do not shrink scaling coefficients
wt_hat(1:nScaling) = wt_data(1:nScaling);

y_hat = idwtr(wt_hat, coarsest, filt);

MSE_one = mean((y_hat - yTrue).^2);

fprintf('One-draw MSE after SURE tuning = %.6g\n', MSE_one);

%% Plot SURE surface
figure('Renderer','painters','Position',[200 200 650 600]);

[KK, LL] = meshgrid(kGrid, lambdaGrid);

contour(KK, LL, SUREsurf, 20, 'LineWidth', 2);
hold on

mask = SUREsurf <= 1.05 * bestSURE;
plot(KK(mask), LL(mask), 'r.', 'MarkerSize', 10);

plot(kHat, lambdaHat, 'ko', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('k');
ylabel('\lambda');
title('Near-Optimal SURE Region');
grid on

%% Plot reconstruction
figure('Renderer','painters','Position',[200 200 1200 700]);

subplot(2,1,1)
plot(t, yNoisy, '.', 'MarkerSize', 8);
hold on
plot(t, y_hat, 'LineWidth', 1.3);
plot(t, yTrue, 'LineWidth', 1.3);
legend('Noisy','SURE estimate','True');
grid on
title(sprintf('%s | SURE-SCOPE | lambda=%.4f, k=%.4f | MSE=%.4g', ...
    fun, lambdaHat, kHat, MSE_one));

subplot(2,1,2)
plot(wt_data, 'LineWidth', 1.0);
hold on
plot(wt_hat, 'LineWidth', 1.0);
legend('WT noisy','WT thresholded');
grid on
title('Wavelet coefficients')
