close all force;
clear;
clc;

set(0,'DefaultFigureVisible','on');

%% ==========================================================
%% PATHS (HPRC)
%% ==========================================================

addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\WaveletSrinkage\LPM\matlab');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\WaveletSrinkage\NewCodes');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\WaveletSrinkage\Codes\codesMatLab5_2');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\wavelab850\Orthogonal');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\wavelab850\Utilities');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\wavelab850\DeNoising');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes\WaveletSrinkage\Block_shrinkage');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Dixon_Codes');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\viji');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Scope_Codes');
addpath('C:\Users\vijin\OneDrive\Desktop\SCOPE\Scope_Codes\MlFunctions');

outroot = 'C:\Users\vijin\OneDrive\Desktop\Wavelets\Project_3_SCOPE\SCOPE_FINAL\scope_output2\';
if ~exist(outroot,'dir'), mkdir(outroot); end

rng(5);

%% ==========================================================
%% SETTINGS
%% ==========================================================

signals   = {'Doppler','HeaviSine','Bumps','Blocks'};
SNR_list  = [3 5 7];
cdf_list  = {'logistic','normal','laplace','sech','uniform','cauchy'};

J = 10;
n = 2^J;
L = 5;
Nrep = 100;

lambdaGrid = linspace(0.01,2,100);
kGrid      = linspace(1,20,100);

%% ==========================================================
%% TABLE 1 COMPUTATION (CORRECT WAVELET INSIDE LOOP)
%% ==========================================================

Results = {};

for s = 1:length(signals)
for sn = 1:length(SNR_list)

    signalName = signals{s};
    SNR = SNR_list(sn);

    %% --- SELECT WAVELET PER SIGNAL ---
    switch signalName
        case 'Blocks'
            wtype = 'Haar'; filtersize = 2;
        case 'Bumps'
            wtype = 'Daubechies'; filtersize = 6;
        otherwise
            wtype = 'Symmlet'; filtersize = 8;
    end

    filt = MakeONFilter(wtype, filtersize);
    coarsest = J - L;

    y0 = MakeSignal(signalName,n);
    yTrue = sqrt(SNR)/std(y0) * y0;

    for c = 1:length(cdf_list)

    distName = cdf_list{c};

    % ---- Correct parameters per distribution ----
    switch distName
        case 'uniform'
            distParams = [-1 1];
        otherwise
            distParams = [0 1];
    end

    bestMSE = inf; 
    bestLam = NaN; 
    bestK   = NaN;

    for i = 1:length(lambdaGrid)
    for j = 1:length(kGrid)

        lam = lambdaGrid(i);
        kk  = kGrid(j);

        mse = objRuleMSE(lam,kk,yTrue,Nrep,...
                         coarsest,filt,distName,distParams);

        if mse < bestMSE
            bestMSE = mse;
            bestLam = lam;
            bestK   = kk;
        end
    end
    end

    Results(end+1,:) = {signalName,SNR,distName,...
                        bestMSE,bestLam,bestK};
end
end
end

ResultsTable = cell2table(Results,...
    'VariableNames',{'Signal','SNR','CDF','AMSE','Lambda','K'});

writetable(ResultsTable,fullfile(outroot,'Table1_results.csv'));
save(fullfile(outroot,'Table1_results.mat'),'ResultsTable');

fprintf('Table 1 saved correctly.\n');

%% ==========================================================
%% CDF GROUPED BOXPLOT (SNR = 5)
%% ==========================================================

SNR = 5;
AllData = [];
GroupSignal = [];
GroupCDF = [];

for s = 1:length(signals)

    signalName = signals{s};

    %% Wavelet again (important)
    switch signalName
        case 'Blocks'
            wtype = 'Haar'; filtersize = 2;
        case 'Bumps'
            wtype = 'Daubechies'; filtersize = 6;
        otherwise
            wtype = 'Symmlet'; filtersize = 8;
    end

    filt = MakeONFilter(wtype, filtersize);
    coarsest = J - L;

    idx = strcmp(ResultsTable.Signal,signalName) & ...
          ResultsTable.SNR==5;

    subTab = ResultsTable(idx,:);

    y0 = MakeSignal(signalName,n);
    yTrue = sqrt(SNR)/std(y0) * y0;

    for c = 1:length(cdf_list)

        distName = cdf_list{c};
        row = strcmp(subTab.CDF,distName);

        lambdaHat = subTab.Lambda(row);
        kHat      = subTab.K(row);

        for rep = 1:Nrep

            y = yTrue + randn(1,n);

            wt_data = dwtr(y,coarsest,filt);
            switch distName
              case 'uniform'
              distParams = [-1 1];  % or [-sqrt(3) sqrt(3)]
              otherwise
              distParams = [0 1];
             end

            wt_hat  = ScopeRule(wt_data,distName,distParams,...
                    lambdaHat,kHat);
            yscope  = idwtr(wt_hat,coarsest,filt);

            mse = mean((yTrue-yscope).^2);

            AllData = [AllData; mse];
            GroupSignal = [GroupSignal; s];
            GroupCDF = [GroupCDF; c];
        end
    end
end

%% ==========================================================
%% IMPROVED GROUPED BOXPLOT WITH SIGNAL SPACING
%% ==========================================================

figure('Color','w','Position',[300 300 1000 550])
hold on

gap = 2.3;      % space between signals (increase if needed)
within = 0.30;  % spacing inside each signal

colors = lines(length(cdf_list));

for s = 1:length(signals)

    base = (s-1)*gap;

    for c = 1:length(cdf_list)

        idx = (GroupSignal==s & GroupCDF==c);
        data = AllData(idx);

        xpos = base + (c-3)*within;

        boxchart(ones(size(data))*xpos, data,...
            'BoxWidth',0.22,...
            'BoxFaceColor',colors(c,:),...
            'MarkerStyle','o');
    end
end

xticks((0:length(signals)-1)*gap)
xticklabels(signals)

set(gca,'FontSize',16)
ylabel('MSE','FontSize',16)

grid on
box on

legend(cdf_list,'Location','northoutside',...
       'Orientation','horizontal','FontSize',14)

exportgraphics(gcf,...
    fullfile(outroot,'CDF_Boxplot_SNR5.png'),...
    'Resolution',300);
close

fprintf('CDF grouped boxplot saved (with spacing).\n');

%% ==========================================================
%% METHOD COMPARISON (DIRECT 2x2 ONLY)
%% ==========================================================

MethodNames = {'SCOPE','BAMS','Decompsh','BlockMed',...
               'BlockMean','Hybrid','BlockJS','VisuShrink','GCV'};

for sn = 1:length(SNR_list)

    SNR = SNR_list(sn);
    AllMethodMSE = cell(length(signals),1);
    BestCDF_perSignal = cell(length(signals),1);

    for s = 1:length(signals)

        signalName = signals{s};

        %% ---- SELECT WAVELET ----
        switch signalName
            case 'Blocks'
                wtype = 'Haar'; filtersize = 2;
            case 'Bumps'
                wtype = 'Daubechies'; filtersize = 6;
            otherwise
                wtype = 'Symmlet'; filtersize = 8;
        end

        filt = MakeONFilter(wtype, filtersize);
        coarsest = J - L;

        %% ---- Get best lambda,k ----
        idx = strcmp(ResultsTable.Signal,signalName) & ...
              ResultsTable.SNR==SNR;

        subTab = ResultsTable(idx,:);
        [~,ind] = min(subTab.AMSE);

        bestCDF   = subTab.CDF{ind};
        lambdaHat = subTab.Lambda(ind);
        kHat      = subTab.K(ind);

        BestCDF_perSignal{s} = bestCDF;

        %% ---- Generate signal ----
        y0 = MakeSignal(signalName,n);
        yTrue = sqrt(SNR)/std(y0) * y0;

        MSE = zeros(9,Nrep);

        for rep = 1:Nrep

            y = yTrue + randn(1,n);

            wt_data = dwtr(y,coarsest,filt);
            switch bestCDF
               case 'uniform'
               distParams = [-1 1];
               otherwise
               distParams = [0 1];
             end

            wt_hat = ScopeRule(wt_data,bestCDF,distParams,...
                   lambdaHat,kHat);
            yscope  = idwtr(wt_hat,coarsest,filt);
            MSE(1,rep)=mean((yTrue-yscope).^2);

            ybams = BAMS(y,filt,7);
            MSE(2,rep)=mean((yTrue-ybams(:)').^2);

            yCV = recdecompsh(y',filt);
            MSE(3,rep)=mean((yTrue-yCV(:)').^2);

            thet0b =[0.5 1 0.7];
            MSE(4,rep)=mean((yTrue-recblockmed('Augment',y,[],thet0b)).^2);
            MSE(5,rep)=mean((yTrue-recblockmean('Augment',y,[],thet0b)).^2);

            theta0s = [0.5 1];
            MSE(6,rep)=mean((yTrue-rechybblockmed('Augment',y,[],theta0s,thet0b)).^2);

            MSE(7,rep)=mean((yTrue-recblockJS('Augment',y,filt)).^2);
            MSE(8,rep)=mean((yTrue-recvisu(y,'H',filt)).^2);
            MSE(9,rep)=mean((yTrue-recgcv(y,filt)).^2);
        end

        AllMethodMSE{s} = MSE;
    end

    %% ======================================================
    %% 2x2 FIGURE ONLY
    %% ======================================================

    figure('Color','w','Position',[200 200 1000 800])

    for s = 1:length(signals)

        subplot(2,2,s)
        boxplot(AllMethodMSE{s}')
        xticks(1:9)
        xticklabels(1:9)
        set(gca,'FontSize',15)
        ylabel('MSE')
        grid on
        box on

        % Title = Signal (bestCDF)
        title([signals{s} ' (' BestCDF_perSignal{s} ')'],...
              'FontSize',16)
    end

    exportgraphics(gcf,...
        fullfile(outroot,...
        ['Method_Boxplot_2x2_SNR' num2str(SNR) '.png']),...
        'Resolution',300);

    close
end

fprintf('2x2 method comparison figures saved.\n');