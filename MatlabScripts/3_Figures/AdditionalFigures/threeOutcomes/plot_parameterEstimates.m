clear; clc
username = getenv('USERNAME');
task = 'threeOutcomes';

load(['C:\Users\' username '\Dropbox\Research\collaborations\ottenheimerDavid\' ...
      'fits\' task '\' task '_MLEfits.mat'])
myColors = importColors_bb;
VP_color = myColors.bluishGreen;
saveLoc = ['C:\Users\' username '\Dropbox\Research\collaborations\ottenheimerDavid\' ...
           'figures\revision\v2\' task];

%% get relevant behavior models
modelCriterion = 'AIC';
plotFlag = false;

models_of_interest_RPE = {'base','curr','mean'};
models_of_interest_V = {'base','mean'};

timePeriod = 'RD';
bm_RD = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag, ...
                          'particularModel', models_of_interest_RPE);

timePeriod = 'cue';
bm_cue = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag, ...
                           'particularModels', models_of_interest_V);
                      
                       
%% fraction of neurons
timePeriod = 'RD';
switch timePeriod
    case 'RD'
        bm = bm_RD;
        mToPlot = models_of_interest_RPE;
    case 'cue'
        bm = bm_cue;
        mToPlot = models_of_interest_V;
    otherwise
        error('timePeriod not found')
end
s_VP = struct();
s_VP.n = length(bm.mask_base);

frac_VP = [];
for m = mToPlot
    m = m{:};
    s_VP.(m) = sum(bm.(['mask_' m]));
    
    frac_VP = [frac_VP s_VP.(m)/s_VP.n];
end

h_frac = figure;
bar([frac_VP; zeros(size(frac_VP))],'stacked')
legend(mToPlot,'interpreter','none')
set(gca,'tickdir','out','xtick',1,'xticklabel',{'VP'})
ylabel('Proportion of neurons')
title(['Activity at ' timePeriod])

saveFigureIteration_ottBari2019(h_frac, saveLoc, ['nFrac_' timePeriod])

%% parameter estimates (asymm model)
timePeriod = 'cue';
model = 'base';

switch model
    case 'base_asymm'
        alphaPPE = getParameters_ott(os,timePeriod,model,'alphaPPE');
        alphaNPE = getParameters_ott(os,timePeriod,model,'alphaNPE');
        asymmInd = alphaPPE ./ (alphaPPE + alphaNPE);
    case 'base'
        alpha = getParameters_ott(os,timePeriod,model,'alpha');
end
slope = getParameters_ott(os,timePeriod,model,'slope');
int = getParameters_ott(os,timePeriod,model,'intercept');
rho = getParameters_ott(os,timePeriod,model,'rho');

model_name = ['mask_' model];
n_VP = bm_RD.(model_name);

h_parameters = figure;
for i = 1:6
    h_sp(i) = subplot(2,3,i); hold on
end

alpha_bins = linspace(0,1,20);
switch model
    case 'base_asymm'
        i = 1;
        % alphaPPE
        subplot(h_sp(i))
        histogram(alphaPPE(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        set(h_sp(1),'tickdir','out')
        xlabel('$\alpha_{PPE}$','Interpreter','latex')
        ylabel('Probability')
        
        i = 2;
        % alphaNPE
        subplot(h_sp(i))
        histogram(alphaNPE(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        set(h_sp(2),'tickdir','out')
        xlabel('$\alpha_{NPE}$','Interpreter','latex')
        ylabel('Probability')
        
        % asymmetry
        i = 3;
        subplot(h_sp(i))
        histogram(asymmInd(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        set(h_sp(2),'tickdir','out')
        xlabel('Asymmetry $(\frac{\alpha_{PPE}}{\alpha_{PPE} + \alpha_{NPE}})$','Interpreter','latex')
        ylabel('Probability')
    case 'base'
        i = 1;
        % alpha
        subplot(h_sp(i))
        histogram(alpha(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        set(h_sp(1),'tickdir','out')
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('Probability')
        i = 3;
end

i = i + 1;
% slope
subplot(h_sp(i))
slope_bins = floor(min(slope(n_VP))) - 1:0.5:ceil(max(slope(n_VP))) + 1;
histogram(slope(n_VP), slope_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
set(gca,'tickdir','out')
xlabel('Slope','Interpreter','latex')
ylabel('Probability')

i = i + 1;
% intercept
subplot(h_sp(i))
int_bins = floor(min(int(n_VP))) - 1:0.5:ceil(max(int(n_VP))) + 1;
histogram(int(n_VP), int_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
set(gca,'tickdir','out')
xlabel('Intercept','Interpreter','latex')
ylabel('Probability')

i = i + 1;
% rho
subplot(h_sp(i))
rho_bins = alpha_bins;
histogram(rho(n_VP),rho_bins,'Normalization','probability','FaceColor',VP_color,'EdgeColor','none')
set(gca,'tickdir','out')
xlabel('$\rho$','Interpreter','latex')
ylabel('Probability')

suptitle(strrep(['Model: ' model ' at ' timePeriod],'_',' '))

saveFigureIteration_ottBari2019(h_parameters, saveLoc, ['params_' timePeriod '_' model])
