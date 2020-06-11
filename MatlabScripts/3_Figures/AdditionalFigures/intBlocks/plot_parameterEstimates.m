clear; clc
load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'intBlocks_MLEfits.mat'));
int_task = [os.Blocks] == 0;
os = os(int_task);
VP_mask = contains({os.Region}, 'VP');
myColors = importColors_bb;
VP_color = myColors.bluishGreen;
NAc_color = myColors.vermillion;

%% get relevant behavior models
modelCriterion = 'AIC';
plotFlag = false;

models_of_interest_RPE = {'base','curr','mean'};
models_of_interest_V = {'base','mean'};

timePeriod = 'RD';
bm_RD = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,...
                          'plotModels_Flag',plotFlag,...
                          'particularModel',models_of_interest_RPE);

timePeriod = 'cue';
bm_cue = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,...
                           'plotModels_Flag',plotFlag,...
                           'particularModel',models_of_interest_V);
                       
                       
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
s_NAc = struct();
s_VP.n = sum(VP_mask);
s_NAc.n = sum(~VP_mask);

frac_VP = [];
frac_NAc = [];
for m = mToPlot
    m = m{:};
    s_VP.(m) = sum(bm.(['mask_' m]) & VP_mask);
    s_NAc.(m) = sum(bm.(['mask_' m]) & ~VP_mask);
    
    frac_VP = [frac_VP s_VP.(m)/s_VP.n];
    frac_NAc = [frac_NAc s_NAc.(m)/s_NAc.n];
end

h_frac = figure;
bar([frac_VP; frac_NAc],'stacked')
legend(mToPlot,'interpreter','none')
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'VP','NAc'})
ylabel('Proportion of neurons')
title(['Activity at ' timePeriod])


%% parameter estimates (asymm model)
timePeriod = 'RD';
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

model_name = ['mask_' model];
n_VP = VP_mask & bm_RD.(model_name);
n_NAc = ~VP_mask & bm_RD.(model_name);

h_parameters = figure;
for i = 1:5
    h_sp(i) = subplot(2,3,i); hold on
end

alpha_bins = linspace(0,1,20);
switch model
    case 'base_asymm'
        i = 1;
        % alphaPPE
        subplot(h_sp(i))
        histogram(alphaPPE(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        histogram(alphaPPE(n_NAc), alpha_bins, 'Normalization', 'probability', 'FaceColor', NAc_color, 'EdgeColor', 'none')
        set(h_sp(1),'tickdir','out')
        xlabel('$\alpha_{PPE}$','Interpreter','latex')
        ylabel('Probability')
        legend('VP','NAc')
        
        i = 2;
        % alphaNPE
        subplot(h_sp(i))
        histogram(alphaNPE(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        histogram(alphaNPE(n_NAc), alpha_bins, 'Normalization', 'probability', 'FaceColor', NAc_color, 'EdgeColor', 'none')
        set(h_sp(2),'tickdir','out')
        xlabel('$\alpha_{NPE}$','Interpreter','latex')
        ylabel('Probability')
        legend('VP','NAc')
        
        % asymmetry
        i = 3;
        subplot(h_sp(i))
        histogram(asymmInd(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        histogram(asymmInd(n_NAc), alpha_bins, 'Normalization', 'probability', 'FaceColor', NAc_color, 'EdgeColor', 'none')
        set(h_sp(2),'tickdir','out')
        xlabel('Asymmetry $(\frac{\alpha_{PPE}}{\alpha_{PPE} + \alpha_{NPE}})$','Interpreter','latex')
        ylabel('Probability')
        legend('VP','NAc')
    case 'base'
        i = 1;
        % alpha
        subplot(h_sp(i))
        histogram(alpha(n_VP), alpha_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
        histogram(alpha(n_NAc), alpha_bins, 'Normalization', 'probability', 'FaceColor', NAc_color, 'EdgeColor', 'none')
        set(h_sp(1),'tickdir','out')
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('Probability')
        legend('VP','NAc')
        i = 3;
end

i = i + 1;
% slope
subplot(h_sp(i))
slope_bins = floor(min([slope(n_VP) slope(n_NAc)])) - 1:0.5:ceil(max([slope(n_VP) slope(n_NAc)])) + 1;
histogram(slope(n_VP), slope_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
histogram(slope(n_NAc), slope_bins, 'Normalization', 'probability', 'FaceColor', NAc_color, 'EdgeColor', 'none')
set(gca,'tickdir','out')
xlabel('Slope','Interpreter','latex')
ylabel('Probability')
legend('VP','NAc')

i = i + 1;
% intercept
subplot(h_sp(i))
int_bins = floor(min([int(n_VP) int(n_NAc)])) - 1:0.5:ceil(max([int(n_VP) int(n_NAc)])) + 1;
histogram(int(n_VP), int_bins, 'Normalization', 'probability', 'FaceColor', VP_color, 'EdgeColor', 'none')
histogram(int(n_NAc), int_bins, 'Normalization', 'probability', 'FaceColor', NAc_color, 'EdgeColor', 'none')
set(gca,'tickdir','out')
xlabel('Intercept','Interpreter','latex')
ylabel('Probability')
legend('VP','NAc')

sgtitle(strrep(['Model: ' model ' at ' timePeriod],'_',' '))
