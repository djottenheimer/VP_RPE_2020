clear; clc
load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'cue_MLEfits.mat'));

myColors = importColors_bb;
VP_color = myColors.bluishGreen;

% get relevant behavior models
modelCriterion = 'AIC';
plotFlag = true;
m_RPE = {'base','base_cue','curr','curr_cue','mean','mean_cue'};
m_V = {'base','base_cue','mean','mean_cue'};

timePeriod = 'RD';
bm_RD = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag,...
                          'particularModel',m_RPE);

timePeriod = 'cue';
bm_cue = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag,...
                           'particularModel',m_V);
                       
%% RD proportion
nTot = length(bm_RD.mask_base);
n_noCue = [sum(bm_RD.mask_base) sum(bm_RD.mask_curr) sum(bm_RD.mask_mean)];
n_cue = [sum(bm_RD.mask_base_cue) sum(bm_RD.mask_curr_cue) sum(bm_RD.mask_mean_cue)];

n_noCue = n_noCue./nTot;
n_cue = n_cue./nTot;


h_RD = figure;
bar([n_noCue; n_cue],'stacked')
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'No cue effect','Cue effect'})
legend('RPE','Curr','Mean','interpreter','none')
ylabel('Proportion')
title('RD')
ylim([0 1])
                            
%% cue proportion
nTot = length(bm_cue.mask_base);
n_noCue = [sum(bm_cue.mask_base) sum(bm_cue.mask_mean)];
n_cue = [sum(bm_cue.mask_base_cue) sum(bm_cue.mask_mean_cue)];

n_noCue = n_noCue./nTot;
n_cue = n_cue./nTot;


h_cue = figure;
bar([n_noCue; n_cue],'stacked')
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'No cue effect','Cue effect'})
legend('Value','Mean','interpreter','none')
ylabel('Proportion')
title('Cue')
ylim([0 1])