clear; clc
load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'cue_MLEfits.mat'));

myColors = importColors_bb;
VP_color = myColors.bluishGreen;

% get relevant behavior models
modelCriterion = 'AIC';
plotFlag = false;
m_RPE = {'base','base_cue','curr','curr_cue','mean','mean_cue'};
m_V = {'base','base_cue','mean','mean_cue'};

timePeriod = 'RD';
bm_RD = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag,...
                          'particularModel',m_RPE);

timePeriod = 'cue';
bm_cue = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag,...
                           'particularModel',m_V);
                       
%% enrichment among cue neurons
predCue_RD = bm_RD.mask_base_cue | bm_RD.mask_curr_cue | bm_RD.mask_mean_cue;
predCue_cue = bm_cue.mask_base_cue | bm_cue.mask_mean_cue;

n_cueRD_cueCue = sum(predCue_cue(predCue_RD));
n_notcueRD_cueCue = sum(predCue_cue(~predCue_RD));
n_cueRD = sum(predCue_RD);
n_notcueRD = sum(~predCue_RD);

[~,p] = prop_test([n_cueRD_cueCue n_notcueRD_cueCue], [n_cueRD n_notcueRD]);
% no enrichment

fprintf('\n------\n')
fprintf('%i of %i RD predictive-cue neurons (%0.2f%%) have predictive-cue response at cue\n', n_cueRD_cueCue, n_cueRD, n_cueRD_cueCue/n_cueRD*100)
fprintf('%i of %i RD non-predictive-cue neurons (%0.2f%%) have predictive-cue response at cue\n', ...
         n_notcueRD_cueCue, n_notcueRD, n_notcueRD_cueCue/n_notcueRD*100)
fprintf('pValue: %0.2f\n',p)
fprintf('------\n')