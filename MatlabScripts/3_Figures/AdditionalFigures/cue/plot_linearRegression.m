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

%%
timePeriod = 'RD';
nanPad = 50;

mod_type = ['mod_' timePeriod];
switch timePeriod
    case 'RD'
        bm = bm_RD;
        mod2 = 'curr';
    case 'cue'
        bm = bm_cue;
        mod2 = 'prev';
    otherwise
        error('timePeriod not found')
end

% VP first
rwd_base_VP = [];
fr_base_VP = [];
fr_base_pred_VP = [];

rwd_mod2_VP = [];
fr_mod2_VP = [];
fr_mod2_pred_VP = [];

rwd_mean_VP = [];
fr_mean_VP = [];
fr_mean_pred_VP = [];

for n = 1:length(os)
    spike_count = os(n).(['spikeCount_' timePeriod]);

    rewards = os(n).rewards(os(n).timeLocked);
    if bm.mask_base(n) == 1 % base neuron
        sign_flip = sign(os(n).(mod_type).base.bestParams(2));
        pred_spike_count = poissrnd(os(n).(mod_type).base.mean_predictedSpikes);
        
        rwd_base_VP = [rwd_base_VP; NaN(nanPad, 1); rewards];
        fr_base_VP = [fr_base_VP; NaN(nanPad, 1); normalize(spike_count)];
        fr_base_pred_VP = [fr_base_pred_VP; NaN(nanPad, 1); normalize(pred_spike_count)];
    elseif bm.(['mask_' mod2])(n) == 1 || bm.(['mask_' mod2 '_cue'])(n) == 1 % curr/prev neuron
        pred_spike_count = poissrnd(os(n).(mod_type).(mod2).mean_predictedSpikes);

        rwd_mod2_VP = [rwd_mod2_VP; NaN(nanPad, 1); rewards];
        fr_mod2_VP = [fr_mod2_VP; NaN(nanPad, 1); normalize(spike_count)];
        fr_mod2_pred_VP = [fr_mod2_pred_VP; NaN(nanPad, 1); normalize(pred_spike_count)];
    elseif bm.mask_mean(n) == 1 || bm.mask_mean_cue(n) == 1 % mean neuron
        pred_spike_count = poissrnd(mean(spike_count), length(spike_count), 1);

        rwd_mean_VP = [rwd_mean_VP; NaN(nanPad, 1); rewards];
        fr_mean_VP = [fr_mean_VP; NaN(nanPad, 1); normalize(spike_count)];
        fr_mean_pred_VP = [fr_mean_pred_VP; NaN(nanPad, 1); normalize(pred_spike_count)];
    end
end

rwdHx_base_VP = [rwd_base_VP generateHistoryMatrix(rwd_base_VP, 10)];
rwdHx_mod2_VP = [rwd_mod2_VP generateHistoryMatrix(rwd_mod2_VP, 10)];
rwdHx_mean_VP = [rwd_mean_VP generateHistoryMatrix(rwd_mean_VP, 10)];
base_mod_VP = fitlm(rwdHx_base_VP, fr_base_VP);
mod2_mod_VP = fitlm(rwdHx_mod2_VP, fr_mod2_VP);
mean_mod_VP = fitlm(rwdHx_mean_VP, fr_mean_VP);
base_pred_mod_VP = fitlm(rwdHx_base_VP, fr_base_pred_VP);
mod2_pred_mod_VP = fitlm(rwdHx_mod2_VP, fr_mod2_pred_VP);
mean_pred_mod_VP = fitlm(rwdHx_mean_VP, fr_mean_pred_VP);


% figure
h_bg = figure; 
h_VP = subplot(121); hold on
h_VP_pred = subplot(122); hold on

subplot(h_VP)
title('VP')
t_VP(1) = plotRegressionWithCI(base_mod_VP, 2:12, h_VP, 'XOffset', -1, 'Color', VP_color);
t_VP(2) = plotRegressionWithCI(mod2_mod_VP, 2:12, h_VP, 'XOffset', -1, 'Color', myColors.darkGray);
t_VP(3) = plotRegressionWithCI(mean_mod_VP, 2:12, h_VP, 'XOffset', -1, 'Color', myColors.lightGray);

subplot(h_VP_pred)
title('VP - simulated neurons')
t_VP_sim(1) = plotRegressionWithCI(base_pred_mod_VP, 2:12, h_VP_pred, 'XOffset', -1, 'Color', VP_color);
t_VP_sim(2) = plotRegressionWithCI(mod2_pred_mod_VP, 2:12, h_VP_pred, 'XOffset', -1, 'Color', myColors.darkGray);
t_VP_sim(3) = plotRegressionWithCI(mean_pred_mod_VP, 2:12, h_VP_pred, 'XOffset', -1, 'Color', myColors.lightGray);

for cP = [h_VP h_VP_pred]
    subplot(cP)
    set(cP,'tickdir','out')
    plot([-0.5 10.5],[0 0],'k--')
    xlim([-0.5 10.5])
    ylim([-0.6 2.0])
    xlabel('Reward n trials back')
    ylabel('Coefficient ($\pm 95\%$ CI)', 'Interpreter', 'latex')
end
legend(t_VP, {'Base', mod2, 'Mean'})
legend(t_VP_sim, {'Base', mod2, 'Mean'})