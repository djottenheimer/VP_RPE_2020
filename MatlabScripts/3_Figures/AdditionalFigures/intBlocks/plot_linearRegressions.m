load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'intBlocks_MLEfits.mat'));
int_task = [os.Blocks] == 0;
os = os(int_task);
VP_mask = contains({os.Region}, 'VP');
myColors = importColors_bb;
VP_color = myColors.bluishGreen;
NAc_color = myColors.vermillion;

% get relevant behavior models
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

%%
timePeriod = 'RD'; % RD or cue
nanPad = 50;

mod_type = ['mod_' timePeriod];
switch timePeriod
    case 'RD'
        latent_var = 'RPEs';
        bm = bm_RD;
        mod2 = 'curr';
    case 'cue'
        latent_var = 'V';
        bm = bm_cue;
    otherwise
        error('timePeriod not found')
end

% VP first
rwd_base_VP = [];
fr_base_VP = [];
fr_base_pred_VP = [];

if strcmp(timePeriod,'RD')
    rwd_mod2_VP = [];
    fr_mod2_VP = [];
    fr_mod2_pred_VP = [];
end

rwd_mean_VP = [];
fr_mean_VP = [];
fr_mean_pred_VP = [];

for n = 1:length(os)
    if VP_mask(n) == 1 % VP neuron
        spike_count = os(n).(['spikeCount_' timePeriod]);
        
        rewards = os(n).rewards(os(n).timeLocked);
        if bm.mask_base(n) == 1 % base neuron
            sign_flip = sign(os(n).(mod_type).base.bestParams(2));
            pred_spike_count = poissrnd(os(n).(mod_type).base.mean_predictedSpikes);

            rwd_base_VP = [rwd_base_VP; NaN(nanPad, 1); rewards];
            fr_base_VP = [fr_base_VP; NaN(nanPad, 1); sign_flip*normalize(spike_count)];
            fr_base_pred_VP = [fr_base_pred_VP; NaN(nanPad, 1); sign_flip*normalize(pred_spike_count)];
        elseif strcmp(timePeriod,'RD') && bm.(['mask_' mod2])(n) == 1 % curr/prev neuron
            sign_flip = sign(os(n).(mod_type).(mod2).bestParams(1));
            pred_spike_count = poissrnd(os(n).(mod_type).(mod2).mean_predictedSpikes);
            
            rwd_mod2_VP = [rwd_mod2_VP; NaN(nanPad, 1); rewards];
            fr_mod2_VP = [fr_mod2_VP; NaN(nanPad, 1); sign_flip*normalize(spike_count)];
            fr_mod2_pred_VP = [fr_mod2_pred_VP; NaN(nanPad, 1); sign_flip*normalize(pred_spike_count)];
        elseif bm.mask_mean(n) == 1 % mean neuron
            pred_spike_count = poissrnd(mean(spike_count), length(spike_count), 1);
            
            rwd_mean_VP = [rwd_mean_VP; NaN(nanPad, 1); rewards];
            fr_mean_VP = [fr_mean_VP; NaN(nanPad, 1); normalize(spike_count)];
            fr_mean_pred_VP = [fr_mean_pred_VP; NaN(nanPad, 1); normalize(pred_spike_count)];
        end
    end
end

rwdHx_base_VP = [rwd_base_VP generateHistoryMatrix(rwd_base_VP, 10)];
rwdHx_mean_VP = [rwd_mean_VP generateHistoryMatrix(rwd_mean_VP, 10)];
base_mod_VP = fitlm(rwdHx_base_VP, fr_base_VP);
mean_mod_VP = fitlm(rwdHx_mean_VP, fr_mean_VP);
base_pred_mod_VP = fitlm(rwdHx_base_VP, fr_base_pred_VP);
mean_pred_mod_VP = fitlm(rwdHx_mean_VP, fr_mean_pred_VP);
if strcmp(timePeriod,'RD')
    rwdHx_mod2_VP = [rwd_mod2_VP generateHistoryMatrix(rwd_mod2_VP, 10)];
    mod2_mod_VP = fitlm(rwdHx_mod2_VP, fr_mod2_VP);
    mod2_pred_mod_VP = fitlm(rwdHx_mod2_VP, fr_mod2_pred_VP);
end

% NAc second
rwd_base_NAc = [];
fr_base_NAc = [];
fr_base_pred_NAc = [];

if strcmp(timePeriod,'RD')
    rwd_mod2_NAc = [];
    fr_mod2_NAc = [];
    fr_mod2_pred_NAc = [];
end

rwd_mean_NAc = [];
fr_mean_NAc = [];
fr_mean_pred_NAc = [];

for n = 1:length(os)
    if VP_mask(n) == 0 % NAc neuron
        spike_count = os(n).(['spikeCount_' timePeriod]);
        rewards = os(n).rewards(os(n).timeLocked);
        if bm.mask_base(n) == 1 % base neuron
            sign_flip = sign(os(n).(mod_type).base.bestParams(2));
            pred_spike_count = poissrnd(os(n).(mod_type).base.mean_predictedSpikes);

            rwd_base_NAc = [rwd_base_NAc; NaN(nanPad, 1); rewards];
            fr_base_NAc = [fr_base_NAc; NaN(nanPad, 1); sign_flip*normalize(spike_count)];
            fr_base_pred_NAc = [fr_base_pred_NAc; NaN(nanPad, 1); sign_flip*normalize(pred_spike_count)];
        elseif strcmp(timePeriod,'RD') && bm.(['mask_' mod2])(n) == 1 % curr/prev neuron
            sign_flip = sign(os(n).(mod_type).(mod2).bestParams(1));
            pred_spike_count = poissrnd(os(n).(mod_type).(mod2).mean_predictedSpikes);
            
            rwd_mod2_NAc = [rwd_mod2_NAc; NaN(nanPad, 1); rewards];
            fr_mod2_NAc = [fr_mod2_NAc; NaN(nanPad, 1); sign_flip*normalize(spike_count)];
            fr_mod2_pred_NAc = [fr_mod2_pred_NAc; NaN(nanPad, 1); sign_flip*normalize(pred_spike_count)];
        elseif bm.mask_mean(n) == 1 % mean neuron
            pred_spike_count = poissrnd(mean(spike_count), length(spike_count), 1);
            
            rwd_mean_NAc = [rwd_mean_NAc; NaN(nanPad, 1); rewards];
            fr_mean_NAc = [fr_mean_NAc; NaN(nanPad, 1); normalize(spike_count)];
            fr_mean_pred_NAc = [fr_mean_pred_NAc; NaN(nanPad, 1); normalize(pred_spike_count)];
        end
    end
end

rwdHx_base_NAc = [rwd_base_NAc generateHistoryMatrix(rwd_base_NAc, 10)];
rwdHx_mean_NAc = [rwd_mean_NAc generateHistoryMatrix(rwd_mean_NAc, 10)];
base_mod_NAc = fitlm(rwdHx_base_NAc, fr_base_NAc);
mean_mod_NAc = fitlm(rwdHx_mean_NAc, fr_mean_NAc);
base_pred_mod_NAc = fitlm(rwdHx_base_NAc, fr_base_pred_NAc);
mean_pred_mod_NAc = fitlm(rwdHx_mean_NAc, fr_mean_pred_NAc);
if strcmp(timePeriod,'RD')
    rwdHx_mod2_NAc = [rwd_mod2_NAc generateHistoryMatrix(rwd_mod2_NAc, 10)];
    mod2_mod_NAc = fitlm(rwdHx_mod2_NAc, fr_mod2_NAc);
    mod2_pred_mod_NAc = fitlm(rwdHx_mod2_NAc, fr_mod2_pred_NAc);
end


% figure
h_bg = figure; 
h_VP = subplot(221); hold on
h_VP_pred = subplot(222); hold on
h_NAc = subplot(223); hold on
h_NAc_pred = subplot(224); hold on

subplot(h_VP)
title('VP')
t_VP(1) = plotRegressionWithCI(base_mod_VP, 2:12, h_VP, 'XOffset', -1, 'Color', VP_color);
if strcmp(timePeriod,'RD')
    t_VP(2) = plotRegressionWithCI(mod2_mod_VP, 2:12, h_VP, 'XOffset', -1, 'Color', myColors.darkGray);
end
t_VP(3) = plotRegressionWithCI(mean_mod_VP, 2:12, h_VP, 'XOffset', -1, 'Color', myColors.lightGray);

subplot(h_NAc)
title('NAc')
t_NAc(1) = plotRegressionWithCI(base_mod_NAc, 2:12, h_NAc, 'XOffset', -1, 'Color', NAc_color);
if strcmp(timePeriod,'RD')
    t_NAc(2) = plotRegressionWithCI(mod2_mod_NAc, 2:12, h_NAc, 'XOffset', -1, 'Color', myColors.darkGray);
end
t_NAc(3) = plotRegressionWithCI(mean_mod_NAc, 2:12, h_NAc, 'XOffset', -1, 'Color', myColors.lightGray);

subplot(h_VP_pred)
title('VP - simulated neurons')
t_VP_sim(1) = plotRegressionWithCI(base_pred_mod_VP, 2:12, h_VP_pred, 'XOffset', -1, 'Color', VP_color);
if strcmp(timePeriod,'RD')
    t_VP_sim(2) = plotRegressionWithCI(mod2_pred_mod_VP, 2:12, h_VP_pred, 'XOffset', -1, 'Color', myColors.darkGray);
end
t_VP_sim(3) = plotRegressionWithCI(mean_pred_mod_VP, 2:12, h_VP_pred, 'XOffset', -1, 'Color', myColors.lightGray);

subplot(h_NAc_pred)
title('NAc - simulated neurons')
t_NAc_sim(1) = plotRegressionWithCI(base_pred_mod_NAc, 2:12, h_NAc_pred, 'XOffset', -1, 'Color', NAc_color);
if strcmp(timePeriod,'RD')
    t_NAc_sim(2) = plotRegressionWithCI(mod2_pred_mod_NAc, 2:12, h_NAc_pred, 'XOffset', -1, 'Color', myColors.darkGray);
end
t_NAc_sim(3) = plotRegressionWithCI(mean_pred_mod_NAc, 2:12, h_NAc_pred, 'XOffset', -1, 'Color', myColors.lightGray);

for cP = [h_VP h_NAc h_VP_pred h_NAc_pred]
    subplot(cP)
    set(cP,'tickdir','out')
    plot([-0.5 10.5],[0 0],'k--')
    xlim([-0.5 10.5])
    ylim([-0.6 1.1])
    xlabel('Reward n trials back')
    ylabel('Coefficient ($\pm 95\%$ CI)', 'Interpreter', 'latex')
end
if strcmp(timePeriod,'RD')
    legend(t_VP, {'Base', mod2, 'Mean'})
    legend(t_VP_sim, {'Base', mod2, 'Mean'})
    legend(t_NAc, {'Base', mod2, 'Mean'})
    legend(t_NAc_sim, {'Base', mod2, 'Mean'})
else
    t_VP(2) = [];
    t_VP_sim(2) = [];
    t_NAc(2) = [];
    t_NAc_sim(2) = [];
    legend(t_VP, {'Base', 'Mean'})
    legend(t_VP_sim, {'Base', 'Mean'})
    legend(t_NAc, {'Base', 'Mean'})
    legend(t_NAc_sim, {'Base', 'Mean'})
end