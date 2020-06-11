clear; clc

load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'threeOutcomes_MLEfits.mat'));
myColors = importColors_bb;
VP_color = myColors.bluishGreen;

% get relevant behavior models
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

%%
timePeriod = 'RD'; % RD, cue
normalization = 'zscore'; % none, zscore, minmax

mod_type = ['mod_' timePeriod];
switch timePeriod
    case 'RD'
        latent_var = 'RPEs';
        bm = bm_RD;
    case 'cue'
        latent_var = 'V';
        bm = bm_cue;
    otherwise
        error('timePeriod not found')
end

nSim = 501; % must be odd

% VP first
all_latent_VP = [];
norm_fr_real_VP = [];
norm_fr_sim_VP = [];

mean_real_VP = [];
mean_pred_VP = [];
corr_spike_count_VP = [];
var_real_VP = [];
var_pred_VP = [];
trialComparison_pred_VP = [];
for n = find(bm.mask_base)
    sign_flip = sign(os(n).(mod_type).base.bestParams(2));
    
    real_spike_count = os(n).(['spikeCount_' timePeriod])';
    
    % generate nSim predicted spike counts
    pred_spike_count = [];
    for i = 1:nSim
        pred_spike_count(i,:) = poissrnd(os(n).(mod_type).base.mean_predictedSpikes)';
    end
    tmp_corr = corr(real_spike_count', pred_spike_count');
    
    % get median correlation
    corr_spike_count_VP = [corr_spike_count_VP median(tmp_corr)];
    
    % save that median neuron
    median_ind = find(tmp_corr == median(tmp_corr), 1);
        
    % get latent variables for plotting
    tmp_latent = os(n).(mod_type).base.(latent_var)';
    % normalize
    switch normalization
        case 'none'
        case 'zscore'
            tmp_latent = normalize(tmp_latent);
        case 'minmax'
            norm_const = 1/max(abs(tmp_latent));
            tmp_latent = norm_const*tmp_latent;
    end
    all_latent_VP = [all_latent_VP tmp_latent];      
    
    % normalize real and predicted spike counts for tuning curves
    norm_fr_real_VP = [norm_fr_real_VP sign_flip*normalize(real_spike_count)];
    norm_fr_sim_VP = [norm_fr_sim_VP sign_flip*normalize(pred_spike_count(median_ind, :))];
    
    % mean spike counts
    mean_real_VP = [mean_real_VP mean(real_spike_count)];
    mean_pred_VP = [mean_pred_VP mean(pred_spike_count(median_ind, :))];
    
    % STD of spike counts; must use simulated spike counts here
    var_real_VP = [var_real_VP var(real_spike_count)];
    var_pred_VP = [var_pred_VP var(pred_spike_count(median_ind, :))];
    
    % save the median simulated spike count
    trialComparison_pred_VP = [trialComparison_pred_VP {pred_spike_count(median_ind, :)}];
end


nBins = 11;
latent_bins = prctile(all_latent_VP, linspace(0, 100, nBins));
spike_bins_real_VP = arrayfun(@(i, j) norm_fr_real_VP(all_latent_VP >= i & all_latent_VP < j), latent_bins(1:end -1), ...
    latent_bins(2:end), 'UniformOutput', false);
spike_bins_sim_VP = arrayfun(@(i, j) norm_fr_sim_VP(all_latent_VP >= i & all_latent_VP < j), latent_bins(1:end -1), ...
    latent_bins(2:end), 'UniformOutput', false);

% figure
x_latent_bins = latent_bins(1:end - 1) + diff(latent_bins)/2;
scatterSize = 15;

h_figure = figure; 
h_lat_VP = subplot(221); hold on
h_mean = subplot(222); hold on
h_corr = subplot(223); hold on
h_var = subplot(224); hold on

subplot(h_lat_VP)
t_lat_VP(1) = plotFilled(x_latent_bins, spike_bins_real_VP, VP_color, h_lat_VP);
t_lat_VP(2) = plotFilled(x_latent_bins, spike_bins_sim_VP, myColors.blue_bright, h_lat_VP);

subplot(h_mean)
scatter(mean_real_VP, mean_pred_VP, scatterSize, 'filled', 'MarkerFaceColor', VP_color)
maxVal = max([mean_real_VP mean_pred_VP]);
plot([0 maxVal],[0 maxVal],'k--')
xlabel('Real spikes (mean)')
ylabel('Predicted spikes (mean)')

subplot(h_corr)
corr_bins = linspace(-1,1,40);
histogram(corr_spike_count_VP, corr_bins, 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceColor', VP_color)
xlabel('Correlation')
ylabel('Probability')
xlim([-0.05 1.05])

subplot(h_var)
scatter(var_real_VP, var_pred_VP, scatterSize, 'filled', 'MarkerFaceColor', VP_color)
maxVal = max([var_real_VP var_pred_VP]);
plot([0 maxVal],[0 maxVal],'k--')
xlabel('Real spikes (variance)')
ylabel('Predicted spikes (variance)')

% clean it up
legend(t_lat_VP, {'VP','Predicted'}, 'location', 'best')
for cP = [h_lat_VP h_mean h_corr h_var]
    subplot(cP)
    set(cP,'tickdir','out')
    if cP == h_lat_VP
        xlabel(latent_var)
        ylabel([timePeriod ' spikes (z-score)'])
        if strcmp(latent_var, 'RPEs')
%             xlim([-1 1])
        else
%             xlim([0 1])
        end
    elseif cP == h_mean
        xlim([0 20]); ylim([0 20])
    elseif cP == h_corr
    elseif cP == h_var
        xlim([0 60]); ylim([0 60])
    end
end


%% plot neurons with particular cross correlations; go with median-correlation simulated spike count
prtile_cutoff = 80;

VP_neur = find(bm.mask_base);
[~, VP_neuron_ind] = min(abs(corr_spike_count_VP - prctile(corr_spike_count_VP, prtile_cutoff)));
VP_neuron_corr = corr_spike_count_VP(VP_neuron_ind);
os_VP = os(VP_neur);

os_VP_neur = os_VP(VP_neuron_ind);

h_figure = figure;
h_VP = subplot(1,1,1); hold on

subplot(h_VP)
plot(os_VP_neur.spikeCount_RD, 'Color', VP_color, 'linewidth', 2)
plot(trialComparison_pred_VP{VP_neuron_ind}, 'Color', myColors.blue_bright, 'linewidth', 2)
title(sprintf('VP\nnind: %i, corr = %0.2f', VP_neuron_ind, VP_neuron_corr))


for cp = h_VP
    subplot(cp)
    xlabel('Trials')
    ylabel('Spike count')
    set(cp,'tickdir','out')
end
