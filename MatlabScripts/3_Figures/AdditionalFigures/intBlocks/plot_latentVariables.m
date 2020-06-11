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

%% real vs predicted tuning curves, correlations
timePeriod = 'RD'; % RD or cue
normalization = 'zscore'; % (latent variable) none, zscore, minmax, zscore_translate (so 0 = 0)

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
for n = find(bm.mask_base & VP_mask)
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
    median_ind = find(tmp_corr == nanmedian(tmp_corr), 1);
        
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
        case 'zscore_translate'
            tmp_latent = normalize(tmp_latent) + mean(tmp_latent)/std(tmp_latent);
    end
    all_latent_VP = [all_latent_VP tmp_latent];
    
    % normalize real and predicted spike counts for tuning curves
    norm_fr_real_VP = [norm_fr_real_VP sign_flip*normalize(real_spike_count)];
    norm_fr_sim_VP = [norm_fr_sim_VP sign_flip*normalize(pred_spike_count(median_ind, :))];
    
    if numel(norm_fr_real_VP) ~= numel(norm_fr_sim_VP)
        fprintf('Returning\n')
        return
    end
    
    % mean spike counts
    mean_real_VP = [mean_real_VP mean(real_spike_count)];
    mean_pred_VP = [mean_pred_VP mean(pred_spike_count(median_ind, :))];
    
    % STD of spike counts; must use simulated spike counts here
    var_real_VP = [var_real_VP var(real_spike_count)];
    var_pred_VP = [var_pred_VP var(pred_spike_count(median_ind, :))];
    
    % save the median simulated spike count
    trialComparison_pred_VP = [trialComparison_pred_VP {pred_spike_count(median_ind, :)}];
end


% NAc second
all_latent_NAc = [];
norm_fr_real_NAc = [];
norm_fr_sim_NAc = [];

mean_real_NAc = [];
mean_pred_NAc = [];
corr_spike_count_NAc = [];
var_real_NAc = [];
var_pred_NAc = [];
trialComparison_pred_NAc = [];
for n = find(bm.mask_base & ~VP_mask)
    sign_flip = sign(os(n).(mod_type).base.bestParams(2));

    real_spike_count = os(n).(['spikeCount_' timePeriod])';
    
    pred_spike_count = [];
    for i = 1:nSim
        pred_spike_count(i,:) = poissrnd(os(n).(mod_type).base.mean_predictedSpikes)';
    end
    tmp_corr = corr(real_spike_count', pred_spike_count');
    
    corr_spike_count_NAc = [corr_spike_count_NAc median(tmp_corr)];
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
    all_latent_NAc = [all_latent_NAc tmp_latent];
    
    norm_fr_real_NAc = [norm_fr_real_NAc sign_flip*normalize(real_spike_count)];
    norm_fr_sim_NAc = [norm_fr_sim_NAc sign_flip*normalize(pred_spike_count(median_ind, :))];
    
    mean_real_NAc = [mean_real_NAc mean(real_spike_count)];
    mean_pred_NAc = [mean_pred_NAc mean(pred_spike_count(median_ind, :))];
    
    var_real_NAc = [var_real_NAc var(real_spike_count)];
    var_pred_NAc = [var_pred_NAc var(pred_spike_count(median_ind, :))];
    
    trialComparison_pred_NAc = [trialComparison_pred_NAc {pred_spike_count(median_ind, :)}];
end

nBins = 11;
latent_bins = prctile([all_latent_VP all_latent_NAc], linspace(0, 100, nBins));
spike_bins_real_VP = arrayfun(@(i, j) norm_fr_real_VP(all_latent_VP >= i & all_latent_VP < j), latent_bins(1:end -1), ...
    latent_bins(2:end), 'UniformOutput', false);
spike_bins_sim_VP = arrayfun(@(i, j) norm_fr_sim_VP(all_latent_VP >= i & all_latent_VP < j), latent_bins(1:end -1), ...
    latent_bins(2:end), 'UniformOutput', false);

spike_bins_real_NAc = arrayfun(@(i, j) norm_fr_real_NAc(all_latent_NAc >= i & all_latent_NAc < j), latent_bins(1:end -1), ...
    latent_bins(2:end), 'UniformOutput', false);
spike_bins_sim_NAc = arrayfun(@(i, j) norm_fr_sim_NAc(all_latent_NAc >= i & all_latent_NAc < j), latent_bins(1:end -1), ...
    latent_bins(2:end), 'UniformOutput', false);


x_latent_bins = latent_bins(1:end - 1) + diff(latent_bins)/2;
scatterSize = 15;

h_figure = figure; 
h_lat_VP = subplot(2,3,1); hold on
h_lat_NAc = subplot(2,3,2); hold on
h_lat_both = subplot(233); hold on
h_mean = subplot(234); hold on
h_corr = subplot(235); hold on
h_var = subplot(236); hold on

subplot(h_lat_VP)
t_lat_VP(1) = plotFilled(x_latent_bins, spike_bins_real_VP, VP_color, h_lat_VP);
t_lat_VP(2) = plotFilled(x_latent_bins, spike_bins_sim_VP, myColors.blue_bright, h_lat_VP);

subplot(h_lat_NAc)
t_lat_NAc(1) = plotFilled(x_latent_bins, spike_bins_real_NAc, NAc_color, h_lat_NAc);
t_lat_NAc(2) = plotFilled(x_latent_bins, spike_bins_sim_NAc, myColors.blue_bright, h_lat_NAc);

subplot(h_lat_both)
t_lat_both(1) = plotFilled(x_latent_bins, spike_bins_real_VP, VP_color, h_lat_VP);
t_lat_both(2) = plotFilled(x_latent_bins, spike_bins_real_NAc, NAc_color, h_lat_NAc);

subplot(h_mean)
scatter(mean_real_VP, mean_pred_VP, scatterSize, 'filled', 'MarkerFaceColor', VP_color)
scatter(mean_real_NAc, mean_pred_NAc, scatterSize, 'filled', 'MarkerFaceColor', NAc_color)
maxVal = max([mean_real_VP mean_pred_VP mean_real_NAc mean_pred_NAc]);
plot([0 maxVal],[0 maxVal],'k--')
xlabel('Real spikes (mean)')
ylabel('Predicted spikes (mean)')

subplot(h_corr)
corr_bins = linspace(-1,1,40);
histogram(corr_spike_count_VP, corr_bins, 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceColor', VP_color)
histogram(corr_spike_count_NAc, corr_bins, 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceColor', NAc_color)
xlabel('Correlation')
ylabel('Probability')
xlim([-0.05 1.05])

subplot(h_var)
scatter(var_real_VP, var_pred_VP, scatterSize, 'filled', 'MarkerFaceColor', VP_color)
scatter(var_real_NAc, var_pred_NAc, scatterSize, 'filled', 'MarkerFaceColor', NAc_color)
maxVal = max([var_real_VP var_pred_VP var_real_NAc var_pred_NAc]);
plot([0 maxVal],[0 maxVal],'k--')
xlabel('Real spikes (variance)')
ylabel('Predicted spikes (variance)')

% clean it up
legend(t_lat_VP, {'VP','Predicted'}, 'location', 'best')
legend(t_lat_NAc, {'NAc','Predicted'}, 'location', 'best')
legend(t_lat_both, {'VP','NAc'}, 'location', 'best')
for cP = [h_lat_VP h_lat_NAc h_lat_both h_mean h_corr h_var]
    subplot(cP)
    set(cP,'tickdir','out')
    if any(cP == [h_lat_VP h_lat_NAc h_lat_both])
        xlabel(latent_var)
        ylabel([timePeriod ' spikes (z-score)'])
        if strcmp(latent_var, 'RPEs')
%             xlim([-1 1])
        else
%             xlim([0 1])
        end
    elseif cP == h_mean
        legend({'VP','NAc'}, 'location', 'best')
        xlim([0 20]); ylim([0 20])
    elseif cP == h_corr
        legend({'VP','NAc'}, 'location', 'best')
    elseif cP == h_var
        legend({'VP','NAc'}, 'location', 'best')
        xlim([0 60]); ylim([0 60])
    end
end


%% plot neurons with particular cross correlations; go with median-correlation simulated spike count
prtile_cutoff = 80;

VP_neur = find(bm.mask_base & VP_mask);
NAc_neur = find(bm.mask_base & ~VP_mask);
[~, VP_neuron_ind] = min(abs(corr_spike_count_VP - prctile(corr_spike_count_VP, prtile_cutoff)));
VP_neuron_corr = corr_spike_count_VP(VP_neuron_ind);
[~, NAc_neuron_ind] = min(abs(corr_spike_count_NAc - prctile(corr_spike_count_NAc, prtile_cutoff)));
NAc_neuron_corr = corr_spike_count_NAc(NAc_neuron_ind);
os_VP = os(VP_neur);
os_NAc = os(NAc_neur);

os_VP_neur = os_VP(VP_neuron_ind);
os_NAc_neur = os_NAc(NAc_neuron_ind);

h_figure = figure;
h_VP = subplot(211); hold on
h_NAc = subplot(212); hold on

subplot(h_VP)
plot(os_VP_neur.spikeCount_RD, 'Color', VP_color, 'linewidth', 2)
plot(trialComparison_pred_VP{VP_neuron_ind}, 'Color', myColors.blue_bright, 'linewidth', 2)
title(sprintf('VP\nnind: %i, corr = %0.2f', VP_neuron_ind, VP_neuron_corr))
legend('VP','Simulated')

subplot(h_NAc)
plot(os_NAc_neur.spikeCount_RD, 'Color', NAc_color, 'linewidth', 2)
plot(trialComparison_pred_NAc{NAc_neuron_ind}, 'Color', myColors.blue_bright, 'linewidth', 2)
title(sprintf('NAc\nnind: %i, corr = %0.2f', NAc_neuron_ind, NAc_neuron_corr))
legend('NAc','Simulated')

for cp = [h_VP h_NAc]
    subplot(cp)
    xlabel('Trials')
    ylabel('Spike count')
    set(cp,'tickdir','out')
end


%% stats for RMSE difference between real and predicted RPE vs spike count plots
mse_VP = mean((norm_fr_real_VP - norm_fr_sim_VP).^2);
mse_CI_VP = bootci(1e3, @mean, (norm_fr_real_VP - norm_fr_sim_VP).^2);

mse_NAc = mean((norm_fr_real_NAc - norm_fr_sim_NAc).^2);
mse_CI_NAc = bootci(1e3, @mean, (norm_fr_real_NAc - norm_fr_sim_NAc).^2);

fprintf('\n---\n')
fprintf('VP: %0.2f [%0.2f - %0.2f]\n', mse_VP, mse_CI_VP(1), mse_CI_VP(2))
fprintf('NAc: %0.2f [%0.2f - %0.2f]\n', mse_NAc, mse_CI_NAc(1), mse_CI_NAc(2))
fprintf('---\n\n')