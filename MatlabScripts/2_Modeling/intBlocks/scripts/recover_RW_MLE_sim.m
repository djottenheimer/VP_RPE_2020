clear os_RPE_temp os_V_temp os_RPE os_V

nStart = 10;

nTrials = 55; % median number of trials
nNeurons = 200;
Vinit = 0.5;

m_RPE = {'base_asymm','base','curr','mean'};
m_V = {'base_asymm','base','mean'};


for ind = 1:nNeurons % for all neurons
    fprintf('n %i of %i\n', ind, nNeurons)
    
    V = NaN(nTrials + 1, 1);
    RPE = NaN(nTrials, 1);
    Vasymm = NaN(nTrials + 1, 1);
    RPEasymm = NaN(nTrials, 1);
    rwd = binornd(1, 0.5, nTrials, 1);
    rwd_prev = [0; rwd(1:end - 1)];

    alpha = rand; % 0 to 1
    alphaPPE = rand; % 0 to 1
    alphaNPE = rand; % 0 to 1
    slope = rand*3 + 1; % 1 to 4
    int = rand*10 - 5; % -5 to 5
    
    V(1) = Vinit;
    Vasymm(1) = Vinit;
    for t = 1:nTrials
        RPE(t) = rwd(t) - V(t);
        V(t + 1) = V(t) + alpha*RPE(t);
        
        RPEasymm(t) = rwd(t) - Vasymm(t);
        if RPEasymm(t) >= 0
            Vasymm(t + 1) = Vasymm(t) + alphaPPE*RPEasymm(t);
        else
            Vasymm(t + 1) = Vasymm(t) + alphaNPE*RPEasymm(t);
        end
    end
    
    % simulate RD neurons
    sp_RPEasymm = poissrnd(exp(slope*RPEasymm + int));
    sp_RPE = poissrnd(exp(slope*RPE + int));
    sp_curr = poissrnd(exp(slope*rwd + int));
    
    % simulate cue neurons
    sp_Vasymm = poissrnd(exp(slope*Vasymm(1:nTrials) + int));
    sp_V = poissrnd(exp(slope*V(1:nTrials) + int));
    sp_prev = poissrnd(exp(slope*rwd_prev + int));
    
    % simulate mean neurons (both for RD and cue)
    sp_mean = poissrnd(exp(int), nTrials, 1);
    
    
    % RD neurons
    os_RPE_temp(ind).params.alpha = alpha;
    os_RPE_temp(ind).params.alphaPPE = alphaPPE;
    os_RPE_temp(ind).params.alphaNPE = alphaNPE;
    os_RPE_temp(ind).params.slope = slope;
    os_RPE_temp(ind).params.int = int;   
    
    os_RPE_temp(ind).rewards = rwd;
    os_RPE_temp(ind).timeLocked = true(size(rwd));
    os_RPE_temp(ind).spikeCount_RPEasymm = sp_RPEasymm;
    os_RPE_temp(ind).spikeCount_RPE = sp_RPE;
    os_RPE_temp(ind).spikeCount_curr = sp_curr;
    os_RPE_temp(ind).spikeCount_mean = sp_mean;
   
    
    % spikeCount is a temporary field
    % fit RPE
    os_RPE_temp(ind).spikeCount = os_RPE_temp(ind).spikeCount_RPEasymm;
    ms = helper_RW_RPE(os_RPE_temp(ind), 'StartingPoints', nStart, 'particularModel', m_RPE);
    os_RPE_temp(ind).mod_RPEasymm = ms;
    
    os_RPE_temp(ind).spikeCount = os_RPE_temp(ind).spikeCount_RPE;
    ms = helper_RW_RPE(os_RPE_temp(ind), 'StartingPoints', nStart, 'particularModel', m_RPE);
    os_RPE_temp(ind).mod_RPE = ms;
    
    os_RPE_temp(ind).spikeCount = os_RPE_temp(ind).spikeCount_curr;
    ms = helper_RW_RPE(os_RPE_temp(ind), 'StartingPoints', nStart, 'particularModel', m_RPE);
    os_RPE_temp(ind).mod_curr = ms;
    
    os_RPE_temp(ind).spikeCount = os_RPE_temp(ind).spikeCount_mean;
    ms = helper_RW_RPE(os_RPE_temp(ind), 'StartingPoints', nStart, 'particularModel', m_RPE);
    os_RPE_temp(ind).mod_mean = ms;
    
    % cue neurons
    os_V_temp(ind).params.alpha = alpha;
    os_V_temp(ind).params.alphaPPE = alphaPPE;
    os_V_temp(ind).params.alphaNPE = alphaNPE;
    os_V_temp(ind).params.slope = slope;
    os_V_temp(ind).params.int = int;   
    
    os_V_temp(ind).rewards = rwd;
    os_V_temp(ind).timeLocked = true(size(rwd));
    os_V_temp(ind).spikeCount_Vasymm = sp_Vasymm;
    os_V_temp(ind).spikeCount_V = sp_V;
    os_V_temp(ind).spikeCount_prev = sp_prev;
    os_V_temp(ind).spikeCount_mean = sp_mean;
    % spikeCount is a temporary field
    % fit V
    os_V_temp(ind).spikeCount = os_V_temp(ind).spikeCount_Vasymm;
    ms = helper_RW_V(os_V_temp(ind), 'StartingPoints', nStart, 'particularModel', m_V);
    os_V_temp(ind).mod_Vasymm = ms;
    
    os_V_temp(ind).spikeCount = os_V_temp(ind).spikeCount_V;
    ms = helper_RW_V(os_V_temp(ind), 'StartingPoints', nStart, 'particularModel', m_V);
    os_V_temp(ind).mod_V = ms;
    
    os_V_temp(ind).spikeCount = os_V_temp(ind).spikeCount_prev;
    ms = helper_RW_V(os_V_temp(ind), 'StartingPoints', nStart, 'particularModel', m_V);
    os_V_temp(ind).mod_prev = ms;
    
    os_V_temp(ind).spikeCount = os_V_temp(ind).spikeCount_mean;
    ms = helper_RW_V(os_V_temp(ind), 'StartingPoints', nStart, 'particularModel', m_V);
    os_V_temp(ind).mod_mean = ms;
    
    % remove spikeCount to avoid future confusion
    os_RPE(ind) = rmfield(os_RPE_temp(ind), 'spikeCount');
    os_V(ind) = rmfield(os_V_temp(ind), 'spikeCount');
end

fprintf('Finished\n')
save(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', ...
              'intBlocks_MLEfits_simulated_offSim.mat'), ...
    'os_RPE', 'os_V');

%%
load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', ...
              'intBlocks_MLEfits_simulated_offSim.mat'))
scoreToUse = 'AIC';
plotModels_Flag = true;
m_RPE = {'base','curr','mean'};
m_V = {'base','mean'};

bm_RPE_aRPE = select_RPEmods(os_RPE, 'RPE', 'scoreToUse', scoreToUse, ...
                                'plotModels_Flag', plotModels_Flag, 'particularModels', m_RPE);
bm_RPE_acurr = select_RPEmods(os_RPE, 'curr', 'scoreToUse', scoreToUse, ...
                                'plotModels_Flag', plotModels_Flag, 'particularModels', m_RPE);
bm_RPE_aamean = select_RPEmods(os_RPE, 'mean', 'scoreToUse', scoreToUse, ...
                                'plotModels_Flag', plotModels_Flag, 'particularModels', m_RPE);
                    
bm_V_aV = select_RPEmods(os_V, 'V', 'particularModels', m_V, 'scoreToUse', scoreToUse, ...
                             'plotModels_Flag', plotModels_Flag);
bm_V_amean = select_RPEmods(os_V, 'mean', 'particularModels', m_V, 'scoreToUse', scoreToUse, ...
                             'plotModels_Flag', plotModels_Flag);
                         
myColors = importColors_bb;

%%
% RPE first
model = 'V';
switch model
    case 'RPE'
        m1 = bm_RPE_aRPE;
        m2 = bm_RPE_acurr;
        m3 = bm_RPE_aamean;
        mod2 = 'curr';
        nNeurons = length(os_RPE);
    case 'V'
        m1 = bm_V_aV;
        m2 = bm_V_amean;
        nNeurons = length(os_V);
    otherwise
        error('model not found')
end

switch model
    case 'RPE'
        aBase_rBase = sum(m1.mask_base)/nNeurons;
        aBase_rMod2 = sum(m1.(['mask_' mod2]))/nNeurons;
        aBase_rMean = sum(m1.mask_mean)/nNeurons;

        aMod2_rBase = sum(m2.mask_base)/nNeurons;
        aMod2_rMod2 = sum(m2.(['mask_' mod2]))/nNeurons;
        aMod2_rMean = sum(m2.mask_mean)/nNeurons;

        aMean_rBase = sum(m3.mask_base)/nNeurons;
        aMean_rMod2 = sum(m3.(['mask_' mod2]))/nNeurons;
        aMean_rMean = sum(m3.mask_mean)/nNeurons;

        mat_for_hmap = [aBase_rBase      aBase_rMod2      aBase_rMean;
                        aMod2_rBase      aMod2_rMod2      aMod2_rMean;
                        aMean_rBase      aMean_rMod2      aMean_rMean];

        h_heatmap = figure;
        axisLabel = {'Base',mod2,'Mean'};
    case 'V'

        aBase_rBase = sum(m1.mask_base)/nNeurons;
        aBase_rMean = sum(m1.mask_mean)/nNeurons;

        aMod2_rBase = sum(m2.mask_base)/nNeurons;
        aMod2_rMean = sum(m2.mask_mean)/nNeurons;

        aMean_rBase = sum(m3.mask_base)/nNeurons;
        aMean_rMean = sum(m3.mask_mean)/nNeurons;

        mat_for_hmap = [aBase_rBase      aBase_rMean;
                        aMean_rBase      aMean_rMean];

        h_heatmap = figure;
        axisLabel = {'Base',mod2,'Mean'};
end
cmap_toUse = cmap_customColors(64, 'whiteBlue');
[hImage, hText, hTick] = heatmap_AD(mat_for_hmap, axisLabel, axisLabel, '%0.2f', ...
    'Colormap', cmap_toUse, ...
    'ShowAllTicks', true, ...
    'UseFigureColormap', false, ...
    'Colorbar', true, ...
    'FontSize', 10, ...
    'MinColorValue', 0, ...
    'MaxColorValue', 1, ...
    'GridLines', '-');
xlabel('Recovered model')
ylabel('True model')
set(gca,'tickdir','out')
title(model)

username = getenv('USERNAME');

%% recover parameters
param_struct = struct();

% RPE neurons
n_RPE = os_RPE(bm_RPE_aRPE.mask_base);
param_struct.RPE.alpha.actual = [];
param_struct.RPE.alpha.recovered = [];
param_struct.RPE.slope.actual = [];
param_struct.RPE.slope.recovered = [];
param_struct.RPE.int.actual = [];
param_struct.RPE.int.recovered = [];

for n = 1:length(n_RPE)
    % alpha
    param_struct.RPE.alpha.actual = [param_struct.RPE.alpha.actual ...
                                             n_RPE(n).params.alpha];
    param_struct.RPE.alpha.recovered = [param_struct.RPE.alpha.recovered ...
                                                n_RPE(n).mod_RPE.base.bestParams(1)];
    % slope
    param_struct.RPE.slope.actual = [param_struct.RPE.slope.actual ...
                                             n_RPE(n).params.slope];
    param_struct.RPE.slope.recovered = [param_struct.RPE.slope.recovered ...
                                                n_RPE(n).mod_RPE.base.bestParams(2)];
    % slope
    param_struct.RPE.int.actual = [param_struct.RPE.int.actual ...
                                             n_RPE(n).params.int];
    param_struct.RPE.int.recovered = [param_struct.RPE.int.recovered ...
                                                n_RPE(n).mod_RPE.base.bestParams(3)];
end

% plot it
binEdges = -1.1:0.2:1.1;

h_paramRecovery = figure;
h(1) = subplot(131); hold on
rec_alpha = param_struct.RPE.alpha.actual - param_struct.RPE.alpha.recovered;
histogram(rec_alpha, binEdges, 'EdgeColor','none','normalization','probability')
xlabel('$\alpha$ (actual - recovered)','interpreter','latex')

h(2) = subplot(132); hold on
rec_slope = param_struct.RPE.slope.actual - param_struct.RPE.slope.recovered;
histogram(rec_slope, binEdges, 'EdgeColor','none','normalization','probability')
xlabel('slope (actual - recovered)','interpreter','latex')

h(3) = subplot(133); hold on
rec_int = param_struct.RPE.int.actual - param_struct.RPE.int.recovered;
histogram(rec_int, binEdges, 'EdgeColor','none','normalization','probability')
xlabel('intercept (actual - recovered)','interpreter','latex')

for curr_h = h
    subplot(curr_h)
    ylim_range = get(curr_h, 'YLim');
    plot([0 0],ylim_range,'--','Color', myColors.gray)
    ylabel('Probability')
    set(curr_h,'tickdir','out')
end

saveFigureIteration_ottBari2019(h_paramRecovery, saveLoc, 'recovery_paramBias','FigureSize','max')