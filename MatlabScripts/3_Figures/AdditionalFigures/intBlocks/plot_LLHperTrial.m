clear; clc

load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'intBlocks_MLEfits.mat'));
int_task = [os.Blocks] == 0;
os = os(int_task);
VP_mask = contains({os.Region}, 'VP');
os = os(VP_mask); % only VP neurons

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

%% likelihood per trial
timePeriod = 'RD'; % RD, cue

LH = struct();
switch timePeriod
    case 'RD'
        mods = {'base','curr','mean'};
        bm = bm_RD;
        mOfInt = 'mod_RD';
    case 'cue'
        mods = {'base','mean'};
        bm = bm_cue;
        mOfInt = 'mod_cue';
end

% generate a blank LH structure
for m1 = mods
    m1 = m1{:};
    for m2 = mods
        m2 = m2{:};
        LH.(m1).(m2) = [];
    end
end

% add for each neuron
LH_cutoff = 10;
for n = 1:length(os)
    nTrials = numel(os(n).spikeCount_RD);
    for m1 = mods
        m1 = m1{:};
        if bm.(['mask_' m1])(n) == 1 % if this is a neuron of interest
            for m2 = mods
                m2 = m2{:};
                tmp_LH = -os(n).(mOfInt).(m2).LH/nTrials;
                if tmp_LH < LH_cutoff
                    LH.(m1).(m2) = [LH.(m1).(m2) tmp_LH];
                else
                    LH.(m1).(m2) = [LH.(m1).(m2) NaN];
                end
            end            
        end
    end
end

% fix names
mods_forLabel = mods;
mods_forLabel = strrep(mods_forLabel,'base','RPE');
mods_forLabel = strrep(mods_forLabel,'curr','Current outcome');
mods_forLabel = strrep(mods_forLabel,'mean','Unmodulated');

% plot it
h_boxplot = figure;
set(h_boxplot, 'units', 'normalized', 'outerposition', [0 0 1 1])
n_plot = length(mods);
for m1_ind = 1:length(mods)
    m1 = mods{m1_ind};
    
    box_mat = [];
    for m2_ind = 1:length(mods)
        m2 = mods{m2_ind};
        box_mat(:,m2_ind) = LH.(m1).(m2);
    end
    
    h_sp_LH(m1_ind) = subplot(2,n_plot,m1_ind); hold on
    boxplot(box_mat,'BoxStyle','filled','PlotStyle','traditional','symbol','')
    set(gca,'tickdir','out','xtick',1:4,'xticklabel',cellfun(@(i) [i ' model'],mods_forLabel,'UniformOutput',false),'xticklabelrotation',60)
    ylabel('LH/trial')
    title([mods_forLabel{m1_ind} ' neurons'],'interpreter','none')
    
    % plot relative to mean model
    box_mat = [];
    for m2_ind = 1:length(mods) - 1 % all except the mean model
        m2 = mods{m2_ind};
        box_mat(:,m2_ind) = LH.(m1).(m2) - LH.(m1).mean;
    end
    h_sp_LHrel(m1_ind) = subplot(2,n_plot,m1_ind + n_plot); hold on
    boxplot(box_mat,'BoxStyle','filled','PlotStyle','traditional','symbol','')
    set(gca,'tickdir','out','xtick',1:3,'xticklabel',cellfun(@(i) [i ' model'],mods_forLabel(1:3),'UniformOutput',false),'xticklabelrotation',60)
    ylabel('\DeltaLH/trial (relative to Unmodulated model)')
    title([mods_forLabel{m1_ind} ' neurons'],'interpreter','none')
    plot([0.5 2.5],[0 0],'k:')
end

pause(0.1)
ylim_LH = [-0.1 max([h_sp_LH.YLim])];
ylim_LHrel = [-1 0.1];
for cp = h_sp_LH
    subplot(cp)
    ylim(ylim_LH)
end
for cp = h_sp_LHrel
    subplot(cp)
    ylim(ylim_LHrel)
end


%%
fprintf('RPE neurons [IQR]: %0.2f [%0.2f - %0.2f]\n', median(LH.base.base), prctile(LH.base.base, 25), prctile(LH.base.base, 75))
fprintf('Current outcome neurons [IQR]: %0.2f [%0.2f - %0.2f]\n', median(LH.curr.curr), prctile(LH.curr.curr, 25), prctile(LH.curr.curr, 75))
fprintf('Unmodulated neurons [IQR]: %0.2f [%0.2f - %0.2f]\n', median(LH.mean.mean), prctile(LH.mean.mean, 25), prctile(LH.mean.mean, 75))