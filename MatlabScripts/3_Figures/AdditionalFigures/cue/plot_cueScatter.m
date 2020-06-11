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
                       
%%
timePeriod = 'cue';
switch timePeriod
    case 'RD'
        bm = bm_RD;
        models_to_cycle = {'base_cue','curr_cue','mean_cue'};
        signFlip = -1; % flip to keep consistent with cue-value
    case 'cue'
        bm = bm_cue;
        models_to_cycle = {'base_cue','mean_cue'};
        signFlip = 1;
end
valueStruct = struct();
valueStruct.vsuc = [];
valueStruct.vmal = [];
for m = models_to_cycle
    m = m{:};
    mOfInt = bm.(['mask_' m]);
    if any(mOfInt)
        valueStruct.(['vsuc_' m]) = getParameters_ott(os(mOfInt),timePeriod,m,'vsuc');
        valueStruct.(['vmal_' m]) = getParameters_ott(os(mOfInt),timePeriod,m,'vmal');
        valueStruct.vsuc = [valueStruct.vsuc getParameters_ott(os(mOfInt),timePeriod,m,'vsuc')];
        valueStruct.vmal = [valueStruct.vmal getParameters_ott(os(mOfInt),timePeriod,m,'vmal')];
    else
        valueStruct.(['vsuc_' m]) = [];
        valueStruct.(['vmal_' m]) = [];
    end
end

valueStruct.vsuc = signFlip*valueStruct.vsuc;
valueStruct.vmal = signFlip*valueStruct.vmal;
h_scatter = figure; hold on
scatter(valueStruct.vsuc, valueStruct.vmal)
plot([-1 1],[0 0],'Color',myColors.gray)
plot([0 0],[-1 1],'Color',myColors.gray)
xlim([-1 1])
ylim([-1 1])
xlabel('Sucrose value','interpreter','latex')
ylabel('Maltodextrin value','interpreter','latex')
title(timePeriod)
set(gca,'tickdir','out','xtick',-1:0.5:1,'ytick',-1:0.5:1)

nTot = numel(valueStruct.vsuc);
q1 = sum(valueStruct.vsuc > 0 & valueStruct.vmal > 0)/nTot * 100;
q2 = sum(valueStruct.vsuc < 0 & valueStruct.vmal > 0)/nTot * 100;
q3 = sum(valueStruct.vsuc < 0 & valueStruct.vmal < 0)/nTot * 100;
q4 = sum(valueStruct.vsuc > 0 & valueStruct.vmal < 0)/nTot * 100;
strFormat = '%0.1f';
text( 0.5,  0.5, [num2str(q1, strFormat) '%'])
text(-0.5,  0.5, [num2str(q2, strFormat) '%'])
text(-0.5, -0.5, [num2str(q3, strFormat) '%'])
text( 0.5, -0.5, [num2str(q4, strFormat) '%'])

p_prop = myBinomTest(sum(valueStruct.vsuc > 0 & valueStruct.vmal < 0), nTot, 0.25, 'two');
fprintf('p_value for quadrant 4 is %0.1e\n', p_prop)

axis square