clear; clc
load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'intBlocks_MLEfits.mat'));
int_task = [os.Blocks] == 0;
os = os(int_task);
VP_mask = contains({os.Region}, 'VP');
os = os(VP_mask);

%% get relevant behavior models
modelCriterion = 'AIC';
plotFlag = true;

models_of_interest_RPE = {'base','adapt','habit','curr','mean'};

timePeriod = 'RD';
bm_RD = select_RPEmods(os, timePeriod,'scoreToUse',modelCriterion,...
                          'plotModels_Flag',plotFlag,...
                          'particularModel',models_of_interest_RPE);


%% pValues
nTot = numel(bm_RD.mask_adapt);
nBase = sum(bm_RD.mask_base);
nAdapt = sum(bm_RD.mask_adapt);
nHabit = sum(bm_RD.mask_habit);

[~,pAdapt] = prop_test([nBase nAdapt],[nTot nTot]);
[~,pHabit] = prop_test([nBase nHabit],[nTot nTot]);

fprintf('\n------\n')
fprintf('Base: %i\nAdapt: %i\nHabit: %i\nTotal: %i\n',nBase,nAdapt,nHabit,nTot);
fprintf('pValue for adapt vs base: %0.2e\n', pAdapt);
fprintf('pValue for habit vs base: %0.2e\n', pHabit);