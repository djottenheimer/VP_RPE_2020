clear; clc
load(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'intBlocks_MLEfits.mat'));
int_task = [os.Blocks] == 0;
os = os(int_task);
VP_mask = contains({os.Region}, 'VP');
os = os(VP_mask);
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
                       
%%
mask_cueV = bm_cue.mask_base;
mask_rdRPE = bm_RD.mask_base;

n_cueV_rdRPE = sum(mask_cueV(mask_rdRPE));
n_cueV_NOTrdRPE = sum(mask_cueV(~mask_rdRPE));

n_rdRPE = sum(mask_rdRPE);
n_NOTrdRPE = sum(~mask_rdRPE);

[~,p] = prop_test([n_cueV_rdRPE n_cueV_NOTrdRPE], [n_rdRPE n_NOTrdRPE]);

fprintf('\n------\n')
fprintf('%i of %i RPE neurons (%0.2f%%) have cue-value response\n', n_cueV_rdRPE, n_rdRPE, n_cueV_rdRPE/n_rdRPE*100);
fprintf('%i of %i non-RPE neurons (%0.2f%%) have cue-value responses \n', n_cueV_NOTrdRPE, n_NOTrdRPE, n_cueV_NOTrdRPE/n_NOTrdRPE*100);
fprintf('pValue: %0.2e\n', p)