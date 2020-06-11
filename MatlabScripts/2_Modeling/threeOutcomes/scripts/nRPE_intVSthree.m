clear; clc
model_root = fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits');

% interspersed task
task = 'intBlocks';
load(fullfile(model_root, [task '_MLEfits.mat']))
int_task = [os.Blocks] == 0;
os = os(int_task);
VP_mask = contains({os.Region}, 'VP');
os = os(VP_mask);

os_int = os;
clear os

task = 'threeOutcomes';
load(fullfile(model_root, [task '_MLEfits.mat']))
os_three = os;
clear os

% get relevant behavior models
modelCriterion = 'AIC';
plotFlag = false;

models_of_interest_RPE = {'base','curr','mean'};

timePeriod = 'RD';
bm_int_RD = select_RPEmods(os_int, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag, ...
                               'particularModel', models_of_interest_RPE);
bm_three_RD = select_RPEmods(os_three, timePeriod,'scoreToUse',modelCriterion,'plotModels_Flag',plotFlag, ...
                               'particularModel', models_of_interest_RPE);

%%
nTot_int = numel(bm_int_RD.mask_base);
nRPE_int = sum(bm_int_RD.mask_base);

nTot_three = numel(bm_three_RD.mask_base);
nRPE_three = sum(bm_three_RD.mask_base);

[~,p] = prop_test([nRPE_int nRPE_three],[nTot_int nTot_three]);

fprintf('\n------\n')
fprintf('Int task: %i RPE of %i tot (%0.2f%%)\n',nRPE_int,nTot_int,nRPE_int/nTot_int*100);
fprintf('Three task: %i RPE of %i tot (%0.2f%%)\n',nRPE_three,nTot_three,nRPE_three/nTot_three*100);
fprintf('pValue: %0.2e\n', p)
fprintf('------\n')