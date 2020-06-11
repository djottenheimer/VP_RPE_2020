clear; clc
task = 'intBlocks';
ott_data = loadData_ott(task);

nStart = 10;
RNG_val = 1;
my_root = fullfile(ottBari2020_root, 'MatlabScripts', '2_Modeling', 'intBlocks', 'fit_model_gauss');
addpath(fullfile(my_root, 'gRWmodels_RPE'))
addpath(fullfile(my_root, 'gRWmodels_value'))
%%
clear os_temp os

models_of_interest_RPE = {'base','base_asymm','curr','mean'};
% models_of_interest_V = {'base','base_asymm','mean'};
all_fits = struct(); % initialize an empty structure 
for ind = 1:length(ott_data.RDHz) % for all neurons
    fprintf('n %i of %i\n', ind, length(ott_data.RDHz))
    
    os_temp(ind).Blocks = ott_data.Blocks(ind);
    os_temp(ind).Blocks12 = ott_data.Blocks12(ind);
    os_temp(ind).Region = ott_data.Region{ind};
    
    os_temp(ind).spikeRate_RD = zscore(round(ott_data.RDHz{ind}*1.2)); % RD period; 1.2s long
    os_temp(ind).spikeRate_cue = zscore(round(ott_data.CueHz{ind}*0.75)); % cue period; 0.75s long
    
    os_temp(ind).rewards = ott_data.AllTrials{ind}(:, 1); % 0 mal, 1 suc
    os_temp(ind).timeLocked = logical(ott_data.AllTrials{ind}(:, 2)); % trials fast enough to have time-locked responses
    
    % spikeRate is a temporary field
    % fit RD
    os_temp(ind).spikeRate = os_temp(ind).spikeRate_RD;
    ms = helper_gRW_RPE(os_temp(ind), 'StartingPoints', nStart, 'RNG', RNG_val, 'ParticularModel', models_of_interest_RPE);
    os_temp(ind).mod_RD = ms;
    
%     % fit cue
%     os_temp(ind).spikeRate = os_temp(ind).spikeRate_cue;
%     ms = helper_gRW_V(os_temp(ind), 'StartingPoints', nStart, 'RNG', RNG_val, 'ParticularModel', models_of_interest_V);
%     os_temp(ind).mod_cue = ms;
    
    % remove spikeRate to avoid future confusion
    os(ind) = rmfield(os_temp(ind), 'spikeRate');
end

fprintf('Finished\n')
C:\Users\bilal\Dropbox\RPE Manuscript copy\Data\Modeling\ModelFits
save(fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits', 'intBlocks_MLEfits_gauss.mat'),'os')