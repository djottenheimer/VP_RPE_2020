clear; clc
task = 'cue';
dir_MLEfit_ott(task)
ott_data = loadData_ott(task);

nStart = 10;
RNG_val = 1;

% just the included neurons
my_n = ott_data.IncludedNeurons;
ott_data.Ninfo = {ott_data.Ninfo{my_n,1}; ott_data.Ninfo{my_n,2}}';
ott_data.RDHz = ott_data.RDHz(my_n);
ott_data.PEHz = ott_data.PEHz(my_n);
ott_data.CueHz = ott_data.CueHz(my_n);
ott_data.Predictors = ott_data.Predictors(my_n);
ott_data.AllTrials = ott_data.AllTrials(my_n);
ott_data.PredCue = ott_data.PredCue(my_n);
ott_data.PredCueAllTrials = ott_data.PredCueAllTrials(my_n);
ott_data.Day = ott_data.Day(my_n);
ott_data.Rat = ott_data.Rat(my_n);
ott_data.IncludedNeurons = ott_data.IncludedNeurons(my_n);
%%
clear os_temp os

models_of_interest_RPE = {'base','base_asymm','base_cue','base_asymm_cue',...
                          'threeValue','threeValue_asymm','threeValue_cue','threeValue_asymm_cue',...
                          'curr','curr_cue','mean','mean_cue'};
models_of_interest_V = {'base','base_asymm','base_cue','base_asymm_cue','mean','mean_cue'};

% models_of_interest_RPE = {'base_asymm','base_asymm_cue',...
%                           'threeValue_asymm','threeValue_asymm_cue'};
% models_of_interest_V = {'base_asymm','base_asymm_cue'};

all_fits = struct(); % initialize an empty structure 
tic
for ind = 1:length(ott_data.AllTrials) % for all neurons
    fprintf('n %i of %i. Elapsed time is %0.2f min\n', ind, length(ott_data.AllTrials), toc/60)
    
    os_temp(ind).include = ott_data.IncludedNeurons(ind);
    os_temp(ind).Ninfo = ott_data.Ninfo(ind, :);
    os_temp(ind).day = ott_data.Day(ind);
    os_temp(ind).rat = ott_data.Rat(ind);
    os_temp(ind).spikeCount_RD = round(ott_data.RDHz{ind}*1.2); % RD period; 1.2s long
    os_temp(ind).spikeCount_cue = round(ott_data.CueHz{ind}*0.75); % cue period; 0.75s long
    os_temp(ind).spikeCount_PE = ott_data.PEHz{ind};
    
    os_temp(ind).rewards = ott_data.AllTrials{ind}(:, 1); % 0 mal, 1 suc
    os_temp(ind).timeLocked = logical(ott_data.AllTrials{ind}(:, 2)); % trials fast enough to have time-locked responses
    os_temp(ind).cueInfo = ott_data.PredCueAllTrials{ind}; % predictive cues logical mask [suc mal none]
    
    % spikeCount is a temporary field
    % fit RD
    os_temp(ind).spikeCount = os_temp(ind).spikeCount_RD;
    ms = helper_RW_RPE(os_temp(ind), 'StartingPoints', nStart, 'RNG', RNG_val, 'ParticularModel', models_of_interest_RPE);
    os_temp(ind).mod_RD = ms;
    
    % fit cue
    os_temp(ind).spikeCount = os_temp(ind).spikeCount_cue;
    ms = helper_RW_V(os_temp(ind), 'StartingPoints', nStart, 'RNG', RNG_val, 'ParticularModel', models_of_interest_V);
    os_temp(ind).mod_cue = ms;
    
    % remove spikeCount to avoid future confusion
    os(ind) = rmfield(os_temp(ind), 'spikeCount');
end

fprintf('Finished\n')
save_MLEfit_ott(task, os);

%%
for n = 1:length(os)
    os(n).mod_RD.base_asymm = os_tmp(n).mod_RD.base_asymm;
    os(n).mod_RD.base_asymm_cue = os_tmp(n).mod_RD.base_asymm_cue;
    os(n).mod_RD.threeValue_asymm = os_tmp(n).mod_RD.threeValue_asymm;
    os(n).mod_RD.threeValue_asymm_cue = os_tmp(n).mod_RD.threeValue_asymm_cue;
    os(n).mod_cue.base_asymm = os_tmp(n).mod_cue.base_asymm;
    os(n).mod_cue.base_asymm_cue = os_tmp(n).mod_cue.base_asymm_cue;
end