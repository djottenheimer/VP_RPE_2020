clear; clc
task = 'intBlocks';
dir_MLEfit_ott(task)
ott_data = loadData_ott(task);

nStart = 10;
RNG_val = 1;
%%
clear os_temp os

models_of_interest_RPE = {'base','base_flipped','base_asymm','base_asymm_flipped','base_rect', ...
                          'adapt','habit','habit_asymm','curr','curr_flipped','mean'};
models_of_interest_V = {'base','base_asymm','mean'};
all_fits = struct(); % initialize an empty structure 
for ind = 1:length(ott_data.RDHz) % for all neurons
    fprintf('n %i of %i\n', ind, length(ott_data.RDHz))
    
    os_temp(ind).Blocks = ott_data.Blocks(ind);
    os_temp(ind).Blocks12 = ott_data.Blocks12(ind);
    os_temp(ind).Region = ott_data.Region{ind};
    
    os_temp(ind).spikeCount_RD = round(ott_data.RDHz{ind}*1.2); % RD period; 1.2s long
    os_temp(ind).spikeCount_cue = round(ott_data.CueHz{ind}*0.75); % cue period; 0.75s long
    os_temp(ind).spikeCount_PE = ott_data.PEHz{ind}; % PE period; 1.0s long
    
    os_temp(ind).rewards = ott_data.AllTrials{ind}(:, 1); % 0 mal, 1 suc
    os_temp(ind).timeLocked = logical(ott_data.AllTrials{ind}(:, 2)); % trials fast enough to have time-locked responses
    
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
