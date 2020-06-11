function dir_MLEfit_ott(task)
% task should be a string
%   standard: 50% suc/malto
%   threeOutcomes: suc/malto/water
%   cued: suc/mal-predictive cue

master_root = ottBari2020_root();

allTasks = {'intBlocks','threeOutcomes','cue'};

taskMask = contains(allTasks, task);
if all(~taskMask) == 1 % task no found
    error('task name not found')
else
    % change directory to appropriate directory
    cd(fullfile(master_root, 'MatlabScripts', '2_Modeling', allTasks{taskMask}, 'fit_model'))
    % add path of appropriate model
    addpath(fullfile(master_root, 'MatlabScripts', '2_Modeling', allTasks{taskMask}, 'fit_model', 'RWmodels_RPE'))
    addpath(fullfile(master_root, 'MatlabScripts', '2_Modeling', allTasks{taskMask}, 'fit_model', 'RWmodels_value'))
    % remove all other paths
    for nonTaskInd = find(~taskMask)
        rpe_path = fullfile(master_root, 'MatlabScripts', '2_Modeling', allTasks{nonTaskInd}, ...
                            'fit_model', 'RWmodels_RPE');
        if contains(path, rpe_path)
            rmpath(rpe_path)
        end
        
        value_path = fullfile(master_root, 'MatlabScripts', '2_Modeling', allTasks{nonTaskInd}, ...
                              'fit_model', 'RWmodels_value');
        if contains(path, value_path)
            rmpath(value_path)
        end
    end
end