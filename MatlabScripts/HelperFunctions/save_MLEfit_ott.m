function save_MLEfit_ott(task, os, varargin)
% task should be a string
%   standard: 50% suc/malto
%   threeOutcomes: suc/malto/water
%   cued: suc/mal-predictive cue
p = inputParser;
p.addParameter('Simulated_Flag', false)
p.parse(varargin{:});

master_root = ottBari2020_root();

saveDir = fullfile(master_root, 'Data', 'Modeling', 'ModelFits');
if p.Results.Simulated_Flag == false
    save(fullfile(saveDir, [task '_MLEfits.mat']), 'os')
    fprintf('Saved os to %s\n', saveDir);
else
    save(fullfile(saveDir, [task '_MLEfits_simulated.mat']), 'os')
    fprintf('Saved simulated os to %s\n', saveDir);
end