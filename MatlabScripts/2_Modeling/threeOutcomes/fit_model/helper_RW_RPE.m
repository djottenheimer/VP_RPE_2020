function ms = helper_RW_RPE(os, varargin)
% helper_RW_RPE    Fits neural models to test RW model predictions ~ RPE
%   ms = helper_RW_RPE(os, varargin)
%   INPUTS
%       os: behavioral data structure
%           .spikeCount: number of spikes in a particular time bin
%           .rewards: 0 for maltodextrin, 1 for sucrose
%           .timeLocked: 1 if the animal licked within 2s, 0 if the animal did not (logical)
%               timeLocked will always be more than or equal to the number of spike trials
%       varargin
%           StartingPoints: determines how many points to optimize from
%           ParticularModel: cell array of strings of models to use
%           RNG: random number generator seed (default = 1)
%   OUTPUTS
%       ms: model structure of fits

p = inputParser;
p.addParameter('StartingPoints', 1)
p.addParameter('ParticularModel', []);
p.addParameter('RNG', []);
p.parse(varargin{:});

if ~isempty(p.Results.RNG)
    rng(p.Results.RNG)
end
              
% Initialize models
if isempty(p.Results.ParticularModel)
    modelNames = {'base', 'curr', 'mean'};
else
    modelNames = p.Results.ParticularModel;
end


% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

% set boundary conditions
alpha_range = [0 1];
slope_range = [0 20]; % reward sensitivity
intercept_range = [-20 20]; % baseline spiking
rho_range = [0 1];

for currMod = modelNames
    currMod = currMod{:};
    
    % initialize output variables
    runs = p.Results.StartingPoints;
    LH = zeros(runs, 1);
    
    if strcmp(currMod, 'base')
        paramNames = {'alpha','slope','intercept','rho'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1) ...
                       rand(runs, 1)*diff(rho_range) + rho_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  slope_range(2);  intercept_range(2);  rho_range(2);
           -alpha_range(1); -slope_range(1); -intercept_range(1); -rho_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@ott_RW_RPE_base, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.rewards, os.timeLocked);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            ott_RW_RPE_base(ms.(currMod).bestParams, os.spikeCount, os.rewards, os.timeLocked);
    elseif strcmp(currMod, 'base_asymm')
        paramNames = {'alphaPPE','alphaNPE','slope','intercept','rho'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1) ...
                       rand(runs, 1)*diff(rho_range) + rho_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  alpha_range(2);  slope_range(2);  intercept_range(2);  rho_range(2);
           -alpha_range(1); -alpha_range(1); -slope_range(1); -intercept_range(1); -rho_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@ott_RW_RPE_base_asymm, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.rewards, os.timeLocked);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            ott_RW_RPE_base_asymm(ms.(currMod).bestParams, os.spikeCount, os.rewards, os.timeLocked);
    elseif strcmp(currMod, 'curr')
        paramNames = {'slope','intercept','rho'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1) ...
                       rand(runs, 1)*diff(rho_range) + rho_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);  rho_range(2);
           -slope_range(1); -intercept_range(1); -rho_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@ott_RW_RPE_curr, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.rewards, os.timeLocked);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            ott_RW_RPE_curr(ms.(currMod).bestParams, os.spikeCount, os.rewards, os.timeLocked);
    elseif strcmp(currMod, 'mean')
        paramNames = {''};
        numParam = 0;        
        [LH, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            ott_RW_RPE_mean([], os.spikeCount);
        bestFit = 1;
        ms.(currMod).bestParams = [];
    else 
        error('RW model: Model name not found')
    end
    ms.(currMod).paramNames = paramNames;
    ms.(currMod).LH = -1 * LH(bestFit, :);
    ms.(currMod).BIC = log(length(os.spikeCount))*numParam  - 2*ms.(currMod).LH;
    ms.(currMod).AIC = 2*numParam - 2*ms.(currMod).LH;
    ms.(currMod).AICc = ms.(currMod).AIC + (2*numParam^2 + 2*numParam)/(length(os.spikeCount) - numParam - 1);
end