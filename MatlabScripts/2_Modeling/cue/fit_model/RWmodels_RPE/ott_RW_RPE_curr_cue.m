function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_curr_cue(startValues, spikeCounts, rewards, timeLocked, cueInfo)

% cued experiment; weight determines how much V_sucrose and V_mal influence spike counts
% cueInfo: N x 3 where the left - suc, mid - mal, right - none

slope = startValues(1);
intercept = startValues(2);
V_sucCue = startValues(3);
V_malCue = startValues(4);

alphaLearn = 0; % only current reward can be encoded
Vinit = 0.5;

trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1); % RPE for trial-by-trial learning
rateParam = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    
    V(t + 1) = V(t) + alphaLearn*RPE(t);
    
    if cueInfo(t, 1) == 1 % sucrose cue
        rateParam(t) = exp(slope*RPE(t) + intercept + V_sucCue);
    elseif cueInfo(t, 2) == 1 % malto cue
        rateParam(t) = exp(slope*RPE(t) + intercept + V_malCue);
    elseif cueInfo(t, 3) == 1 % nonpredictive cue
        rateParam(t) = exp(slope*RPE(t) + intercept);
    else
        error('cueInfo is 0 for all columns\n');
    end
end

probSpike = poisspdf(spikeCounts, rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);
V = V(1:trials);
V = V(timeLocked);
RPE = RPE(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end