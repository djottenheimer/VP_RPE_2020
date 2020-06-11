function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_threeValue(startValues, spikeCounts, rewards, timeLocked, cueInfo)

% cued experiment
alphaLearn = startValues(1);
slope = startValues(2);
intercept = startValues(3);
Vinit = 0.5;
Vsuc = 1;
Vmal = 0;

trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);
RPEforObs = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials
    if cueInfo(t, 1) == 1 % sucrose cue
        RPEforObs(t) = rewards(t) - Vsuc;
        V(t + 1) = V(t); % carry forward the value function for the nonpredictive cue
    elseif cueInfo(t, 2) == 1 % malto cue
        RPEforObs(t) = rewards(t) - Vmal;
        V(t + 1) = V(t); % carry forward the value function
    elseif cueInfo(t, 3) == 1 % nonpredictive cue
        RPE(t) = rewards(t) - V(t);
        V(t + 1) = V(t) + alphaLearn*RPE(t);
        RPEforObs(t) = RPE(t);
    end
end
rateParam = exp(slope*RPEforObs + intercept);

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