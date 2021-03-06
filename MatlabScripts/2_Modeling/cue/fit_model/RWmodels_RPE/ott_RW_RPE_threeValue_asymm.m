function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_threeValue_asymm(startValues, spikeCounts, rewards, timeLocked, cueInfo)

% cued experiment
alphaPPE = startValues(1);
alphaNPE = startValues(2);
slope = startValues(3);
intercept = startValues(4);
Vinit = alphaPPE / (alphaPPE + alphaNPE);
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
        if RPE(t) >= 0
            V(t + 1) = V(t) + alphaPPE*RPE(t);
        else
            V(t + 1) = V(t) + alphaNPE*RPE(t);
        end
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