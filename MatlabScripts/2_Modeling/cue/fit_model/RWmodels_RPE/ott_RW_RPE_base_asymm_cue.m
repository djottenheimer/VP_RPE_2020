function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_base_asymm_cue(startValues, spikeCounts, rewards, timeLocked, cueInfo)

% cued experiment

alphaPPE = startValues(1);
alphaNPE = startValues(2);
slope = startValues(3);
intercept = startValues(4);
V_sucCue = startValues(5);
V_malCue = startValues(6);
Vinit = alphaPPE / (alphaPPE + alphaNPE);


trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);
rateParam = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    if RPE(t) >= 0
        V(t + 1) = V(t) + alphaPPE*RPE(t);
    else
        V(t + 1) = V(t) + alphaNPE*RPE(t);
    end
    
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