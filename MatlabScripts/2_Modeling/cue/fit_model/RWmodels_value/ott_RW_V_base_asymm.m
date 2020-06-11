function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_base_asymm(startValues, spikeCounts, rewards, timeLocked)

% cued experiment

alphaPPE = startValues(1);
alphaNPE = startValues(2);
slope = startValues(3);
intercept = startValues(4);
Vinit = alphaPPE / (alphaPPE + alphaNPE);


trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    if RPE(t) >= 0
        V(t + 1) = V(t) + alphaPPE*RPE(t);
    else
        V(t + 1) = V(t) + alphaNPE*RPE(t);
    end
end
rateParam = exp(slope*V(1:trials) + intercept);

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