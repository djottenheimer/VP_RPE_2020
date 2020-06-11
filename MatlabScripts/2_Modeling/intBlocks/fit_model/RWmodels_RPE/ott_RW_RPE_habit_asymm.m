function [LH, probSpike, V, mean_predictedSpikes, RPE, slope_vec] = ott_RW_RPE_habit_asymm(startValues, spikeCounts, rewards, timeLocked)

alphaPPE = startValues(1);
alphaNPE = startValues(2);
slope_max = startValues(3); % slope positive
slope_min = startValues(4);
int = startValues(5);
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
slope_vec = (slope_min - slope_max)*V(1:trials) + slope_min; % slope modulated by value
% right side (log_mean_spikes - slope_vec*0.5) constrains firing rate to mean firing rate at r = 0.5
rateParam = exp(slope_vec.*rewards + int - slope_vec*0.5); % firing rate modulated by rewards
probSpike = poisspdf(spikeCounts, rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);
V = V(1:trials);
V = V(timeLocked);
RPE = RPE(timeLocked);
slope_vec = slope_vec(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end