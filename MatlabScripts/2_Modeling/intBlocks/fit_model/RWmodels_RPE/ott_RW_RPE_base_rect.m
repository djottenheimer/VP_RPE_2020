function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_base_rect(startValues, spikeCounts, rewards, timeLocked)

alphaLearn = startValues(1);
RPErect = startValues(2); % if RPE is below this amount, fix it
slope = startValues(3);
intercept = startValues(4);
Vinit = 0.5;

trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);
RPE_forRateParam = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    V(t + 1) = V(t) + alphaLearn*RPE(t);
    if RPE(t) < RPErect % if we need to rectify this RPE
        RPE_forRateParam(t) = RPErect;
    else
        RPE_forRateParam(t) = RPE(t);
    end
end
rateParam = exp(slope*RPE_forRateParam + intercept);

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