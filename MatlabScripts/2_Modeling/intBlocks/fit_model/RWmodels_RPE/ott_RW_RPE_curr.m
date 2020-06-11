function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_curr(startValues, spikeCounts, rewards, timeLocked)
% firing rate only correlates with current reward; no real learning happens

slope = startValues(1);
intercept = startValues(2);

rateParam = exp(slope*rewards + intercept);
rateParam(rateParam < 0) = 0.1;

probSpike = poisspdf(spikeCounts, rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;