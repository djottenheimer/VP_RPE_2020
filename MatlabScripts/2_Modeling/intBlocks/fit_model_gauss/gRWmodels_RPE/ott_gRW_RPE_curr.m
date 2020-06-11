function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_gRW_RPE_curr(startValues, spikeRate, rewards, timeLocked)
% firing rate only correlates with current reward; no real learning happens

slope = startValues(1);
intercept = startValues(2);
sigma = startValues(3);

rateParam = slope*rewards + intercept;
probSpike = normpdf(spikeRate, rateParam(timeLocked), sigma); % normal distribution

mean_predictedSpikes = rateParam(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;