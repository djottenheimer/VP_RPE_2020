function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_gRW_RPE_mean(startValues, spikeRate, rewards, timeLocked)

probSpike = normpdf(spikeRate, mean(spikeRate), std(spikeRate));

mean_predictedSpikes = mean(spikeRate);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;