function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_gRW_V_mean(startValues, spikeRate, rewards, timeLocked)

probSpike = normpdf(spikeRate, mean(spikeRate), std(spikeRate)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = mean(spikeRate);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;