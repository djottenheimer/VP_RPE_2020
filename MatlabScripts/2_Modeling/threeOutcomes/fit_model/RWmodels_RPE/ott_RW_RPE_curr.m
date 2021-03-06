function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_RPE_curr(startValues, spikeCounts, rewards, timeLocked)
% firing rate only correlates with current reward; no real learning happens

% reward is 0: mal, 1: suc, 2: water

slope = startValues(1);
intercept = startValues(2);
rho = startValues(3); % how valuable is maltodextrin on a water -> malto -> sucrose scale

water_ind = rewards == 2;
mal_ind = rewards == 0;
rewards(water_ind) = 0;
rewards(mal_ind) = rho; % scale mal between 0 and 1

rateParam = exp(slope*rewards + intercept);

probSpike = poisspdf(spikeCounts, rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;