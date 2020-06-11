function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_base(startValues, spikeCounts, rewards, timeLocked)

% reward is 0: mal, 1: suc, 2: water

alphaLearn = startValues(1);
slope = startValues(2);
intercept = startValues(3);
rho = startValues(4); % how valuable is maltodextrin on a water -> malto -> sucrose scale
Vinit = (1 + 0 + rho)/3;

water_ind = rewards == 2;
mal_ind = rewards == 0;
rewards(water_ind) = 0;
rewards(mal_ind) = rho; % scale mal between 0 and 1

trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    V(t + 1) = V(t) + alphaLearn*RPE(t);
end
rateParam = exp(slope*V(1:trials) + intercept);
rateParam(rateParam < 0) = 0.1; % set rate param to zero if it goes below; might consider a better rule in the future

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