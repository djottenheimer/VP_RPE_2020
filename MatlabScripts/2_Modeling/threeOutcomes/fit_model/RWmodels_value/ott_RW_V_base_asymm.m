function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_base_asymm(startValues, spikeCounts, rewards, timeLocked)

% reward is 0: mal, 1: suc, 2: water

alphaPPE = startValues(1);
alphaNPE = startValues(2);
slope = startValues(3);
intercept = startValues(4);
rho = startValues(5); % how valuable is maltodextrin on a water -> malto -> sucrose scale
Vinit = 2/3*(alphaPPE / (alphaPPE + alphaNPE)) + 1/3*rho;

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
    if RPE(t) >= 0
        V(t + 1) = V(t) + alphaPPE*RPE(t);
    else
        V(t + 1) = V(t) + alphaNPE*RPE(t);
    end
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