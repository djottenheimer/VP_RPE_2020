function [LH, probSpike, V, mean_predictedSpikes, RPE, int_vec] = ott_RW_RPE_adapt(startValues, spikeCounts, rewards, timeLocked)

alphaLearn = startValues(1);
slope = startValues(2);
int_max = startValues(3);
int_min = startValues(4);
Vinit = 0.5;

trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    V(t + 1) = V(t) + alphaLearn*RPE(t);
end
int_vec = (int_min - int_max)*V(1:trials) + int_min; % intercept modulated by value
rateParam = exp(slope*rewards + int_vec); % firing rate modulated by rewards
probSpike = poisspdf(spikeCounts, rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);
V = V(1:trials);
V = V(timeLocked);
RPE = RPE(timeLocked);
int_vec = int_vec(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end