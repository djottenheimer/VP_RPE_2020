function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_gRW_V_base(startValues, spikeRate, rewards, timeLocked)

alphaLearn = startValues(1);
slope = startValues(2);
intercept = startValues(3);
sigma = startValues(4);
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
rateParam = slope*V(1:trials) + intercept;
probSpike = normpdf(spikeRate, rateParam(timeLocked), sigma); % normal distribution

mean_predictedSpikes = rateParam(timeLocked);
V = V(1:trials);
V = V(timeLocked);
RPE = RPE(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end