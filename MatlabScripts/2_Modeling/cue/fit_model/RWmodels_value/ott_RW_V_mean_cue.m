function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_mean_cue(startValues, spikeCounts, rewards, timeLocked, cueInfo)

V_sucCue = startValues(1);
V_malCue = startValues(2);

trials = length(rewards);

nonpred_trials = cueInfo(:, 3) == 1; % trials without predictive cues
nonpred_trials = nonpred_trials(timeLocked); % subselect the ones with responses
mean_nonpredSpikes = mean(spikeCounts(nonpred_trials)); % mean firing of non-predicted trials
log_mean_nonpredSpikes = log(mean_nonpredSpikes); % log transform so exp(log...) = mean_nonpredSpikes

rateParam = zeros(trials, 1);

for t = 1:trials
    if cueInfo(t, 1) == 1 % sucrose cue
        rateParam(t) = exp(log_mean_nonpredSpikes + V_sucCue);
    elseif cueInfo(t, 2) == 1 % malto cue
        rateParam(t) = exp(log_mean_nonpredSpikes + V_malCue);
    elseif cueInfo(t, 3) == 1 % nonpredictive cue
        rateParam(t) = exp(log_mean_nonpredSpikes);
    else
        error('cueInfo is 0 for all columns\n');
    end
end

probSpike = poisspdf(spikeCounts, rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;