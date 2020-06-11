clear; clc
Vinit = 0.5;
alpha = 0.7;
blockLen = 5;
tTotal = blockLen*2;
nRepeat_max = 1e5;

myColors = importColors_bb();
c_suc = myColors.bluishGreen;
c_mal = myColors.blue;

h_sim = figure;
h_random = subplot(131); hold on
h_suc_first = subplot(132); hold on
h_mal_first = subplot(133); hold on


% random
V = NaN(nRepeat_max, blockLen + 2);
RPE = NaN(nRepeat_max, blockLen + 1);
for nRepeat = 1:nRepeat_max
    V(nRepeat,1) = Vinit;
    for t = 1:blockLen + 1
        r = binornd(1, 0.5);
        RPE(nRepeat, t) = r - V(nRepeat, t);
        V(nRepeat, t+1) = V(nRepeat, t) + alpha*RPE(nRepeat, t);
    end
end

subplot(h_random)
suc_trials = RPE > 0;
mal_trials = RPE < 0;
suc_RPE = RPE;
mal_RPE = RPE;
suc_RPE(~suc_trials) = NaN;
mal_RPE(~mal_trials) = NaN;
plot(1:5, nanmean(suc_RPE(:,2:end)), 'linewidth',2,'Color',myColors.orange)
plot(1:5, nanmean(mal_RPE(:,2:end)), 'linewidth',2,'Color',myColors.reddishPurple)

clear V RPE
% sucrose first
V(1) = Vinit;
for t = 1:tTotal
    if t <= blockLen
        r = 1;
    else
        r = 0;
    end
    RPE(t) = r - V(t);
    V(t+1) = V(t) + alpha*RPE(t);
end

subplot(h_suc_first); hold on
plot(1:blockLen, RPE(1:blockLen), 'linewidth', 2, 'Color', c_suc)
plot(blockLen+1:length(RPE), RPE(blockLen+1:end), 'linewidth', 2, 'Color', c_mal)

clear V RPE
% mal first
V(1) = Vinit;
for t = 1:tTotal
    if t <= blockLen
        r = 0;
    else
        r = 1;
    end
    RPE(t) = r - V(t);
    V(t+1) = V(t) + alpha*RPE(t);
end

subplot(h_mal_first); hold on
plot(1:blockLen, RPE(1:blockLen), 'linewidth', 2, 'Color', c_mal)
plot(blockLen+1:length(RPE), RPE(blockLen+1:end), 'linewidth', 2, 'Color', c_suc)

% clean up figure
for cp = [h_random h_suc_first h_mal_first]
    subplot(cp)
    ylim([-1 1])
    set(cp,'tickdir','out')
    ylabel('RPE')
    switch cp
        case h_random
            title('Random')
            legend('Sucrose','Maltodextrin')
        case h_suc_first
            title('Blocked (sucrose first)')
            legend('Sucrose','Maltodextrin')
        case h_mal_first
            title('Blocked (maltodextrin first)')
    end
    xlabel('Session progress')
end