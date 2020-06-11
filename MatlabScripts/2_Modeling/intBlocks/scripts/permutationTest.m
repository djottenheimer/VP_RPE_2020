clear; clc
model_root = fullfile(ottBari2020_root, 'Data', 'Modeling', 'ModelFits');
task = 'intBlocks';
load(fullfile(model_root, [task '_MLEfits.mat']))
os_int = os;
clear os
int_task = [os_int.Blocks] == 0;
os_int = os_int(int_task);
VP_mask = contains({os_int.Region}, 'VP');
os_int = os_int(VP_mask); % only look at VP neurons

task = 'threeOutcomes';
load(fullfile(model_root, [task '_MLEfits.mat']))
os_three = os;
clear os

task = 'cue';
load(fullfile(model_root, [task '_MLEfits.mat']))
os_cue = os;
clear os

%% 
% figure out which sessions are unique
timePeriod = 'RD'; % cue, RD
switch timePeriod
    case 'RD'
        latentVar = 'RPEs';
        spOfInt = 'spikeCount_RD';
        mOfInt = 'mod_RD';
        models_of_interest_int_three = {'base','curr','mean'};
        models_of_interest_cue = {'base','base_cue','curr','curr_cue','mean','mean_cue'};
    case 'cue'
        latentVar = 'V';
        spOfInt = 'spikeCount_cue';
        mOfInt = 'mod_cue';
        models_of_interest_int_three = {'base','mean'};
        models_of_interest_cue = {'base','base_cue','mean','mean_cue'};
end

nTrials = [];
session = 1;
os_int(1).session = 1;
nTrials = [nTrials length(os_int(1).(spOfInt))];
for i = 2:length(os_int)
    nTrials = [nTrials length(os_int(i).(spOfInt))];
    if length(os_int(i).rewards) ~= length(os_int(i-1).rewards) % different number of trials
        session = session + 1;
    elseif ~all(os_int(i).rewards == os_int(i-1).rewards) % same trials but different rewards
        session = session + 1;
    end
    os_int(i).session = session;
end
session = session + 1;
os_three(1).session = session;
nTrials = [nTrials length(os_three(1).(spOfInt))];
for i = 2:length(os_three)
    nTrials = [nTrials length(os_three(i).(spOfInt))];
    if length(os_three(i).rewards) ~= length(os_three(i-1).rewards) % different number of trials
        session = session + 1;
    elseif ~all(os_three(i).rewards == os_three(i-1).rewards) % same trials but different rewards
        session = session + 1;
    end
    os_three(i).session = session;
end
session = session + 1;
os_cue(1).session = session;
nTrials = [nTrials length(os_cue(1).(spOfInt))];
for i = 2:length(os_cue)
    nTrials = [nTrials length(os_cue(i).(spOfInt))];
    if length(os_cue(i).rewards) ~= length(os_cue(i-1).rewards) % different number of trials
        session = session + 1;
    elseif ~all(os_cue(i).rewards == os_cue(i-1).rewards) % same trials but different rewards
        session = session + 1;
    end
    os_cue(i).session = session;
end

modelCriterion = 'AIC';
plotFlag = false;
bm_int = select_RPEmods(os_int, timePeriod,'scoreToUse',modelCriterion,...
                          'plotModels_Flag',plotFlag,...
                          'particularModel',models_of_interest_int_three);
bm_three = select_RPEmods(os_three, timePeriod,'scoreToUse',modelCriterion,...
                          'plotModels_Flag',plotFlag,...
                          'particularModel',models_of_interest_int_three);
bm_cue = select_RPEmods(os_cue,timePeriod,'scoreToUse',modelCriterion,...
                          'plotModels_Flag',plotFlag,...
                          'particularModel',models_of_interest_cue);

% generate permutation structure to make life easier
permStruct = struct();
nInd = 1;
for n = 1:length(os_int)
    session = ['s' int2str(os_int(n).session)];
    neuron = ['n' int2str(nInd)];
    if bm_int.mask_base(n) == 1 % if this is an RPE neuron
        permStruct.(session).(neuron).task = 'int';
        permStruct.(session).(neuron).spikes = os_int(n).(spOfInt);
        permStruct.(session).(neuron).(latentVar) = os_int(n).(mOfInt).base.(latentVar);

        nInd = nInd + 1;
    end
end
for n = 1:length(os_three)
    session = ['s' int2str(os_three(n).session)];
    neuron = ['n' int2str(nInd)];
    if bm_three.mask_base(n) == 1 % if this is an RPE neuron
        permStruct.(session).(neuron).task = 'three';
        permStruct.(session).(neuron).spikes = os_three(n).(spOfInt);
        permStruct.(session).(neuron).(latentVar) = os_three(n).(mOfInt).base.(latentVar);

        nInd = nInd + 1;
    end
end
for n = 1:length(os_cue)
    session = ['s' int2str(os_cue(n).session)];
    neuron = ['n' int2str(nInd)];
    if bm_cue.mask_base(n) == 1 % if this is an RPE neuron
        permStruct.(session).(neuron).task = 'cue';
        permStruct.(session).(neuron).spikes = os_cue(n).(spOfInt);
        permStruct.(session).(neuron).(latentVar) = os_cue(n).(mOfInt).base.(latentVar);

        nInd = nInd + 1;
    elseif bm_cue.mask_base_cue(n) == 1 % if this is an RPE neuron w/ cue response
        permStruct.(session).(neuron).task = 'cue';
        permStruct.(session).(neuron).spikes = os_cue(n).(spOfInt);
        permStruct.(session).(neuron).(latentVar) = os_cue(n).(mOfInt).base_cue.(latentVar);

        nInd = nInd + 1;
    end
end

%% generate distribution of RPEs to generate t-statistics
all_s = fields(permStruct)';
maxTrials = min(nTrials);
for s_ind = 1:length(all_s)
    s = all_s{s_ind};
    fprintf('On session %s\n',s)
    all_n = fields(permStruct.(s))';
    % go through each neuron
    for n = all_n
        n = n{:};
        
        % generate RPEmat for other neurons; used for t-stat
        latentMat = [];
        s_exclude_list = all_s;
        s_exclude_list(s_ind) = [];
        for s_other = s_exclude_list
            s_other = s_other{:};
            
            all_n_other = fields(permStruct.(s_other))';
            rand_n_ind = randi(length(all_n_other)); % select a random neuron
            latentMat = [latentMat; permStruct.(s_other).(all_n_other{rand_n_ind}).(latentVar)(1:maxTrials)'];
        end
        permStruct.(s).(n).([latentVar 'mat']) = latentMat;
        
        % generate a neuron-specific t-distribution
        tDist = [];
        for t_ind = 1:size(permStruct.(s).(n).([latentVar 'mat']), 1)
            fr = normalize(permStruct.(s).(n).spikes(1:maxTrials));
            tDist_lm = fitlm(fr, permStruct.(s).(n).([latentVar 'mat'])(t_ind,:));
            tDist = [tDist tDist_lm.Coefficients.tStat(2)];
        end
        tReal_lm = fitlm(fr, permStruct.(s).(n).(latentVar)(1:maxTrials));
        tReal = tReal_lm.Coefficients.tStat(2);
        permStruct.(s).(n).tDist = tDist;
        permStruct.(s).(n).tReal = tReal;
    end
end

%% determine significant neurons
pCutoff = 0.05;
pCutoff = pCutoff * 100;
all_sig_n = [];
for s = all_s
    s = s{:};
    for n = fields(permStruct.(s))'
        n = n{:};
        
        % determine empirical significance
        threshold = prctile(permStruct.(s).(n).tDist, [pCutoff/2 100-pCutoff/2]);
        if permStruct.(s).(n).tReal >= 0
            if permStruct.(s).(n).tReal > threshold(2)
                all_sig_n = [all_sig_n true];
            else
                all_sig_n = [all_sig_n false];
            end
        else
            if permStruct.(s).(n).tReal < threshold(1)
                all_sig_n = [all_sig_n true];
            else
                all_sig_n = [all_sig_n false];
            end
        end
    end
end

fprintf('%i of %i neurons (%0.1f%%) are significant\n', ...
        sum(all_sig_n), length(all_sig_n), 100*mean(all_sig_n))
% 180 of 196 neurons (91.8%) are significant

%% how correlated is the latent variable of one session with others
all_corr = {};
noise_corr = {};
for s = all_s
    s = s{:};
    all_n = fields(permStruct.(s))';
    
    for n = all_n
        n = n{:};
        
        v1 = permStruct.(s).(n).(latentVar)(1:maxTrials);
        v2 = permStruct.(s).(n).([latentVar 'mat'])';
        all_corr = [all_corr {corr(v1, v2)}];
        
        v1_noise = normrnd(0,1,size(v1));
        v2_noise = normrnd(0,1,size(v2));
        noise_corr = [noise_corr {corr(v1_noise, v2_noise)}];
    end
end

real_median = cellfun(@(i) median(abs(i)), all_corr);
noise_median = cellfun(@(i) median(abs(i)), noise_corr);

real_CI = bootci(1e3, @median, real_median);
noise_CI = bootci(1e3, @median, noise_median);

fprintf('Real: median %0.3f (%0.3f - %0.3f)\n', median(real_median), real_CI(1), real_CI(2));
fprintf('Noise: median %0.3f (%0.3f - %0.3f)\n', median(noise_median), noise_CI(1), noise_CI(2));

% RPE neurons
% Real: median 0.114 (0.109 - 0.118)
% Noise: median 0.106 (0.102 - 0.109)

