%% Evolution Analysis
%This script analyses the evolution of firing activity over the course of
%the two blocks
%% Clear and load
% clearvars;
% close all;
load('RAWintBlocks.mat');
RAW=RAWblocks;

%% Generate spike rates per trial (Takes some time)


events = {'Cue';'PECue';'RD'};
event_intervals = [0. 0.5;
                  -0.5 0.5;
                   0.75 1.95];
base_intervals = [-11. -1.;
                  -11. -1.;
                  -11. -1.]; %all relative to cue
RD1 = strcmp('RD1',RAW(1).Einfo(:,2));
RD2 = strcmp('RD2',RAW(1).Einfo(:,2));
Cue = strcmp('Cue',RAW(1).Einfo(:,2));

num_trials = 60;
num_neurons = 0;
for session = 1:length(RAW)
    num_neurons = num_neurons + length(RAW(session).Nrast);
end
trial_info(num_neurons).blocks = [];

for event = 1:length(events) 
    counter = 1; 
    neuron_spike_rate = zeros(num_neurons,num_trials);
    baseline_spike_rate = zeros(num_neurons,num_trials);
    for session = 1:length(RAW)
        event_index = strcmp(events{event,1},RAW(session).Einfo(:,2));
        event_onsets = RAW(session).Erast{event_index, 1};
        baseline_onsets = RAW(session).Erast{Cue, 1};
        if event == 2 % PECue
            RD_onsets = RAW(session).Erast{strcmp('RD',RAW(session).Einfo(:,2)), 1};
            index = zeros(length(event_onsets),1);
            for i = 1:length(event_onsets)
                index(i) = any(abs(RD_onsets-event_onsets(i))<2);
            end
            event_onsets = event_onsets(logical(index));
        end
        intervals = event_onsets + event_intervals(event,:);
        baseline_intervals = baseline_onsets + base_intervals(event,:);
        
        for neuron = 1:length(RAW(session).Nrast)
            spike_times = RAW(session).Nrast{neuron,1};
            spike_count = NaN(num_trials,1);
            baseline_count = NaN(num_trials,1);
            for interval = 1:length(intervals)
                spike_count(interval,1) = sum(spike_times>=intervals(interval,1) & ...
                                              spike_times<=intervals(interval,2));

            end
            for interval = 1:length(baseline_intervals)
                baseline_count(interval,1) = sum(spike_times>=baseline_intervals(interval,1) & ...
                                 spike_times<=baseline_intervals(interval,2));
            end
            neuron_spike_rate(counter,:) = spike_count'./(event_intervals(event,2)-event_intervals(event,1));
            baseline_spike_rate(counter,:) = baseline_count'./(base_intervals(event,2)-base_intervals(event,1));
            name = RAW(session).Ninfo{neuron,1};
            if event == 1 % On the first pass make the trial info struct
                    if name(5) == 'n'
                        trial_info(counter).blocks = 0;
                    else
                        trial_info(counter).blocks = str2double(name(5));
                    end
                    trial_info(counter).session = session;
                    trial_info(counter).RD_suc = length(RAW(session).Erast{RD1});
                    trial_info(counter).RD_mal = length(RAW(session).Erast{RD2});
                    trial_info(counter).area = name(1:2);
                    reward_times = sortrows([RAW(session).Erast{RD1},ones(length(RAW(session).Erast{RD1}),1);
                                               RAW(session).Erast{RD2},zeros(length(RAW(session).Erast{RD2}),1)],1);
                    trial_info(counter).current_reward = reward_times(:,2);
                    % Find previous rewards for cues:
                    cue_times = RAW(session).Erast{Cue}(2:end);
                    cue_previous_reward = zeros(length(cue_times),1);
                    for cue = 1:length(cue_times)
                        [~,index] = max(reward_times(reward_times(:,1) < cue_times(cue),1)); %index of most recent reward consumed
                        if isempty(index)
                            cue_previous_reward(cue,1) = NaN; 
                        else
                            cue_previous_reward(cue,1) = reward_times(index,2); 
                        end
                    end
                    trial_info(counter).cue_previous_reward = cue_previous_reward;

            end
            counter = counter + 1;
        end
    end
    baseline_mean = nanmean(baseline_spike_rate,2);
    baseline_std = nanstd(baseline_spike_rate,0,2);
    events{event, 2} = (neuron_spike_rate-baseline_mean)./baseline_std; %z-score
    events{event, 3} = (neuron_spike_rate-baseline_mean); %just baseline substracted for linear models
    events{event, 4} = (neuron_spike_rate-baseline_mean)./(max(neuron_spike_rate,[],2)-min(neuron_spike_rate,[],2)); %baseline subtracted and normalized to total range

end
%% Define responsive logical vectors

interspersed = [trial_info.blocks]' == 0;
blocks1 = [trial_info.blocks]' == 1;
blocks2 = [trial_info.blocks]' == 2;
area_vectors = [strcmp({trial_info.area},'NA')',strcmp({trial_info.area},'VP')'];
blocks_vectors = [[trial_info.blocks]'==0,[trial_info.blocks]'>0]; %Interspersed, blocks
session_type = {interspersed,blocks1,blocks2};
session_titles = {'Interspersed','Sucrose First','Malt First'};


%% five samplings of 3 trials, equally space through each block of rewards

fives = NaN(num_neurons,10);
for event = 1:size(events,1)
    for neuron = 1:num_neurons
        trial_firing = events{event,2}(neuron,:); % Pull up the trial by trial firing
        trial_firing = trial_firing(~isnan(trial_firing)); % Remove any NaN values
        
        first1 = 1;
        last2 = length(trial_firing);
        switch trial_info(neuron).blocks
            case 0 % Interspersed
                suc_trials = find(trial_info(neuron).current_reward == 1);
                malt_trials = find(trial_info(neuron).current_reward == 0);
                suc_first = suc_trials(1:3);
                suc_quarter = ceil(length(suc_trials)/4);
                suc_second = suc_trials(suc_quarter-1:suc_quarter+1);
                suc_center = ceil(length(suc_trials)/2);
                suc_mid = suc_trials(suc_center-1:suc_center+1);
                suc_3rdquarter = ceil(3*length(suc_trials)/4);
                suc_fourth = suc_trials(suc_3rdquarter-1:suc_3rdquarter+1);
                suc_last = suc_trials(end-2:end);  
                
                malt_first = malt_trials(1:3);
                malt_quarter = ceil(length(malt_trials)/4);
                malt_second = malt_trials(malt_quarter-1:malt_quarter+1);
                malt_center = ceil(length(malt_trials)/2);
                malt_mid = malt_trials(malt_center-1:malt_center+1);
                malt_3rdquarter = ceil(3*length(malt_trials)/4);
                malt_fourth = malt_trials(malt_3rdquarter-1:malt_3rdquarter+1);
                malt_last = malt_trials(end-2:end);  
                
                slices = [suc_first';malt_first';suc_second';malt_second';suc_mid';malt_mid';suc_fourth';malt_fourth';suc_last';malt_last'];
            case 1 % Suc first
                first2 = trial_info(neuron).RD_suc+1;
                last1 = first2-1;
                mid1 = ceil((last1-first1)/2)+first1;
                mid2 = ceil((last2-first2)/2)+first2;
                second1 = ceil((last1-first1)/4)+first1;
                second2 = ceil((last2-first2)/4)+first2;
                fourth1 = ceil(3*(last1-first1)/4)+first1;
                fourth2 = ceil(3*(last2-first2)/4)+first2;
                slices = [first1:first1+2;...
                          second1-1:second1+1;...
                          mid1-1:mid1+1;...
                          fourth1-1:fourth1+1;...
                          last1-2:last1;...
                          first2:first2+2;...
                          second2-1:second2+1;...
                          mid2-1:mid2+1;...
                          fourth2-1:fourth2+1;...
                          last2-2:last2];
            case 2 % Malt first
                first2 = trial_info(neuron).RD_mal+1;
                last1 = first2-1;
                mid1 = ceil((last1-first1)/2)+first1;
                mid2 = ceil((last2-first2)/2)+first2;
                second1 = ceil((last1-first1)/4)+first1;
                second2 = ceil((last2-first2)/4)+first2;
                fourth1 = ceil(3*(last1-first1)/4)+first1;
                fourth2 = ceil(3*(last2-first2)/4)+first2;
                slices = [first1:first1+2;...
                          second1-1:second1+1;...
                          mid1-1:mid1+1;...
                          fourth1-1:fourth1+1;...
                          last1-2:last1;...
                          first2:first2+2;...
                          second2-1:second2+1;...
                          mid2-1:mid2+1;...
                          fourth2-1:fourth2+1;...
                          last2-2:last2];
        end

        for slice = 1:size(slices,1)
            fives(neuron,slice) = mean(trial_firing(slices(slice,:)));
        end
    end
    events{event,5} = fives; %Just finished changing fives here
end

save('evolution_data.mat','events')
save('evolution_info.mat','trial_info')
disp('done')