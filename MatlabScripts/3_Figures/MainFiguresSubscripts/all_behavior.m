function BehavStats=all_behavior(RAWblocks, RAWCued, plot_region, plot_sub_figures, perform_stats, Colors)

%% interspersed vs blocks

%subselect sessions of interest
if strcmp(plot_region,'Both')
    RAW=RAWblocks;
else
    
    for i=1:length(RAWblocks)
        name=char(RAWblocks(i).Region);
        region{i,1}=name(1:2);
    end

    region_sessions = strcmp(plot_region,region);
    RAW=RAWblocks(region_sessions);
end

% Latency
% Find the indicies of all blocks and interspersed sessions
blocked_sessions = strcmp('B',{RAW.Blocks})';
interspersed_sessions = strcmp('I',{RAW.Blocks})';
num_sessions = length(RAW);
num_trials = 60;

% preallocate arrays
latencyR1 = nan(num_sessions, num_trials);
latencyR2 = nan(num_sessions, num_trials);
fraction_completed = zeros(num_sessions,2); %suc,malt
rat_name = cell(num_sessions,1);

for session = 1:num_sessions
    
    % Make logical index vectors to call every event type
    PER1 = strcmp('PER1',RAW(session).Einfo(:,2)); %finds the index of the event
    PER2 = strcmp('PER2',RAW(session).Einfo(:,2)); %finds the index of the event
    CueR1 = strcmp('CueR1',RAW(session).Einfo(:,2)); %finds the index of the event
    CueR2 = strcmp('CueR2',RAW(session).Einfo(:,2)); %finds the index of the event
    PECue = strcmp('PECue',RAW(session).Einfo(:,2)); %finds the index of the event
    
    % Latency Stuff
    suc_latencies = cell2mat(MakePSR04(RAW(session).Erast(PECue),RAW(session).Erast{CueR1},[0 20],{2,'first'}));
    malt_latencies = cell2mat(MakePSR04(RAW(session).Erast(PECue),RAW(session).Erast{CueR2},[0 20],{2,'first'}));
    latencyR1(session,1:length(suc_latencies)) = suc_latencies;
    latencyR2(session,1:length(malt_latencies)) = malt_latencies;
    
    % Trial Completion Stuff
    suc_completed = length(RAW(session).Erast{PER1,1});
    suc_total = length(RAW(session).Erast{CueR1,1});
    malt_completed = length(RAW(session).Erast{PER2,1});
    malt_total = length(RAW(session).Erast{CueR2,1});
    fraction_completed(session,:) = [suc_completed/suc_total,malt_completed/malt_total];
    
    % Rat names
    name=char(RAW(session).Ninfo(1,1));
    rat_name{session,1} = name(1:3); 
end

% Latency stuff
median_latencyR1 = median(latencyR1,2,'omitnan');
median_latencyR2 = median(latencyR2,2,'omitnan');

% Average for rat
unique_rat_names = unique(rat_name);
rat_latency_data = zeros(length(unique_rat_names),4);
rat_completion_data = zeros(length(unique_rat_names),4);
rat_latency_data_sem = zeros(length(unique_rat_names),4);
rat_completion_data_sem = zeros(length(unique_rat_names),4);
for rat = 1:length(unique_rat_names)
    rat_index = strcmp(unique_rat_names{rat}, rat_name);
    rat_latency_data(rat,:) = [mean(median_latencyR1(interspersed_sessions & rat_index)),...  % Sucrose inter
                               mean(median_latencyR2(interspersed_sessions & rat_index)),...  % Malto inter
                               mean(median_latencyR1(blocked_sessions & rat_index)),...  % Sucrose block
                               mean(median_latencyR2(blocked_sessions & rat_index))];  % Malto block
    rat_completion_data(rat,:) = [mean(fraction_completed(interspersed_sessions & rat_index, 1)),...  % Sucrose inter
                                  mean(fraction_completed(interspersed_sessions & rat_index, 2)),...  % Malto inter
                                  mean(fraction_completed(blocked_sessions & rat_index, 1)),...  % Sucrose block
                                  mean(fraction_completed(blocked_sessions & rat_index, 2))];  % Malto block    
    rat_latency_data_sem(rat,:) = [nanste(median_latencyR1(interspersed_sessions & rat_index),1),...  % Sucrose inter
                                   nanste(median_latencyR2(interspersed_sessions & rat_index),1),...  % Malto inter
                                   nanste(median_latencyR1(blocked_sessions & rat_index),1),...  % Sucrose block
                                   nanste(median_latencyR2(blocked_sessions & rat_index),1)];  % Malto block
    rat_completion_data_sem(rat,:) = [nanste(fraction_completed(interspersed_sessions & rat_index, 1),1),...  % Sucrose inter
                                      nanste(fraction_completed(interspersed_sessions & rat_index, 2),1),...  % Malto inter
                                      nanste(fraction_completed(blocked_sessions & rat_index, 1),1),...  % Sucrose block
                                      nanste(fraction_completed(blocked_sessions & rat_index, 2),1)];  % Malto block
end

mean_latency_data = [mean(median_latencyR1(interspersed_sessions)),...  % Sucrose inter
                     mean(median_latencyR2(interspersed_sessions)),...  % Malto inter
                     mean(median_latencyR1(blocked_sessions)),...  % Sucrose block
                     mean(median_latencyR2(blocked_sessions));...  % Malto block
                     nanste(median_latencyR1(interspersed_sessions),1),...  % Sucrose inter
                     nanste(median_latencyR2(interspersed_sessions),1),...  % Malto inter
                     nanste(median_latencyR1(blocked_sessions),1),...  % Sucrose block
                     nanste(median_latencyR2(blocked_sessions),1)];  % Malto block
mean_completion_data = [mean(fraction_completed(interspersed_sessions, 1)),...  % Sucrose inter
                        mean(fraction_completed(interspersed_sessions, 2)),...  % Malto inter
                        mean(fraction_completed(blocked_sessions, 1)),...  % Sucrose block
                        mean(fraction_completed(blocked_sessions, 2));...  % Malto block
                        nanste(fraction_completed(interspersed_sessions, 1),1),...  % Sucrose inter
                        nanste(fraction_completed(interspersed_sessions, 2),1),...  % Malto inter
                        nanste(fraction_completed(blocked_sessions, 1),1),...  % Sucrose block
                        nanste(fraction_completed(blocked_sessions, 2),1)];  % Malto block

% make figure
figure('Units','inches','Position',[3,3,8,5]);
% Latency subplot
mean_data = {mean_latency_data, mean_completion_data};
rat_data = {rat_latency_data, rat_completion_data, rat_latency_data_sem, rat_completion_data_sem};
color = {Colors('sucrose'),Colors('maltodextrin'),Colors('blockrose'),Colors('blockrodextrin')};
titles = {'Latency','Completion'};

for plot_type = 1:2
    subplot(3,4,[5+(plot_type-1)*2 6+(plot_type-1)*2])
    hold on
    for rat = 1:length(rat_latency_data)
        xvals=[0.9 1.9]+(rat-1)*0.2/(length(unique_rat_names)-1);
        errorbar(xvals,rat_data{plot_type}(rat,1:2),rat_data{plot_type+2}(rat,1:2),'Color',[.7 .7 .7])
        errorbar(xvals+2,rat_data{plot_type}(rat,3:4),rat_data{plot_type+2}(rat,3:4),'Color',[.7 .7 .7])
    end
    for i = 1:4
        errorbar(i,mean_data{plot_type}(1,i),mean_data{plot_type}(2,i),'o','Color',color{i},'MarkerFaceColor',color{i})
    end

    if plot_type == 1
%         plot([1,2],[5,5],'k')
%         plot([3,4],[4.5,4.5],'k')
%         plot([1.5,3.5],[5.5,5.5],'k')
%         plot([1.5,1.5],[5.5,5],'k')
%         plot([3.5,3.5],[5.5,4.5],'k')
%         text(2.5, 5.7, '*','Color','k','FontSize',16)
%         text(2.98, 1.77, '#','Color','k') %This is to note that the two sucrose latencies are significantly different, which should be explained in the figure legend.
        axis([0.5 9.5 0 10]);
        ylabel('Latency (sec)');
    else
        axis([0.5 9.5 0 1]);
        ylabel('Fraction of trials completed')
    end

end

subplot(3,4,1:3)
sessions_to_use = [1 8];
% figure;

session_end = strcmp('MedEnd',RAW(sessions_to_use(1)).Einfo(:,2)); %finds the index of the event
stop(1) = RAW(sessions_to_use(1)).Erast{session_end,1};

session_end = strcmp('MedEnd',RAW(sessions_to_use(2)).Einfo(:,2)); %finds the index of the event
stop(2) = RAW(sessions_to_use(2)).Erast{session_end,1};
line([0,1+stop(1)/60],[4.3,4.3],'Color','k')
line([0,1+stop(2)/60],[6.3,6.3],'Color','k')
for session = 1:length(sessions_to_use)
    
    %get events for this session
    RD1 = strcmp('RD1',RAW(sessions_to_use(session)).Einfo(:,2)); %finds the index of the event
    RD2 = strcmp('RD2',RAW(sessions_to_use(session)).Einfo(:,2)); %finds the index of the event
    
    suc_delivery = RAW(sessions_to_use(session)).Erast{RD1,1}/60;
    malt_delivery = RAW(sessions_to_use(session)).Erast{RD2,1}/60;

    for trial = 1:length(suc_delivery)
        line([suc_delivery(trial),suc_delivery(trial)],[3+2*(session-1)+1.1,3+2*(session-1)+1.7],'Color',color{3-2*(session-1)})
    end
    for trial = 1:length(malt_delivery)
        line([malt_delivery(trial),malt_delivery(trial)],[3+2*(session-1)+0.9,3+2*(session-1)+1.5],'Color',color{4-2*(session-1)})
    end
end

% title('Two Example Sessions','FontSize',11)

if perform_stats 
    % Latency
    suc_interspersed = median_latencyR1(interspersed_sessions);
    malt_interspersed = median_latencyR2(interspersed_sessions);
    suc_blocked = median_latencyR1(blocked_sessions);
    malt_blocked = median_latencyR2(blocked_sessions);
    
    num_interspersed = length(suc_interspersed);
    num_blocked = length(suc_blocked);
    rat_interspersed = rat_name(interspersed_sessions);
    rat_blocked = rat_name(blocked_sessions);
    
%     t=table(cat(1,rat_interspersed,rat_blocked),categorical(cat(1,zeros(num_interspersed,1),ones(num_blocked,1))),cat(1,suc_interspersed,suc_blocked),cat(1,malt_interspersed,malt_blocked),'variablenames',{'subject','blocked','sucrose','maltodextrin'});
%     vars=table({'Suc'; 'Mal'},'variablenames',{'Reward'});
%     rm=fitrm(t,'sucrose-maltodextrin~subject+blocked','WithinDesign',vars);
%     BehavStats.IntblockLatency.ranovatbl=ranova(rm);
%     BehavStats.IntblockLatency.blocked=multcompare(rm,'blocked');
%     BehavStats.IntblockLatency.reward=multcompare(rm,'Reward');
%     BehavStats.IntblockLatency.blockedbyreward=multcompare(rm,'blocked','by','Reward');
    
    
%     figure;
    [~,BehavStats.IntblockLatency.anovatbl,stats] = anovan(cat(1,suc_interspersed,malt_interspersed,suc_blocked,malt_blocked),...
        {cat(1,zeros(num_interspersed*2,1),ones(num_blocked*2,1)),...
        cat(1,zeros(num_interspersed,1),ones(num_interspersed,1),zeros(num_blocked,1),ones(num_blocked,1)),...
        cat(1,rat_interspersed,rat_interspersed,rat_blocked,rat_blocked)},...
        'varnames',{'Int/Block','Reward','Subject'},'model',[1 0 0;1 1 0;0 0 1],... %[1 0 0;0 1 0;1 1 0;0 0 1],...
        'random',3,'display','off');
    [BehavStats.IntblockLatency.multcompare] = multcompare(stats,'dimension',[1 2],'display','off');
    BehavStats.IntblockLatency.stats=stats;
    
    % Trials completed
    suc_inter = fraction_completed(interspersed_sessions,1);
    malt_inter = fraction_completed(interspersed_sessions,2);
    suc_block = fraction_completed(blocked_sessions,1);
    malt_block = fraction_completed(blocked_sessions,2);
    
%     t=table(cat(1,rat_interspersed,rat_blocked),categorical(cat(1,zeros(num_interspersed,1),ones(num_blocked,1))),cat(1,suc_inter,suc_block),cat(1,malt_inter,malt_block),'variablenames',{'subject','blocked','sucrose','maltodextrin'});
%     vars=table({'Suc'; 'Mal'},'variablenames',{'Reward'});
%     rm=fitrm(t,'sucrose-maltodextrin~subject+blocked','WithinDesign',vars);
%     BehavStats.IntblockCompletion.ranovatbl=ranova(rm);
%     BehavStats.IntblockCompletion.blocked=multcompare(rm,'blocked');
%     BehavStats.IntblockCompletion.reward=multcompare(rm,'Reward');
%     BehavStats.IntblockCompletion.blockedbyreward=multcompare(rm,'blocked','by','Reward');
    
%     figure;
    [~,BehavStats.IntblockCompletion.anovatbl,stats] = anovan(cat(1,suc_inter,malt_inter,suc_block,malt_block),...
        {cat(1,zeros(num_interspersed*2,1),ones(num_blocked*2,1)),...
        cat(1,zeros(num_interspersed,1),ones(num_interspersed,1),zeros(num_blocked,1),ones(num_blocked,1)),...
        cat(1,rat_interspersed,rat_interspersed,rat_blocked,rat_blocked)},...
        'varnames',{'Int/Block','Reward','Subject'},'model',[1 0 0;1 1 0;0 0 1],... %[1 0 0;0 1 0;1 1 0;0 0 1],...
        'random',3,'display','off');
    [BehavStats.IntblockCompletion.multcompare] = multcompare(stats,'dimension',[1 2],'display','off');
    BehavStats.IntblockCompletion.stats=stats;
end

%lick plots

global Dura Tm BSIZE %Tbin
BSIZE=0.05; %what's the bin for the PSTH?
Dura=[-4 15]; Tm=Dura(1):BSIZE:Dura(2);
%Tbin=-0.5:0.005:0.5; %window used to determine the optimal binsize

%Smoothing PSTHs
%smoothing filter
smoothbins=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',8); %std=3.98
filterweights=pdf(halfnormal,[0:smoothbins]);

%generate lick PSTHs
for session = 1:num_sessions
    
    % Make logical index vectors to call every event type
    R1 = strcmp('RD1',RAW(session).Einfo(:,2)); %finds the index of the event
    R2 = strcmp('RD2',RAW(session).Einfo(:,2)); %finds the index of the event
    Licks = strcmp('Licks',RAW(session).Einfo(:,2)); %finds the index of the event
    
    EvList{1}=RAW(session).Erast{R1};
    EvList{2}=RAW(session).Erast{R2};

    for k=1:2
        [PSR1,N1]=MakePSR04(RAW(session).Erast(Licks),EvList{k},Dura,{1});% makes collpased rasters. PSR1 is a cell(neurons)
        [PTH1,~,~]=MakePTH07(PSR1,repmat(N1,size(RAW(session).Erast(Licks),1),1),{2,0,BSIZE});%-----Fixed bin used here    

        PTH1smooth=[];
        for l=1:length(Tm)
           PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 smoothbins]):l).*fliplr(filterweights(1:min([l smoothbins+1]))))/sum(filterweights(1:min([l smoothbins+1])));                  
        end
        Lick.IntBlock(k).PSTHraw(session,:)=PTH1smooth;

    end
end

%plotting
Xaxis=[-2 12];
Ishow=find(Tm>=Xaxis(1) & Tm<=Xaxis(2));
time1=Tm(Ishow);

for type=1:2
    
    if type==1 Sel=interspersed_sessions; end
    if type==2 Sel=blocked_sessions; end
    
    psth1=nanmean(Lick.IntBlock(1).PSTHraw(Sel,Ishow),1);
    sem1=nanste(Lick.IntBlock(1).PSTHraw(Sel,Ishow),1); %calculate standard error of the mean
    up1=psth1+sem1;
    down1=psth1-sem1;

    psthE=nanmean(Lick.IntBlock(2).PSTHraw(Sel,Ishow),1);
    semE=nanste(Lick.IntBlock(2).PSTHraw(Sel,Ishow),1); %calculate standard error of the mean
    upE=psthE+semE;
    downE=psthE-semE;

    %plotting
    subplot(4,4,12+type);
    hold on;
    plot(time1,psth1,'Color',color{1+(type-1)*2},'linewidth',1);
    plot(time1,psthE,'Color',color{2+(type-1)*2},'linewidth',1);
    patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],color{1+(type-1)*2},'EdgeColor','none');alpha(0.5);
    patch([time1,time1(end:-1:1)],[upE,downE(end:-1:1)],color{2+(type-1)*2},'EdgeColor','none');alpha(0.5);
    plot([0 0],[-2 8],':','color','k','linewidth',0.75);
    axis([Xaxis(1) Xaxis(2) 0 8]);
    xlabel('Seconds from reward delivery');
    if type==1 ylabel('Mean lick rate (licks/s)'); end
    if type==1 title('Interspersed'); end
    if type==2 title('Blocked'); end
    legend('Suc trials','Mal trials');

end

%% cued sessions

%subselect sessions of interest
for i=1:length(RAWCued)
    name=char(RAWCued(i).Ninfo(1,1));
    sessrat(i,1)=str2num(name(2:3));
    sessday(i,1)=str2num(name(8:9));
end

%decide which sessions to include
inc_sessions = (sessrat==3 | sessrat==4 | sessrat==9 | sessrat==10) & sessday>10;
RAW=RAWCued(inc_sessions);
%RAW=RAWCued; %all sessions

% Latency
% Find the indicies of all blocks and interspersed sessions
num_sessions = length(RAW);
num_trials = 60; %just for pre-allocating

% preallocate arrays
latencyNo1 = nan(num_sessions, num_trials); %non-predictive sucrose
latencyNo2 = nan(num_sessions, num_trials); %non-predictive maltodextrin
latencyPr1 = nan(num_sessions, num_trials); %predictive sucrose
latencyPr2 = nan(num_sessions, num_trials); %predictive maltodextrin
fraction_completed = zeros(num_sessions,4); %Predictive suc, pm, ns, nm
rat_name = cell(num_sessions,1);

for session = 1:num_sessions
    
    % Make logical index vectors to call every event type
    PENo1 = strcmp('PECue31',RAW(session).Einfo(:,2)); %finds the index of the event
    PENo2 = strcmp('PECue32',RAW(session).Einfo(:,2)); %finds the index of the event  
    PEPr1 = strcmp('PECue1',RAW(session).Einfo(:,2)); %finds the index of the event
    PEPr2 = strcmp('PECue2',RAW(session).Einfo(:,2)); %finds the index of the event
    Cue1 = strcmp('Cue1',RAW(session).Einfo(:,2)); %finds the index of the event
    Cue2 = strcmp('Cue2',RAW(session).Einfo(:,2)); %finds the index of the event
    Cue31 = strcmp('Cue31',RAW(session).Einfo(:,2)); %finds the index of the event
    Cue32 = strcmp('Cue32',RAW(session).Einfo(:,2)); %finds the index of the event
    PECue = strcmp('PECue',RAW(session).Einfo(:,2)); %finds the index of the event
    
    % Latency Stuff
    No1_latencies = cell2mat(MakePSR04(RAW(session).Erast(PECue),RAW(session).Erast{Cue31},[0 20],{2,'first'}));
    latencyNo1(session,1:length(No1_latencies)) = No1_latencies;

    No2_latencies = cell2mat(MakePSR04(RAW(session).Erast(PECue),RAW(session).Erast{Cue32},[0 20],{2,'first'}));
    latencyNo2(session,1:length(No2_latencies)) = No2_latencies;
    
    Pr1_latencies = cell2mat(MakePSR04(RAW(session).Erast(PECue),RAW(session).Erast{Cue1},[0 20],{2,'first'}));
    latencyPr1(session,1:length(Pr1_latencies)) = Pr1_latencies;

    Pr2_latencies = cell2mat(MakePSR04(RAW(session).Erast(PECue),RAW(session).Erast{Cue2},[0 20],{2,'first'}));
    latencyPr2(session,1:length(Pr2_latencies)) = Pr2_latencies;


    
    % Trial Completion Stuff
    No1_completed = length(RAW(session).Erast{PENo1,1});
    No1_total = length(RAW(session).Erast{Cue31,1});
    No2_completed = length(RAW(session).Erast{PENo2,1});
    No2_total = length(RAW(session).Erast{Cue32,1});
    Pr1_completed = length(RAW(session).Erast{PEPr1,1});
    Pr1_total = length(RAW(session).Erast{Cue1,1});
    Pr2_completed = length(RAW(session).Erast{PEPr2,1});
    Pr2_total = length(RAW(session).Erast{Cue2,1});
    fraction_completed(session,:) = [No1_completed/No1_total,No2_completed/No2_total,Pr1_completed/Pr1_total,Pr2_completed/Pr2_total];
    
    % Rat names
    name=char(RAW(session).Ninfo(1,1));
    rat_name{session,1} = name(1:3); 
end

% Latency stuff
median_latencyNo1 = median(latencyNo1,2,'omitnan');
median_latencyNo2 = median(latencyNo2,2,'omitnan');
median_latencyPr1 = median(latencyPr1,2,'omitnan');
median_latencyPr2 = median(latencyPr2,2,'omitnan');


% Average for rat
unique_rat_names = unique(rat_name);
rat_latency_data = zeros(length(unique_rat_names),4);
rat_completion_data = zeros(length(unique_rat_names),4);
rat_latency_data_sem = zeros(length(unique_rat_names),4);
rat_completion_data_sem = zeros(length(unique_rat_names),4);
for rat = 1:length(unique_rat_names)
    rat_index = strcmp(unique_rat_names{rat}, rat_name);
    rat_latency_data(rat,:) = [mean(median_latencyNo1(rat_index)),...  % Sucrose non-pr
                               mean(median_latencyNo2(rat_index)),...  % Malto non-pr
                               mean(median_latencyPr1(rat_index)),...  % Sucrose pr
                               mean(median_latencyPr2(rat_index))];  % Malto pr
    rat_completion_data(rat,:) = [mean(fraction_completed(rat_index, 1)),...  % Sucrose non-pr
                                  mean(fraction_completed(rat_index, 2)),...  % Malto non-pr
                                  mean(fraction_completed(rat_index, 3)),...  % Sucrose pr
                                  mean(fraction_completed(rat_index, 4))];  % Malto pr    
    rat_latency_data_sem(rat,:) = [nanste(median_latencyNo1(rat_index),1),...  % Sucrose non-pr
                                   nanste(median_latencyNo2(rat_index),1),...  % Malto non-pr
                                   nanste(median_latencyPr1(rat_index),1),...  % Sucrose pr
                                   nanste(median_latencyPr2(rat_index),1)];  % Malto pr
    rat_completion_data_sem(rat,:) = [nanste(fraction_completed(rat_index, 1),1),...  
                                      nanste(fraction_completed(rat_index, 2),1),...  
                                      nanste(fraction_completed(rat_index, 3),1),...  
                                      nanste(fraction_completed(rat_index, 4),1)];  
end

mean_latency_data = [mean(median_latencyNo1),...  % Sucrose non-pr
                     mean(median_latencyNo2),...  % Malto non-pr
                     mean(median_latencyPr1),...  % Sucrose pr
                     mean(median_latencyPr2);...  % Malto pr
                     nanste(median_latencyNo1,1),...  
                     nanste(median_latencyNo2,1),...  
                     nanste(median_latencyPr1,1),...  
                     nanste(median_latencyPr2,1)];  
mean_completion_data = [mean(fraction_completed(:, 1)),... 
                        mean(fraction_completed(:, 2)),...  
                        mean(fraction_completed(:, 3)),...  
                        mean(fraction_completed(:, 4));...  
                        nanste(fraction_completed(:, 1),1),...  
                        nanste(fraction_completed(:, 2),1),...  
                        nanste(fraction_completed(:, 3),1),...  
                        nanste(fraction_completed(:, 4),1)];  


% Latency subplot
mean_data = {mean_latency_data, mean_completion_data};
rat_data = {rat_latency_data, rat_completion_data, rat_latency_data_sem, rat_completion_data_sem};
color = {Colors('sucrose'),Colors('maltodextrin'),Colors('blockrose'),Colors('blockrodextrin')};
titles = {'Latency','Completion'};

for plot_type = 1:2
    subplot(3,4,[5+(plot_type-1)*2 6+(plot_type-1)*2])
    hold on
    for rat = 1:length(rat_latency_data)
        xvals=[5.9 6.9]+(rat-1)*0.2/(length(unique_rat_names)-1);
        errorbar(xvals,rat_data{plot_type}(rat,1:2),rat_data{plot_type+2}(rat,1:2),'Color',[.7 .7 .7])
        errorbar(xvals+2,rat_data{plot_type}(rat,3:4),rat_data{plot_type+2}(rat,3:4),'Color',[.7 .7 .7])
    end
    for i = 1:4
        errorbar(i+5,mean_data{plot_type}(1,i),mean_data{plot_type}(2,i),'o','Color',color{i},'MarkerFaceColor',color{i})
    end

    if plot_type == 1
%         plot([6,7],[5,5],'k')
%         plot([8,9],[4.5,4.5],'k')
%         plot([6.5,8.5],[5.5,5.5],'k')
%         plot([6.5,6.5],[5.5,5],'k')
%         plot([8.5,8.5],[5.5,4.5],'k')
%         text(7.5, 5.7, '*','Color','k','FontSize',16)
%         text(2.98, 1.77, '#','Color','k') %This is to note that the two sucrose latencies are significantly different, which should be explained in the figure legend.
        axis([0.5 9.5 0 10]);
        ylabel('Port entry latency (sec)');
        text(0.5,11,'Interspersed');
        text(3,11,'Blocked');
        text(5.5,11,'Non-predictive');
        text(8,11,'Predictive');
    else
        axis([0.5 9.5 0 1]);
        ylabel('Fraction of trials completed')
        text(0.5,1.1,'Interspersed');
        text(3,1.1,'Blocked');
        text(5.5,1.1,'Non-predictive');
        text(8,1.1,'Predictive');
    end
    conditions = {"Suc","Mal","Suc","Mal","Suc","Mal","Suc","Mal"};
    set(gca,'XTick',[1:4 6:9],'XTickLabel',conditions)
    
end

subplot(3,4,1:3)
session_to_use = 2;

session_end = strcmp('MedEnd',RAW(session_to_use(1)).Einfo(:,2)); %finds the index of the event
stop = RAW(session_to_use(1)).Erast{session_end,1};

line([0,1+stop(1)/60],[2.3,2.3],'Color','k')
    
%get events for this session
Cue1RD = strcmp('Cue1RD',RAW(session_to_use).Einfo(:,2)); %finds the index of the event
Cue2RD = strcmp('Cue2RD',RAW(session_to_use).Einfo(:,2)); %finds the index of the event
Cue3RD1 = strcmp('Cue3RD1',RAW(session_to_use).Einfo(:,2)); %finds the index of the event
Cue3RD2 = strcmp('Cue3RD2',RAW(session_to_use).Einfo(:,2)); %finds the index of the event

suc_deliveryN = RAW(session_to_use).Erast{Cue3RD1,1}/60;
malt_deliveryN = RAW(session_to_use).Erast{Cue3RD2,1}/60;
suc_deliveryP = RAW(session_to_use).Erast{Cue1RD,1}/60;
malt_deliveryP = RAW(session_to_use).Erast{Cue2RD,1}/60;

for trial = 1:length(suc_deliveryN)
    a=line([suc_deliveryN(trial),suc_deliveryN(trial)],[2.1,2.7],'Color',color{1});
end
for trial = 1:length(malt_deliveryN)
    b=line([malt_deliveryN(trial),malt_deliveryN(trial)],[1.9,2.5],'Color',color{2});
end
for trial = 1:length(suc_deliveryP)
    c=line([suc_deliveryP(trial),suc_deliveryP(trial)],[2.1,2.7],'Color',color{3});
end
for trial = 1:length(malt_deliveryP)
    d=line([malt_deliveryP(trial),malt_deliveryP(trial)],[1.9,2.5],'Color',color{4});
end

text(60,6.3,'Interspersed');
text(60,4.3,'Blocked');
text(105,2.3,'Cued');


text(0,3,'Sucrose, unpredictable','Color',color{1});
text(30,3,'Sucrose, predicted by cue','Color',color{3});
text(0,1.6,'Maltodextrin, unpredictable','Color',color{2});
text(30,1.6,'Maltodextrin, predicted by cue','Color',color{4});


XTick=[0 20 40 60 80 100];
axis([0 110 0.9 7]);
ax1 = gca;   
ax1.YAxis.Visible = 'off'; 
xlabel('Session time (minutes)')



if perform_stats
   
    % Latency
    
%     t=table(rat_name,median_latencyNo1,median_latencyNo2,median_latencyPr1,median_latencyPr2,'variablenames',{'subject','NPSuc','NPMal','PrSuc','PrMal'});
%     vars=table({'NPSuc';'NPMal';'PrSuc';'PrMal'},{'Non'; 'Non'; 'Pr' ;'Pr'},{'Suc'; 'Mal' ;'Suc' ;'Mal'},'variablenames',{'TrialType','Predictive','Reward'});
%     rm=fitrm(t,'NPSuc-PrMal~subject','WithinDesign',vars);
%     BehavStats.CuedLatency.ranovatbl=ranova(rm);
%     BehavStats.CuedLatency.predictive=multcompare(rm,'Predictive');
%     BehavStats.CuedLatency.reward=multcompare(rm,'Reward');
%     BehavStats.CuedLatency.predictivebyreward=multcompare(rm,'Predictive','by','Reward');
    
    suc_interspersed = median_latencyNo1;
    malt_interspersed = median_latencyNo2;
    suc_blocked = median_latencyPr1;
    malt_blocked = median_latencyPr2;

    [~,BehavStats.CuedLatency.anovatbl,stats] = anovan(cat(1,suc_interspersed,malt_interspersed,suc_blocked,malt_blocked),...
        {cat(1,zeros(num_sessions*2,1),ones(num_sessions*2,1)),...
        cat(1,zeros(num_sessions,1),ones(num_sessions,1),zeros(num_sessions,1),ones(num_sessions,1)),...
        cat(1,rat_name,rat_name,rat_name,rat_name)},...
        'varnames',{'Non-pr/Pr','Reward','Subject'},'model',[1 0 0;1 1 0;0 0 1],... %[1 0 0;0 1 0;1 1 0;0 0 1],...
        'random',3,'display','off');
    [BehavStats.CuedLatency.multcompare,BehavStats.CuedLatency.multcompare_means] = multcompare(stats,'dimension',[1],'display','off');
    BehavStats.CuedLatency.stats=stats;

    

    % Trials completed
    
%     t=table(rat_name,fraction_completed(:,1),fraction_completed(:,2),fraction_completed(:,3),fraction_completed(:,4),'variablenames',{'subject','NPSuc','NPMal','PrSuc','PrMal'});
%     vars=table({'NPSuc';'NPMal';'PrSuc';'PrMal'},{'Non'; 'Non'; 'Pr' ;'Pr'},{'Suc'; 'Mal' ;'Suc' ;'Mal'},'variablenames',{'TrialType','Predictive','Reward'});
%     rm=fitrm(t,'NPSuc-PrMal~subject','WithinDesign',vars);
%     BehavStats.CuedCompletion.ranovatbl=ranova(rm);
%     BehavStats.CuedCompletion.predictive=multcompare(rm,'Predictive');
%     BehavStats.CuedCompletion.reward=multcompare(rm,'Reward');
%     BehavStats.CuedCompletion.predictivebyreward=multcompare(rm,'Predictive','by','Reward');   
    
    
    suc_inter = fraction_completed(:,1);
    malt_inter = fraction_completed(:,2);
    suc_block = fraction_completed(:,3);
    malt_block = fraction_completed(:,4);
%     
%     figure;
    [~,BehavStats.CuedCompletion.anovatbl,stats] = anovan(cat(1,suc_inter,malt_inter,suc_block,malt_block),...
        {cat(1,zeros(num_sessions*2,1),ones(num_sessions*2,1)),...
        cat(1,zeros(num_sessions,1),ones(num_sessions,1),zeros(num_sessions,1),ones(num_sessions,1)),...
        cat(1,rat_name,rat_name,rat_name,rat_name)},...
        'varnames',{'Non-pr/Pr','Reward','Subject'},'model',[1 0 0;1 1 0;0 0 1],... %[1 0 0;0 1 0;1 1 0;0 0 1],...
        'random',3,'display','off');
    [BehavStats.CuedCompletion.multcompare] = multcompare(stats,'dimension',[1 2],'display','off');
    BehavStats.CuedCompletion.stats=stats;
end


%lick plots

%generate lick PSTHs
for session = 1:num_sessions
    
    % Make logical index vectors to call every event type
    R1 = strcmp('Cue1RD',RAW(session).Einfo(:,2)); %finds the index of the event
    R2 = strcmp('Cue2RD',RAW(session).Einfo(:,2)); %finds the index of the event
    R31 = strcmp('Cue3RD1',RAW(session).Einfo(:,2)); %finds the index of the event
    R32 = strcmp('Cue3RD2',RAW(session).Einfo(:,2)); %finds the index of the event
    Licks = strcmp('Licks',RAW(session).Einfo(:,2)); %finds the index of the event
    
    EvList{1}=RAW(session).Erast{R1};
    EvList{2}=RAW(session).Erast{R2};
    EvList{3}=RAW(session).Erast{R31};
    EvList{4}=RAW(session).Erast{R32};
    
    for k=1:4
        [PSR1,N1]=MakePSR04(RAW(session).Erast(Licks),EvList{k},Dura,{1});% makes collpased rasters. PSR1 is a cell(neurons)
        [PTH1,~,~]=MakePTH07(PSR1,repmat(N1,size(RAW(session).Erast(Licks),1),1),{2,0,BSIZE});%-----Fixed bin used here    

        PTH1smooth=[];
        for l=1:length(Tm)
           PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 smoothbins]):l).*fliplr(filterweights(1:min([l smoothbins+1]))))/sum(filterweights(1:min([l smoothbins+1])));                  
        end
        Lick.Cued(k).PSTHraw(session,:)=PTH1smooth;

    end
end

%plotting

for type=1:2
      
    psth1=nanmean(Lick.Cued(1+(type-1)*2).PSTHraw(:,Ishow),1);
    sem1=nanste(Lick.Cued(1+(type-1)*2).PSTHraw(:,Ishow),1); %calculate standard error of the mean
    up1=psth1+sem1;
    down1=psth1-sem1;

    psthE=nanmean(Lick.Cued(2+(type-1)*2).PSTHraw(:,Ishow),1);
    semE=nanste(Lick.Cued(2+(type-1)*2).PSTHraw(:,Ishow),1); %calculate standard error of the mean
    upE=psthE+semE;
    downE=psthE-semE;

    %plotting
    subplot(4,4,14+type);
    hold on;
    plot(time1,psth1,'Color',color{1+(type-1)*2},'linewidth',1);
    plot(time1,psthE,'Color',color{2+(type-1)*2},'linewidth',1);
    patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],color{1+(type-1)*2},'EdgeColor','none');alpha(0.5);
    patch([time1,time1(end:-1:1)],[upE,downE(end:-1:1)],color{2+(type-1)*2},'EdgeColor','none');alpha(0.5);
    plot([0 0],[-2 8],':','color','k','linewidth',0.75);
    axis([Xaxis(1) Xaxis(2) 0 8]);
    xlabel('Seconds from reward delivery');
    if type==1 title('Cued, predictive'); end
    if type==2 title('Cued, non-predictive'); end
    %legend('Suc trials','Mal trials');

end


%% preferences

%data
x=[1 2];

CPref{1}=[83.8 92.7]; %C3
CPref{2}=[68.7 75.7]; %C4
CPref{3}=[66.7 90.1]; %C9
CPref{4}=[96.0 70.7]; %C10
%CPPref{5}=[82.7  86.9]; %C2
%CPPref{6}=[77.4 90.2]; %C7

%not included in ephys
if strcmp(plot_region,'Both')
    IBPref{1}=[81.7 81.3];
    IBPref{2}=[77 73.9];
    IBPref{3}=[37.4 74.6];
    IBPref{4}=[71 95.6];
    IBPref{5}=[89 83.6];
    IBPref{6}=[94 93];
    IBPref{7}=[59.6 78.8];
    IBPref{8}=[88.8 85.8];
    IBPref{9}=[71.2 89];
    IBPref{10}=[87 76.6];
    IBPref{11}=[85.2 84.1];
elseif strcmp(plot_region,'VP')
    IBPref{1}=[59.6 78.8];
    IBPref{2}=[88.8 85.8];
    IBPref{3}=[71.2 89];
    IBPref{4}=[87 76.6];
    IBPref{5}=[85.2 84.1];
elseif strcmp(plot_region,'NA')
    IBPref{1}=[81.7 81.3];
    IBPref{2}=[77 73.9];
    IBPref{3}=[37.4 74.6];
    IBPref{4}=[71 95.6];
    IBPref{5}=[89 83.6];
    IBPref{6}=[94 93];
end

subplot(3,5,5);
hold on;
%do the first 2 first in order to get the legend right
plot(x,IBPref{1},'Marker','o','MarkerFaceColor','w','LineWidth',1,'Color','k');
plot(x,CPref{1},'Marker','o','MarkerFaceColor','k','LineWidth',1,'Color','k');

for i=2:length(IBPref)
    plot(x,IBPref{i},'Marker','o','MarkerFaceColor','w','LineWidth',1,'Color','k');
end

for i=2:length(CPref)
    plot(x,CPref{i},'Marker','o','MarkerFaceColor','k','LineWidth',1,'Color','k');
end

axis([0.7 2.2 0 100]);
title('Two-bottle choice');
xlabel('Initial       Final');
ylabel('% sucrose consumption');
plot([0 3],[50 50],':','color','k','linewidth',0.75);
set(gca,'xtick',[])
legend('Int / block','Cued','location','southeast');




end
