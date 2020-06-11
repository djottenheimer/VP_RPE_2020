clear all;


%parameters
window=[-5 20]; %seconds relative to reward delivery
binsize=0.2; %in seconds
bins=(window(2)-window(1))/binsize;
xaxis=linspace(window(1)+binsize/2,window(1)+binsize/2+(bins-1)*binsize,bins);

WOI=[0 10];
boi=(WOI(1)-window(1))/binsize+1:(WOI(2)-window(1))/binsize; 

address=['C:\Users\dottenh2\Dropbox\MATLAB\David\2R_blocks\ChR2PPBehavior']; %where are the running files located
%address=['/Users/David/Dropbox/MATLAB/David/2R_blocks/ChR2PPBehavior'];
files=dir([address,'\\*!*']);
%files=dir([address,'//*!*']);

for n=1:length(files)
    filename=fullfile(files(n).folder,files(n).name);
    file=fopen(filename);
    L=textscan(file,'%s','Delimiter',':');
    fclose('all');
    x=str2num(files(n).name(30:31)); %x is rat #
    
    %find start of A, start of C, and end of C
    Jstrt=find(strcmp('J',L{1,1})); %Cue
    Kstrt=find(strcmp('K',L{1,1})); %Sucrose Delivery
    Ostrt=find(strcmp('O',L{1,1})); %Port entry times
    Pstrt=find(strcmp('P',L{1,1})); %Port durations
    Qstrt=find(strcmp('Q',L{1,1})); %Laser stimulations
    Xstrt=find(strcmp('X',L{1,1})); %List (marks end of laser stims)
    
    %Sucrose delivery times
    sucrosedel=(L{1,1}(Kstrt+2:2:Ostrt-1));
    sucrosedeltime=[];
    for i=1:length(sucrosedel)
        medtext=textscan(char(sucrosedel(i,1)),'%f','Delimiter',' ');
        sdtime=medtext{1,1}(:,1);
        sucrosedeltime=(cat(1,sucrosedeltime,sdtime(isnan(sdtime)==0)));
    end
    sucrosedeltime(sucrosedeltime==0)=[]; %delete 0s
    
    %Port entry times
    pet=(L{1,1}(Ostrt+2:2:Pstrt-1));
    petime=[];
    for i=1:length(pet)
        medtext=textscan(char(pet(i,1)),'%f','Delimiter',' ');
        portentry=medtext{1,1}(:,1);
        petime=(cat(1,petime,portentry(isnan(portentry)==0)));
    end
    removes=petime==0;
    petime(petime==0)=[]; %delete 0s
    
    %Port entry durations
    durationt=(L{1,1}(Pstrt+2:2:Qstrt-1));
    durationtime=[];
    for i=1:length(durationt)
        medtext=textscan(char(durationt(i,1)),'%f','Delimiter',' ');
        duration=medtext{1,1}(:,1);
        durationtime=(cat(1,durationtime,duration(isnan(duration)==0)));
    end
    durationtime(removes)=[]; %delete 0s
    
    
    %Laser stims
    ST=(L{1,1}(Qstrt+2:2:Xstrt-1));
    STtime=[];
    for i=1:length(ST)
        medtext=textscan(char(ST(i,1)),'%f','Delimiter',' ');
        stime=medtext{1,1}(:,1);
        STtime=(cat(1,STtime,stime(isnan(stime)==0)));
    end
    
    STtime(STtime==0)=[]; %delete 0s
    
    stim_trials = [];
    for trl = 1:length(STtime)
        stim_trials(trl,1) =max(sucrosedeltime(sucrosedeltime<STtime(trl)));
    end
    
    nost_trials = sucrosedeltime(ismember(sucrosedeltime,stim_trials)==0);
    
    
    
    %port occupancy matrix stim
    stim_matrix=[];
    for trl=1:length(stim_trials)
        for bin=1:bins
            reltime=window(1)+(bin-1)*binsize;
            time=stim_trials(trl)+reltime;
            
            %check if this time point was during a port entry
            if sum(petime<=time)
                [entry_time,index] = max(petime(petime<=time));
                stim_matrix(trl,bin) = time <= (entry_time + durationtime(index));
            else
                stim_matrix(trl,bin) = 0;
            end
        end
    end
    
    %port occupancy matrix no stim
    nost_matrix=[];
    for trl=1:length(nost_trials)
        for bin=1:bins
            reltime=window(1)+(bin-1)*binsize;
            time=nost_trials(trl)+reltime;
            
            %check if this time point was during a port entry
            if sum(petime<=time)

            [entry_time,index] = max(petime(petime<=time));
            nost_matrix(trl,bin) = time <= (entry_time + durationtime(index));
            else
                nost_matrix(trl,bin) = 0;
            end
        end
    end
    
    stim_PSTH(x,:) = mean(stim_matrix);
    nost_PSTH(x,:) = mean(nost_matrix);

    
    
    
end

%% plotting
figure;

yfp=[4 6 10 11 13 14 21]; %at least 30 trials, both fibers in general vicinity
chr=[1 2 7 8 9 16 17 18 22 23 24]; %good placement, at least 30 trials

colors{2}=[0 0.3 1];
colors{1}=[0.4 0.4 0.4];
groups{1,1}=yfp;
groups{2,1}=chr;

for group=1:2
    subplot(2,2,2+group);
    hold on;
    psth1 = mean(nost_PSTH(groups{group,1},:));
    up1 = psth1 + nanste(nost_PSTH(groups{group,1},:),1);
    down1 = psth1 - nanste(nost_PSTH(groups{group,1},:),1);
    
    psth2 = mean(stim_PSTH(groups{group,1},:));
    up2 = psth2 + nanste(stim_PSTH(groups{group,1},:),1);
    down2 = psth2 - nanste(stim_PSTH(groups{group,1},:),1);
    
    patch([0 0 2 2],[0 1 1 0],[0.7 0.7 1],'edgecolor','none');
    a=plot(xaxis,psth1,'Color',colors{1},'linewidth',1);
    b=plot(xaxis,psth2,'Color',colors{2},'linewidth',1);

    patch([xaxis,xaxis(end:-1:1)],[up1,down1(end:-1:1)],colors{1},'EdgeColor','none');alpha(0.5);
    patch([xaxis,xaxis(end:-1:1)],[up2,down2(end:-1:1)],colors{2},'EdgeColor','none');alpha(0.5);  

    ylabel('Fraction of trials in port');
    xlabel('Seconds from reward delivery');
    plot([0 0],[0 1],':','color','k','linewidth',1);
    legend([a b],'No laser','Laser','location','northeast');
    if group==1 title('GFP'); end
    if group==2 title('ChR2'); end
end


%difference in time in port following 10s?
boi = (xaxis>= 15) & (xaxis<= 20);
nost_boi = mean(nost_PSTH(:,boi),2);
stim_boi = mean(stim_PSTH(:,boi),2);
p = signrank(nost_boi(chr),stim_boi(chr));
%% inhibition

clear all;


%parameters
window=[-5 20]; %seconds relative to reward delivery
binsize=0.2; %in seconds
bins=(window(2)-window(1))/binsize;
xaxis=linspace(window(1)+binsize/2,window(1)+binsize/2+(bins-1)*binsize,bins);

WOI=[0 10];
boi=(WOI(1)-window(1))/binsize+1:(WOI(2)-window(1))/binsize; 

address=['C:\Users\dottenh2\Dropbox\MATLAB\David\2R_blocks\OptoPPBehavior']; %where are the running files located
%address=['/Users/David/Dropbox/MATLAB/David/2R_blocks/OptoPPBehavior'];
files=dir([address,'\\*!*']);
%files=dir([address,'//*!*']);

for n=1:length(files)
    filename=fullfile(files(n).folder,files(n).name);
    file=fopen(filename);
    L=textscan(file,'%s','Delimiter',':');
    fclose('all');
    x=str2num(files(n).name(30:31)); %x is rat #
    
    %find start of A, start of C, and end of C
    Jstrt=find(strcmp('J',L{1,1})); %Cue
    Kstrt=find(strcmp('K',L{1,1})); %Sucrose Delivery
    Ostrt=find(strcmp('O',L{1,1})); %Port entry times
    Pstrt=find(strcmp('P',L{1,1})); %Port durations
    Qstrt=find(strcmp('Q',L{1,1})); %Laser stimulations
    Xstrt=find(strcmp('X',L{1,1})); %List (marks end of laser stims)
    
    %Sucrose delivery times
    sucrosedel=(L{1,1}(Kstrt+2:2:Ostrt-1));
    sucrosedeltime=[];
    for i=1:length(sucrosedel)
        medtext=textscan(char(sucrosedel(i,1)),'%f','Delimiter',' ');
        sdtime=medtext{1,1}(:,1);
        sucrosedeltime=(cat(1,sucrosedeltime,sdtime(isnan(sdtime)==0)));
    end
    sucrosedeltime(sucrosedeltime==0)=[]; %delete 0s
    
    %Port entry times
    pet=(L{1,1}(Ostrt+2:2:Pstrt-1));
    petime=[];
    for i=1:length(pet)
        medtext=textscan(char(pet(i,1)),'%f','Delimiter',' ');
        portentry=medtext{1,1}(:,1);
        petime=(cat(1,petime,portentry(isnan(portentry)==0)));
    end
    removes=petime==0;
    petime(petime==0)=[]; %delete 0s
    
    %Port entry durations
    durationt=(L{1,1}(Pstrt+2:2:Qstrt-1));
    durationtime=[];
    for i=1:length(durationt)
        medtext=textscan(char(durationt(i,1)),'%f','Delimiter',' ');
        duration=medtext{1,1}(:,1);
        durationtime=(cat(1,durationtime,duration(isnan(duration)==0)));
    end
    durationtime(removes)=[]; %delete 0s
    
    
    %Laser stims
    ST=(L{1,1}(Qstrt+2:2:Xstrt-1));
    STtime=[];
    for i=1:length(ST)
        medtext=textscan(char(ST(i,1)),'%f','Delimiter',' ');
        stime=medtext{1,1}(:,1);
        STtime=(cat(1,STtime,stime(isnan(stime)==0)));
    end
    
    STtime(STtime==0)=[]; %delete 0s
    
    stim_trials = [];
    for trl = 1:length(STtime)
        stim_trials(trl,1) =max(sucrosedeltime(sucrosedeltime<STtime(trl)));
    end
    
    nost_trials = sucrosedeltime(ismember(sucrosedeltime,stim_trials)==0);
    
    
    
    %port occupancy matrix stim
    stim_matrix=[];
    for trl=1:length(stim_trials)
        for bin=1:bins
            reltime=window(1)+(bin-1)*binsize;
            time=stim_trials(trl)+reltime;
            
            %check if this time point was during a port entry
            if sum(petime<=time)
                [entry_time,index] = max(petime(petime<=time));
                stim_matrix(trl,bin) = time <= (entry_time + durationtime(index));
            else
                stim_matrix(trl,bin) = 0;
            end
        end
    end
    
    %port occupancy matrix no stim
    nost_matrix=[];
    for trl=1:length(nost_trials)
        for bin=1:bins
            reltime=window(1)+(bin-1)*binsize;
            time=nost_trials(trl)+reltime;
            
            %check if this time point was during a port entry
            if sum(petime<=time)

            [entry_time,index] = max(petime(petime<=time));
            nost_matrix(trl,bin) = time <= (entry_time + durationtime(index));
            else
                nost_matrix(trl,bin) = 0;
            end
        end
    end
    
    stim_PSTH(x,:) = mean(stim_matrix);
    nost_PSTH(x,:) = mean(nost_matrix);

    
    
    
end

%% plotting

yfp=[1 2 6 9 10 15 23]; %at least 30 trials, both fibers in general vicinity
%arc=[3 4 7 8 11 13 14 17 18 21 22];
%arc=[3 4 11 13 18 22]; %latency effects
arc=[3 4 8 11 13 18 22]; %good placement, at least 30 trials

colors{2}=[0 0.6 0.4];
colors{1}=[0.4 0.4 0.4];
groups{1,1}=yfp;
groups{2,1}=arc;

for group=1:2
    subplot(2,2,group);
    hold on;
    psth1 = mean(nost_PSTH(groups{group,1},:));
    up1 = psth1 + nanste(nost_PSTH(groups{group,1},:),1);
    down1 = psth1 - nanste(nost_PSTH(groups{group,1},:),1);
    
    psth2 = mean(stim_PSTH(groups{group,1},:));
    up2 = psth2 + nanste(stim_PSTH(groups{group,1},:),1);
    down2 = psth2 - nanste(stim_PSTH(groups{group,1},:),1);
    
    patch([0 0 5 5],[0 1 1 0],[0.7 1 0.7],'edgecolor','none');
    a=plot(xaxis,psth1,'Color',colors{1},'linewidth',1);
    b=plot(xaxis,psth2,'Color',colors{2},'linewidth',1);

    patch([xaxis,xaxis(end:-1:1)],[up1,down1(end:-1:1)],colors{1},'EdgeColor','none');alpha(0.5);
    patch([xaxis,xaxis(end:-1:1)],[up2,down2(end:-1:1)],colors{2},'EdgeColor','none');alpha(0.5);  

    ylabel('Fraction of trials in port');
    xlabel('Seconds from reward delivery');
    plot([0 0],[0 1],':','color','k','linewidth',1);
    legend([a b],'No laser','Laser','location','northeast');
    if group==1 title('YFP'); end
    if group==2 title('ArchT'); end

end

