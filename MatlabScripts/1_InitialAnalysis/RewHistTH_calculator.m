%Looking at the average firing rate for a given window in each of 4
%current/previous reward conditions
MaltValue=0.75;

load ('RAWTH.mat');
RAW=RAWTH;

%run linear model and stats? 1 is yes, 0 is no
runanalysis=1;


%get parameters
trialsback=10; %how many trials back to look
Baseline=[-11 -1]; %For normalizing activity
RDWindow=[0.75 1.95];
PEWindow=[-1.5 -0.5]; %relative to RD
CueWindow=[0 0.75];

%save variables
RewHistTH.trialsback=trialsback;
RewHistTH.Baseline=Baseline;
RewHistTH.RDWindow=RDWindow;
RewHistTH.PEWindow=PEWindow;
RewHistTH.CueWindow=CueWindow;

%reset
NN=0;


if runanalysis==1  
    for i=1:length(RAW) %loops through sessions
        %events
        RD=strmatch('RD',RAW(i).Einfo(:,2),'exact');
        R1=strmatch('R1withanylick',RAW(i).Einfo(:,2),'exact');
        R2=strmatch('R2withanylick',RAW(i).Einfo(:,2),'exact');
        R3=strmatch('R3withanylick',RAW(i).Einfo(:,2),'exact');
        Cue=strmatch('Cue',RAW(i).Einfo(:,2),'exact');

        %% linear model for impact of previous rewards
        %reset
        X=NaN(length(RAW(i).Erast{RD,1}(:,1)),trialsback+1);
        Y=[];
        AllTrials=[];

        %set up the matrix with outcome identity for each session
        rewards1=cat(2,RAW(i).Erast{R1,1}(:,1),ones(length(RAW(i).Erast{R1,1}(:,1)),1));
        rewards2=cat(2,RAW(i).Erast{R2,1}(:,1),MaltValue*ones(length(RAW(i).Erast{R2,1}(:,1)),1));
        rewards3=cat(2,RAW(i).Erast{R3,1}(:,1),zeros(length(RAW(i).Erast{R3,1}(:,1)),1));
        rewards=cat(1,rewards1,rewards2,rewards3);
        [rewards(:,1),ind]=sort(rewards(:,1));
        rewards(:,2)=rewards(ind,2);

        AllTrials(:,1)=rewards(:,2);
        AllTrials(:,2)=0;

        for k=1:length(RAW(i).Erast{RD,1}(:,1))
            time=RAW(i).Erast{RD,1}(k,1);
            entry=find(round(rewards(:,1))==round(time));
            for m=1:trialsback+1
                if entry+1-m>0
                    X(k,m)=rewards(entry+1-m,2);
                end
            end

            AllTrials(entry,2)=1;

        end



        for j= 1:length(RAW(i).Nrast) %Number of neurons within sessions

            NN=NN+1;
            rewspk=0;
            basespk=0;

            %get mean baseline firing for all reward trials
            [Bcell1,B1n]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{Cue},Baseline,{2});% makes trial by trial rasters for baseline
            for y= 1:B1n
                basespk(1,y)=sum(Bcell1{1,y}>Baseline(1));
            end

            Bhz=basespk/(Baseline(1,2)-Baseline(1,1));
            Bmean=nanmean(Bhz);
            Bstd=nanstd(Bhz);

            %get trial by trial firing rate for all reward trials
            Window=RDWindow;
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{RD},Window,{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>Window(1));
            end       
            Y=((rewspk(1:end,1)/(Window(1,2)-Window(1,1)))-Bmean)/Bstd; %subtract baseline firing
            RewHistTH.RDHz{NN,1}=Y(trialsback+1:end);

            %get trial by trial firing rate for all PE trials based on fixed window
            rewspk=0;
            Window=PEWindow;
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{RD},Window,{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>Window(1));
            end       
            Y=((rewspk(1:end,1)/(Window(1,2)-Window(1,1)))-Bmean)/Bstd; %z-score
            RewHistTH.PEHz{NN,1}=Y(trialsback+1:end);               

            %get trial by trial firing rate for all cue trials based on fixed window
            CueList=[];
            for k=1:length(RAW(i).Erast{RD})
                CueList(k,1)=RAW(i).Erast{Cue}(max(find(RAW(i).Erast{Cue}<RAW(i).Erast{RD}(k,1))),1);
            end

            Window=CueWindow;
            rewspk=0;
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),CueList,CueWindow,{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>CueWindow(1));
            end       
            Y=((rewspk(1:end,1)/(CueWindow(1,2)-CueWindow(1,1)))-Bmean)/Bstd; %subtract baseline to keep it linear
            RewHistTH.CueHz{NN,1}=Y(trialsback+1:end);

            RewHistTH.Predictors{NN,1}=X(trialsback+1:end,:);
            RewHistTH.AllTrials{NN,1}=AllTrials;

            fprintf('Neuron # %d\n',NN);
        end
    end


end


save('RewHistTH.mat','RewHistTH');