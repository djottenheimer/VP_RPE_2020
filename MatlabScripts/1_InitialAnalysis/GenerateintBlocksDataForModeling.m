load ('RAWintBlocks.mat');
RAW=RAWblocks;
runanalysis=1;

%get parameters
trialsback=10; %how many trials back to look
Baseline=[-11 -1]; %For normalizing activity
RDWindow=[0.75 1.95];
PEWindow=[-1.5 -0.5]; %relative to RD
CueWindow=[0 0.75];


%reset
NN=0;
EvMeanz=0;

Ninfo=[];
Nneurons=0;
%Finds the total number of neurons in 2R and marks them by region/session
for i=1:length(RAW)
    Ninfo=cat(1,Ninfo,RAW(i).Ninfo);
    Nneurons=Nneurons+size(RAW(i).Nrast,1);
end
for i=1:Nneurons
    Session=string(Ninfo(i,1));
    Name=char(Session);
    Region(i,1)=cellstr(Name(1:2)); %an array with each neurons region
    Rat(i,1)=cellstr(Name(1:3)); %an array with each neuron's rat
    Blocks(i,1)=cellstr(Name(4));
    Blocks12(i,1)=cellstr(Name(5));
end
CS.Blocks=strcmp('B',Blocks); %neurons from blocks sessions are marked 1, int is 0
CS.Blocks12=zeros(length(Blocks12),1); %start with all 0s
CS.Blocks12(strcmp('1',Blocks12))=1; %neurons from sessions starting with sucrose are 1
CS.Blocks12(strcmp('2',Blocks12))=2; %neurons from sessions starting with malt are 2


if runanalysis==1  
    for i=1:length(RAW) %loops through sessions
        
        %ADD IN STUFF
        
        %events
        EV1=strmatch('RD1P2',RAW(i).Einfo(:,2),'exact');
        EV2=strmatch('RD1P1',RAW(i).Einfo(:,2),'exact');
        EV3=strmatch('RD2P2',RAW(i).Einfo(:,2),'exact');
        EV4=strmatch('RD2P1',RAW(i).Einfo(:,2),'exact');
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
        rewards2=cat(2,RAW(i).Erast{R2,1}(:,1),zeros(length(RAW(i).Erast{R2,1}(:,1)),1));
        rewards=cat(1,rewards1,rewards2);
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
            Y=(rewspk(1:end,1)/(Window(1,2)-Window(1,1)));%-Bmean)/Bstd; %subtract baseline firing
            CS.RDHz{NN,1}=Y;

            %get trial by trial firing rate for all PE trials based on fixed window
            rewspk=0;
            Window=PEWindow;
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{RD},Window,{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>Window(1));
            end       
            Y=(rewspk(1:end,1)/(Window(1,2)-Window(1,1)));%-Bmean)/Bstd; %z-score
            CS.PEHz{NN,1}=Y;               

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
            Y=(rewspk(1:end,1)/(CueWindow(1,2)-CueWindow(1,1)));%-Bmean)/Bstd; %subtract baseline to keep it linear
            CS.CueHz{NN,1}=Y;
            
            
            %All Cues
            AllCues = RAW(i).Erast{Cue};
            rewspk=0;
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),AllCues,CueWindow,{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>CueWindow(1));
            end       
            Y=(rewspk(1:end,1)/(CueWindow(1,2)-CueWindow(1,1)));%-Bmean)/Bstd; %subtract baseline to keep it linear
            CS.CueHzAll{NN,1}=Y;

            CS.Predictors{NN,1}=X;
            CS.AllTrials{NN,1}=AllTrials;

            fprintf('Neuron # %d\n',NN);
        end
        

    end
end

%% Other data

CS.Rat=Rat;
CS.Region=Region;

save('ModData_intBlocks.mat','CS');