%run analysis?
runanalysis=1;

%which sessions go in
%there's a logical later that highlights included neurons
IncRats=[2 3 4 7 9 10];
IncDays=8:20;

%get parameters
trialsback=10; %how many trials back to look
Baseline=[-11 -1]; %Relative to cue
RDWindow=[0.75 1.95];
PEWindow=[-1.5 -0.5]; %relative to RD
CueWindow=[0 0.3];

%%

%reset
NN=0;
EvMeanz=0;

if runanalysis==1  
    load ('RAWCued.mat');
    RAW=RAWCued;
    CD.Ninfo={};
    for i=1:length(RAW) %loops through sessions 
        if any(IncRats==RAW(i).Rat) && any(IncDays==RAW(i).Day) %only look at included sessions
            
            %store info about each neuron's channel
            CD.Ninfo=cat(1,CD.Ninfo,RAW(i).Ninfo);
            
            %events
            Cues=strmatch('Cues',RAW(i).Einfo(:,2),'exact');
            Cue1=strmatch('Cue1',RAW(i).Einfo(:,2),'exact');
            Cue2=strmatch('Cue2',RAW(i).Einfo(:,2),'exact');
            Cue3=strmatch('Cue3',RAW(i).Einfo(:,2),'exact');
            RDC=strmatch('RDC',RAW(i).Einfo(:,2),'exact');
            RDA=strmatch('RDA',RAW(i).Einfo(:,2),'exact');
            R1=strmatch('R1withanylick',RAW(i).Einfo(:,2),'exact');
            R2=strmatch('R2withanylick',RAW(i).Einfo(:,2),'exact');

            %% linear model for impact of previous rewards
            %reset
            Y=[];
            AllTrials=[];
            PrCue=[];
            

            %set up the matrix with outcome identity for each session
            rewards1=cat(2,RAW(i).Erast{R1,1}(:,1),ones(length(RAW(i).Erast{R1,1}(:,1)),1));
            rewards2=cat(2,RAW(i).Erast{R2,1}(:,1),zeros(length(RAW(i).Erast{R2,1}(:,1)),1));
            rewards=cat(1,rewards1,rewards2);
            [rewards(:,1),ind]=sort(rewards(:,1));
            rewards(:,2)=rewards(ind,2);
            
            %put both conditions together
            RDCkey=cat(2,RAW(i).Erast{RDC,1}(:,1),ones(length(RAW(i).Erast{RDC,1}(:,1)),1)); %1 if predictive trial
            RDAkey=cat(2,RAW(i).Erast{RDA,1}(:,1),zeros(length(RAW(i).Erast{RDA,1}(:,1)),1)); %0 if not
            RDall=cat(1,RDCkey,RDAkey);
            [RDall(:,1),ind]=sort(RDall(:,1));
            RDall(:,2)=RDall(ind,2);
            
            AllTrials(:,1)=rewards(:,2);
            AllTrials(:,2)=0;
            
            
            X=NaN(length(RDall),trialsback+1);
            %create results key
            for k=1:length(RDall(:,1))
                time=RDall(k,1);
                entry=find(round(rewards(:,1))==round(time));
                for m=1:trialsback+1
                    if entry+1-m>0
                        X(k,m)=rewards(entry+1-m,2);
                    end
                end
                
                %add a 1 to next column if predicted sucrose, column after
                %that if predicted maltodextrin
                if RDall(k,2)==1
                    if X(k,1)==1
                        X(k,trialsback+2)=1;
                        X(k,trialsback+3)=0;
                        PrCue(k,1)=1;
                        PrCue(k,2)=0;
                    else
                        X(k,trialsback+2)=0;
                        X(k,trialsback+3)=1;
                        PrCue(k,1)=0;
                        PrCue(k,2)=1;
                    end
                else
                    X(k,trialsback+2:trialsback+3)=0;
                    PrCue(k,1:2)=0;
                end
                      
                AllTrials(entry,2)=1;
                
            end
            
            %get cue for each trial
            cue1s=RAW(i).Erast{Cue1,1}(:,1);
            cue2s=RAW(i).Erast{Cue2,1}(:,1);
            cue3s=RAW(i).Erast{Cue3,1}(:,1);
            
            PrCueAllTrials=zeros(length(rewards),3);
            for j=1:length(rewards)
                mostrecentcue(1)=max([cue1s(cue1s<rewards(j,1))' 0]);
                mostrecentcue(2)=max([cue2s(cue2s<rewards(j,1))' 0]);
                mostrecentcue(3)=max([cue3s(cue3s<rewards(j,1))' 0]);
                [~,mostrecent]=max(mostrecentcue);
                PrCueAllTrials(j,mostrecent)=1;
            end
            
            
             for j= 1:length(RAW(i).Nrast) %Number of neurons within sessions
 
                 NN=NN+1;
                basespk=0;
                
                %get mean baseline firing for all trials
                [Bcell1,B1n]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{Cues},Baseline,{2});% makes trial by trial rasters for baseline
                for y= 1:B1n
                    basespk(1,y)=sum(Bcell1{1,y}>Baseline(1));
                end

                Bhz=basespk/(Baseline(1,2)-Baseline(1,1));
                Bmean=nanmean(Bhz);
                Bstd=nanstd(Bhz);

                %get trial by trial firing rate for all reward trials
                rewspk=[];
                [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RDall(:,1),RDWindow,{2});% makes trial by trial rasters for event
                for y= 1:EVn
                    rewspk(y,1)=sum(EVcell{1,y}>RDWindow(1));
                end       
                Y=(rewspk(1:end,1)/(RDWindow(1,2)-RDWindow(1,1)));%-Bmean)/Bstd; %normalize the activity to baseline
                
                
                CD.RDHz{NN,1}=Y;
                
                %get trial by trial firing rate for all port entries
                rewspk=[];
                [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RDall(:,1),PEWindow,{2});% makes trial by trial rasters for event
                for y= 1:EVn
                    rewspk(y,1)=sum(EVcell{1,y}>PEWindow(1));
                end       
                Y=(rewspk(1:end,1)/(PEWindow(1,2)-PEWindow(1,1)));%-Bmean)/Bstd; %normalize the activity to baseline
                
                CD.PEHz{NN,1}=Y;
              
                
                %get trial by trial firing rate for all cues
                CueList=[];
                rewspk=[];
                
                for k=1:length(RDall(:,1))
                    CueList(k,1)=RAW(i).Erast{Cues}(max(find(RAW(i).Erast{Cues}<RDall(k,1))),1);
                end
                
                [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),CueList,CueWindow,{2});% makes trial by trial rasters for event
                for y= 1:EVn
                    rewspk(y,1)=sum(EVcell{1,y}>CueWindow(1));
                end       
                Y=(rewspk(1:end,1)/(CueWindow(1,2)-CueWindow(1,1)));%-Bmean)/Bstd; %subtract baseline to keep it linear
                CD.CueHz{NN,1}=Y;
                
                
                
                %save the predictors
                CD.Predictors{NN,1}=X;
                CD.AllTrials{NN,1}=AllTrials;
                CD.PredCue{NN,1}=PrCue;
                CD.PredCueAllTrials{NN,1}=PrCueAllTrials;
                
                CD.Day(NN,1)=RAW(i).Day;
                CD.Rat(NN,1)=RAW(i).Rat;
                
                


                fprintf('Neuron # %d\n',NN);
            end
        end

    end
end

%% Other data

%which neurons are plotted
DOI = 11:20;
ROI = [3 4 9 10];
CD.IncludedNeurons=ChooseNs(CD,DOI,ROI);

save('ModData_Cued.mat','CD');