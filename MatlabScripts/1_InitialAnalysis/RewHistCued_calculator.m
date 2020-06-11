%Looking at the average firing rate for a given window in each of 4
%current/previous reward conditions

%run linear model and stats? 1 is yes, 0 is no
runanalysis=1;

%which sessions are we looking at
IncRats=[3 4 9 10]; %needs to match 'Cued.m'
IncDays=11:20; %needs to match 'Cued.m'

%get parameters
trialsback=10; %how many trials back to look
Baseline=[-11 -1]; %Relative to cue
RDWindow=[0.75 1.95];
PEWindow=[-1.5 -0.5]; %relative to RD
CueWindow=[0 0.75];

%reset
NN=0;
EvMeanz=0;

if runanalysis==1 
    load('RAWCued.mat')
    RAW=RAWCued;
    for i=1:length(RAW) %loops through sessions
        if any(IncRats==RAW(i).Rat) && any(IncDays==RAW(i).Day) %only look at included sessions
            %events
            Cues=strmatch('Cues',RAW(i).Einfo(:,2),'exact');
            RDC=strmatch('RDC',RAW(i).Einfo(:,2),'exact');
            RDA=strmatch('RDA',RAW(i).Einfo(:,2),'exact');
            R1=strmatch('R1withanylick',RAW(i).Einfo(:,2),'exact');
            R2=strmatch('R2withanylick',RAW(i).Einfo(:,2),'exact');

            %% linear model for impact of previous rewards
            %reset
            X=[];
            Y=[];

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
            
            %create results key
            for k=trialsback+1:length(RDall(:,1))
                time=RDall(k,1);
                entry=find(rewards(:,1)==time);
                for m=1:trialsback+1
                    X(k-trialsback,m)=rewards(entry+1-m,2);
                end
                
                %add a 1 to next column if predicted sucrose, column after
                %that if predicted maltodextrin
                if RDall(k,2)==1
                    if X(k-trialsback,1)==1
                        X(k-trialsback,trialsback+2)=1;
                        X(k-trialsback,trialsback+3)=0;
                    else
                        X(k-trialsback,trialsback+2)=0;
                        X(k-trialsback,trialsback+3)=1;
                    end
                else
                    X(k-trialsback,trialsback+2:trialsback+3)=0;
                end
                %X(k-trialsback,trialsback+4)=(rewards(entry,1)-rewards(entry-1,1))/60; %time since last reward, in minutes
                        
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
                Y=((rewspk(trialsback+1:end,1)/(RDWindow(1,2)-RDWindow(1,1)))-Bmean)/Bstd; %normalize the activity to baseline
                
                
                RewHistCued.RDHz{NN,1}=Y;
                
                %get trial by trial firing rate for all port entries
                rewspk=[];
                [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RDall(:,1),PEWindow,{2});% makes trial by trial rasters for event
                for y= 1:EVn
                    rewspk(y,1)=sum(EVcell{1,y}>PEWindow(1));
                end       
                Y=((rewspk(trialsback+1:end,1)/(PEWindow(1,2)-PEWindow(1,1)))-Bmean)/Bstd; %normalize the activity to baseline
                
                RewHistCued.PEHz{NN,1}=Y;
              
                
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
                Y=((rewspk(trialsback+1:end,1)/(CueWindow(1,2)-CueWindow(1,1)))-Bmean)/Bstd; %subtract baseline to keep it linear
                RewHistCued.CueHz{NN,1}=Y;
                
                
                
                %save the predictors
                RewHistCued.Predictors{NN,1}=X;


                fprintf('Neuron # %d\n',NN);
            end
        end

    end
end

save('RewHistCued.mat','RewHistCued');

