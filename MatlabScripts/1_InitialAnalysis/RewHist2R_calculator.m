%Looking at the average firing rate for a given window in each of 4
%current/previous reward conditions

%run linear model and stats? 1 is yes, 0 is no
runanalysis=1;

%get parameters
trialsback=10; %how many trials back to look
Baseline=[-11 -1]; %Relative to cue
RDWindow=[0.75 1.95];
PEWindow=[-1.5 -0.5]; %relative to RD
CueWindow=[0 0.75];



%for looking across time
bins=45;
binstart=-1.5; %in seconds
binsize=0.1; %in seconds

%save variables
RewHist2R.trialsback=trialsback;
RewHist2R.Baseline=Baseline;
RewHist2R.RDWindow=RDWindow;
RewHist2R.PEWindow=PEWindow;
RewHist2R.CueWindow=CueWindow;
RewHist2R.bins=bins;
RewHist2R.binstart=binstart;
RewHist2R.binsize=binsize;

%reset
NN=0;
EvMeanz=0;

if runanalysis==1  
    load ('RAWintBlocks.mat');
    RAW=RAWblocks;
    for i=1:length(RAW) %loops through sessions
        if strcmp('VP',RAW(i).Region(1:2)) && strcmp('I',RAW(i).Blocks) %only look at VP interspersed sessions
            %events
            EV1=strmatch('RD1P2',RAW(i).Einfo(:,2),'exact');
            EV2=strmatch('RD1P1',RAW(i).Einfo(:,2),'exact');
            EV3=strmatch('RD2P2',RAW(i).Einfo(:,2),'exact');
            EV4=strmatch('RD2P1',RAW(i).Einfo(:,2),'exact');
            RD=strmatch('RD',RAW(i).Einfo(:,2),'exact');
            Cue=strmatch('Cue',RAW(i).Einfo(:,2),'exact');
            R1=strmatch('R1withanylick',RAW(i).Einfo(:,2),'exact');
            R2=strmatch('R2withanylick',RAW(i).Einfo(:,2),'exact');

            % linear model for impact of previous rewards
            %reset
            X=[];
            Y=[];

            %set up the matrix with outcome identity for each session
            rewards1=cat(2,RAW(i).Erast{R1,1}(:,1),ones(length(RAW(i).Erast{R1,1}(:,1)),1));
            rewards2=cat(2,RAW(i).Erast{R2,1}(:,1),zeros(length(RAW(i).Erast{R2,1}(:,1)),1));
            rewards=cat(1,rewards1,rewards2);
            [rewards(:,1),ind]=sort(rewards(:,1));
            rewards(:,2)=rewards(ind,2);

            for k=trialsback+1:length(RAW(i).Erast{RD,1}(:,1))
                time=RAW(i).Erast{RD,1}(k,1);
                entry=find(rewards(:,1)==time);
                for m=1:trialsback+1
                    X(k-trialsback,m)=rewards(entry+1-m,2);
                end
            end

            for j= 1:length(RAW(i).Nrast) %Number of neurons within sessions

                NN=NN+1;
                basespk=0;

                %get mean baseline firing for all reward trials
                [Bcell1,B1n]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{Cue},Baseline,{2});% makes trial by trial rasters for baseline
                for y= 1:B1n
                    basespk(1,y)=sum(Bcell1{1,y}>Baseline(1));
                end

                Bhz=basespk/(Baseline(1,2)-Baseline(1,1));
                Bmean=nanmean(Bhz);
                Bstd=nanstd(Bhz);
                
                for k=1:bins
                    
                    Window=[binstart+(k-1)*binsize binstart+k*binsize];
                    rewspk=[];
                    
                    %get trial by trial firing rate for all reward trials
                    [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{RD},Window,{2});% makes trial by trial rasters for event
                    for y= 1:EVn
                        rewspk(y,1)=sum(EVcell{1,y}>Window(1));
                    end       
                    Y=((rewspk(trialsback+1:end,1)/(Window(1,2)-Window(1,1)))-Bmean)/Bstd; %normalize the activity to baseline
                                      
                    RewHist2R.BinHz{NN,k}=Y;


                end
                
                %get trial by trial firing rate for all reward trials based on fixed window
                rewspk=0;
                Window=RDWindow;
                [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{RD},Window,{2});% makes trial by trial rasters for event
                for y= 1:EVn
                    rewspk(y,1)=sum(EVcell{1,y}>Window(1));
                end       
                Y=((rewspk(trialsback+1:end,1)/(Window(1,2)-Window(1,1)))-Bmean)/Bstd; %z-score
                RewHist2R.RDHz{NN,1}=Y;
                
                %get trial by trial firing rate for all PE trials based on fixed window
                rewspk=0;
                Window=PEWindow;
                [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{RD},Window,{2});% makes trial by trial rasters for event
                for y= 1:EVn
                    rewspk(y,1)=sum(EVcell{1,y}>Window(1));
                end       
                Y=((rewspk(trialsback+1:end,1)/(Window(1,2)-Window(1,1)))-Bmean)/Bstd; %z-score
                RewHist2R.PEHz{NN,1}=Y;               
                
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
                Y=((rewspk(trialsback+1:end,1)/(CueWindow(1,2)-CueWindow(1,1)))-Bmean)/Bstd; %subtract baseline to keep it linear
                RewHist2R.CueHz{NN,1}=Y;
                
                %save the predictors
                RewHist2R.Predictors{NN,1}=X;


                fprintf('Neuron # %d\n',NN);
            end
        end

    end

end

save('RewHist2R.mat','RewHist2R');