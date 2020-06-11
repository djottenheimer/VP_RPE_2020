%PSTHs and some simple analysis

clear all; clc;
global Dura Baseline Tm Tbase BSIZE Tbin
tic
if ~exist('RAWCued'), load ('RAWCued.mat'); end
RAW=RAWCued;

%Main settings
SAVE_FLAG=1;
BSIZE=0.01; %Do not change
Dura=[-4 10]; Tm=Dura(1):BSIZE:Dura(2);
%Baseline=[-22 0]; Tbase=Baseline(1):BSIZE:Baseline(2); %now defined line 98
Tbin=-0.5:0.005:0.5; %window used to determine the optimal binsize
Tbase=-11:BSIZE:-1; %for fixed baseline for z-score
PStat=0.05; %for comparing pre versus post windows, or event A versus event B
MinNumTrials=5;
R=[];R.Ninfo={};NN=0;Nneurons=0;

%for fixed bin size
BinSize=0.01; %in seconds

%which sessions included
IncRats=[3 4 9 10];
IncDays=11:20;

%Smoothing PSTHs
%smoothing filter
smoothbins=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',6.6); %std=3.98
filterweights=pdf(halfnormal,[0:smoothbins]);

%old smoothing
% Smoothing=1; %0 for raw and 1 for smoothing
% SmoothTYPE='lowess';
% SmoothSPAN=7;
% if Smoothing~=1, SmoothTYPE='NoSmoothing';SmoothSPAN=NaN; end

% List of events to analyze and analysis windows EXTRACTED from excel file
spreadsheet=importdata('Cued.xlsx');
Erefnames=spreadsheet.textdata(3:11,1);
prewin =  spreadsheet.data(1:9,1:2);
postwin = spreadsheet.data(1:9,3:4);

%%
%Finds the total number of neurons accross all sessions
for i=1:length(RAW)
    if any(IncRats==RAW(i).Rat) && any(IncDays==RAW(i).Day) %only look at included sessions
        R.Ninfo=cat(1,R.Ninfo,RAW(i).Ninfo);
        Nneurons=Nneurons+size(RAW(i).Nrast,1);
    end
end

R.Ninfo=cat(2, R.Ninfo, mat2cell([1:Nneurons]',ones(1,Nneurons),1));
R.Erefnames= Erefnames;


% preallocating
R.Param.Tm=Tm;
R.Param.Tbin=Tbin;
R.Param.Dura=Dura;
R.Param.Baseline=Baseline;
R.Param.PStat=PStat;
R.Param.MinNumTrials=MinNumTrials;
R.Param.path=path;
R.Param.prewin=prewin;
R.Param.postwin=postwin;
R.Param.SmoothTYPE='HalfNormal';
R.Param.SmoothSPAN=std(halfnormal);

for k=1:length(Erefnames)
    R.Ev(k).PSTHraw(1:Nneurons,1:length(Tm))=NaN(Nneurons,length(Tm));
    R.Ev(k).PSTHz(1:Nneurons,1:length(Tm))=NaN(Nneurons,length(Tm));
    R.Ev(k).Meanraw(1:Nneurons,1)=NaN;
    R.Ev(k).Meanz(1:Nneurons,1)=NaN;
    R.Ev(k).BW(1:Nneurons,1)=NaN;
    R.Ev(k).signrank(1:Nneurons,1)=NaN;
    R.Ev(k).RespDir(1:Nneurons,1)=NaN;
    R.Ev(k).NumberTrials(1:Nneurons,1)=NaN;
end

%% runs the main routine

figure;

for i=1:length(RAW) %loops through sessions
    if any(IncRats==RAW(i).Rat) && any(IncDays==RAW(i).Day) %only look at included sessions
        for j= 1:size(RAW(i).Nrast,1) %Number of neurons per session
            NN=NN+1; %neuron counter
            R.Day(NN,1)=RAW(i).Day;
            R.Rat(NN,1)=RAW(i).Rat;

            %use the same baseline for z-scoring every event
            %baseline period is before every cue in the whole session
            bltrlhz=[];
            EvInd=strcmp('Cues',RAW(i).Einfo(:,2));
            [PSRbl,N0]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{EvInd},[-11 -1],{2});% makes collapsed rasters for baseline for use in normalizing
            for l=1:length(PSRbl)
                bltrlhz(l,1)=sum((PSRbl{l}<-1 & PSRbl{l}>-11)/10);
            end
            bmean=nanmean(bltrlhz(2:2:end)); %2:2:end is to match the number of trials in int/blocks sessions
            bstd=nanstd(bltrlhz(2:2:end));

            for k=1:length(Erefnames) %loops thorough the events
                EvInd=strcmp(Erefnames(k),RAW(i).Einfo(:,2)); %find the event id number from RAW

                if sum(EvInd)==0
                    fprintf('HOWDY, CANT FIND EVENTS FOR ''%s''\n',Erefnames{k});
                end

                R.Ev(k).NumberTrials(NN,1)=length(RAW(i).Erast{EvInd});
                if  ~isempty(EvInd) && R.Ev(k).NumberTrials(NN,1)>=MinNumTrials %avoid analyzing sessions where that do not have enough trials
                    [PSR1,N1]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{EvInd},Dura,{1});% makes collpased rasters. PSR1 is a cell(neurons)
                    [PSR2,N2]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{EvInd},Dura,{2});% makes trial by trial rasters. PSR1 is a cell(neurons, trials)
                    [PSR3,~]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{EvInd},prewin(k,:),{2});% makes trial by trial rasters for baseline for use in response detection
                    if ~isempty(PSR1{1}) %to avoid errors, added on 12/28 2011
                        %optimal bin size
                        %[PTH1,BW1,~]=MakePTH07(PSR1,repmat(N1, size(RAW(i).Nrast{j},1),1),{2,1});%-----DP used here
                        %fixed bin size
                        [PTH1,BW1,~]=MakePTH07(PSR1,repmat(N1, size(RAW(i).Nrast{j},1),1),{2,0,BinSize});%fix bin used here

                        %------------- Fills the R.Ev(k) fields --------------
                        R.Ev(k).BW(NN,1)=BW1;

                        %new smoothing
                        for l=1:length(Tm)
                           PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 smoothbins]):l).*fliplr(filterweights(1:min([l smoothbins+1]))))/sum(filterweights(1:min([l smoothbins+1])));                  
                        end
                        R.Ev(k).PSTHraw(NN,:)=PTH1smooth;

                        %normalize already smoothed activity
                        for l=1:length(PTH1smooth)
                            R.Ev(k).PSTHz(NN,l)=(PTH1smooth(l)-bmean)/bstd;
                        end

                        %------------------ firing (in Hz) per trial in pre/post windows ------------------
                        %used to make the between events comparisons and Response detection in a single window----
                        ev(k).pre=NaN(size(RAW(i).Erast{EvInd},1),1);
                        ev(k).post=NaN(size(RAW(i).Erast{EvInd},1),1);
                        for m=1:size(RAW(i).Erast{EvInd},1) %loops through trials
                            ev(k).pre(m)=sum(PSR3{m}<prewin(k,2) & PSR3{m}>prewin(k,1))/(prewin(k,2)-prewin(k,1)); %CHANGED FROM PSR2 to PSR0 here 10/24/17
                            ev(k).post(m)=sum(PSR2{m}<postwin(k,2) & PSR2{m}>postwin(k,1))/(postwin(k,2)-postwin(k,1));
                        end

                        R.Ev(k).Meanraw(NN,1)=nanmean(ev(k).post);
                        R.Ev(k).Meanz(NN,1)=(nanmean(ev(k).post)-bmean)/bstd;

                        %---------------------------Response detection w/ SignRank on pre/post windows---------------------------


                        R.Ev(k).signrank(NN,1)=signrank(ev(k).pre, ev(k).post); %Signrank used here because it is a dependant sample test
                        if R.Ev(k).signrank(NN,1)<PStat
                            R.Ev(k).RespDir(NN,1)=sign(mean(ev(k).post)-mean(ev(k).pre));
                        else R.Ev(k).RespDir(NN,1)=0;
                        end


                    end %if ~isempty(PSR0{1}) || ~isempty(PSR1{1})
                end %if EvInd=0 OR n(trials) < MinNumTrials fills with NaN
            end %Events

            fprintf('Neuron ID # %d\n',NN);
        end %neurons: FOR j= 1:size(RAW(i).Nrast,1)
    end
end %sessions: FOR i=1:length(RAW)

R_cued=R;

if SAVE_FLAG
    save('R_cued_all.mat','R_cued')
end

toc

