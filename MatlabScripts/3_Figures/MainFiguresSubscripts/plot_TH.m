function stats=plot_TH(Colors, R, RewHist, bm_RD, bm_cue, RAW, os, plot_sub_figures, bm_RD_BIC)
%% behavior

if plot_sub_figures

    figure;

    subplot(3,1,1)
    session_to_use = 3;

    session_end = strcmp('MedEnd',RAW(session_to_use(1)).Einfo(:,2)); %finds the index of the event
    stop = RAW(session_to_use(1)).Erast{session_end,1};

    line([0,3+stop(end)/60],[2.3,2.3],'Color','k')

    %get events for this session
    RD1 = strcmp('RD1',RAW(session_to_use).Einfo(:,2)); %finds the index of the event
    RD2 = strcmp('RD2',RAW(session_to_use).Einfo(:,2)); %finds the index of the event
    RD3 = strcmp('RD3',RAW(session_to_use).Einfo(:,2)); %finds the index of the event

    suc_delivery = RAW(session_to_use).Erast{RD1,1}/60;
    malt_delivery = RAW(session_to_use).Erast{RD2,1}/60;
    wat_delivery = RAW(session_to_use).Erast{RD3,1}/60;

    for trial = 1:length(suc_delivery)
        a=line([suc_delivery(trial),suc_delivery(trial)],[2.2,2.8],'Color',Colors('sucrose'));
    end
    for trial = 1:length(malt_delivery)
        b=line([malt_delivery(trial),malt_delivery(trial)],[2,2.6],'Color',Colors('maltodextrin'));
    end
    for trial = 1:length(wat_delivery)
        c=line([wat_delivery(trial),wat_delivery(trial)],[1.8,2.4],'Color',Colors('water'));
    end


    text(90,3,'Sucrose, unpredictable','Color',Colors('sucrose'));
    text(90,2.25,'Maltodextrin, unpredictable','Color',Colors('maltodextrin'));
    text(90,1.5,'Water, unpredictable','Color',Colors('water'));




    xticks([0 20 40 60 80]);
    axis([0 110 0.9 4]);
    ax1 = gca;   
    ax1.YAxis.Visible = 'off'; 
    xlabel('Session time (minutes)')



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
    for session = 1:length(RAW)

        % Make logical index vectors to call every event type
        R1 = strcmp('RD1',RAW(session).Einfo(:,2)); %finds the index of the event
        R2 = strcmp('RD2',RAW(session).Einfo(:,2)); %finds the index of the event
        R3 = strcmp('RD3',RAW(session).Einfo(:,2)); %finds the index of the event
        Licks = strcmp('Licks',RAW(session).Einfo(:,2)); %finds the index of the event

        EvList{1}=RAW(session).Erast{R1};
        EvList{2}=RAW(session).Erast{R2};
        EvList{3}=RAW(session).Erast{R3};

        for k=1:3
            [PSR1,N1]=MakePSR04(RAW(session).Erast(Licks),EvList{k},Dura,{1});% makes collpased rasters. PSR1 is a cell(neurons)
            [PTH1,~,~]=MakePTH07(PSR1,repmat(N1,size(RAW(session).Erast(Licks),1),1),{2,0,BSIZE});%-----Fixed bin used here    

            PTH1smooth=[];
            for l=1:length(Tm)
               PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 smoothbins]):l).*fliplr(filterweights(1:min([l smoothbins+1]))))/sum(filterweights(1:min([l smoothbins+1])));                  
            end
            Lick(k).PSTHraw(session,:)=PTH1smooth;

        end
    end

    %plotting
    Xaxis=[-2 12];
    Ishow=find(Tm>=Xaxis(1) & Tm<=Xaxis(2));
    time1=Tm(Ishow);



    Sel=ones(length(RAW),1)==1;

    psth1=nanmean(Lick(1).PSTHraw(Sel,Ishow),1);
    sem1=nanste(Lick(1).PSTHraw(Sel,Ishow),1); %calculate standard error of the mean
    up1=psth1+sem1;
    down1=psth1-sem1;

    psth2=nanmean(Lick(2).PSTHraw(Sel,Ishow),1);
    sem2=nanste(Lick(2).PSTHraw(Sel,Ishow),1); %calculate standard error of the mean
    up2=psth2+sem2;
    down2=psth2-sem2;

    psth3=nanmean(Lick(3).PSTHraw(Sel,Ishow),1);
    sem3=nanste(Lick(3).PSTHraw(Sel,Ishow),1); %calculate standard error of the mean
    up3=psth3+sem3;
    down3=psth3-sem3;

    %plotting
    subplot(3,3,9);
    hold on;
    plot(time1,psth1,'Color',Colors('sucrose'),'linewidth',1);
    plot(time1,psth2,'Color',Colors('maltodextrin'),'linewidth',1);
    plot(time1,psth3,'Color',Colors('water'),'linewidth',1);
    patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(0.5);
    patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(0.5);
    patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('water'),'EdgeColor','none');alpha(0.5);
    plot([0 0],[-2 8],':','color','k','linewidth',0.75);
    axis([Xaxis(1) Xaxis(2) 0 8]);
    xlabel('Seconds from reward delivery');
    ylabel('Mean lick rate (licks/s)');
    legend('Suc trials','Mal trials','Wat trials');

end

%% main figure

xaxis=[-2 5];
time_indicies=find(R.Param.Tm>=xaxis(1) & R.Param.Tm<=xaxis(end)); %returns the indicies where these are both true
activity_xaxis=R.Param.Tm(time_indicies);
suc_delivery=strcmp('RD1', R.Erefnames);
malt_delivery=strcmp('RD2', R.Erefnames);
wat_delivery=strcmp('RD3', R.Erefnames);

THcolors{1,1}=[1  0.7  0.4];
THcolors{2,1}=[.8  0.6  0.3];
THcolors{3,1}=[.5  0.3  0];
THcolors{4,1}=[1  0.5  0.8];
THcolors{5,1}=[.8  0.4  .7];
THcolors{6,1}=[.5  0.1  .3];
THcolors{7,1}=[0.2  0.8  0.9];
THcolors{8,1}=[0.1  0.7  0.8];
THcolors{9,1}=[0 0.5  0.6];

HistWindow=RewHist.RDWindow;
CueWindow=RewHist.CueWindow;

figure;

%% raster
alph=0.3;
%raster
window=[0 3];
sample_neuron_number=228; %237, 228

subplot(1,4,1);
hold on;

sessions=unique(R.Ninfo(:,1));
session=R.Ninfo(sample_neuron_number,1);
sample_neuron_session=strcmp(session,sessions);
neurons=R.Ninfo(strcmp(session,R.Ninfo(:,1)),2);
sample_neuron_session_number=strcmp(R.Ninfo(sample_neuron_number,2),neurons);


RD=strcmp('RD', RAW(sample_neuron_session).Einfo(:,2));
RDtimes=RAW(sample_neuron_session).Erast{RD,1};
[~,order]=sort(os(sample_neuron_number).mod_RD.base.RPEs);
RDtimes=RDtimes(order);

patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[0 length(RDtimes)*11 length(RDtimes)*11 0],Colors('gray'),'edgecolor','none');alpha(alph);
PlotRaster(RAW(sample_neuron_session).Nrast{sample_neuron_session_number,1},RDtimes,window);
ylabel('Trials, sorted by RPE (low to high)');
yticks([]);
xticks(0:3);
xticklabels(0:3);
% plot([HistWindow(1) HistWindow(1)],[0 length(RDtimes)*11],'color',[0 0.5 0.3],'linewidth',0.85);
% plot([HistWindow(2) HistWindow(2)],[0 length(RDtimes)*11],'color',[0 0.5 0.3],'linewidth',0.85);
xlabel('Seconds from reward delivery');
xlim(window);

disp('Raster')
%% activity plots
ymax=5;
alph=0.3;
Sel=R.Rat>0; %all neurons

%plot sucrose preferring neurons to sucrose
psth1=nanmean(R.Ev(suc_delivery).PSTHz(Sel,time_indicies),1); 
sem1=nanste(R.Ev(suc_delivery).PSTHz(Sel,time_indicies),1); %calculate standard error of the mean
up1=psth1+sem1;
down1=psth1-sem1;

%plot sucrose preferring neurons to malt
psth2=nanmean(R.Ev(malt_delivery).PSTHz(Sel,time_indicies),1); 
sem2=nanste(R.Ev(malt_delivery).PSTHz(Sel,time_indicies),1); %calculate standard error of the mean
up2=psth2+sem2;
down2=psth2-sem2;

%plot sucrose preferring neurons to malt
psth3=nanmean(R.Ev(wat_delivery).PSTHz(Sel,time_indicies),1); 
sem3=nanste(R.Ev(wat_delivery).PSTHz(Sel,time_indicies),1); %calculate standard error of the mean
up3=psth3+sem3;
down3=psth3-sem3;

subplot(3,8,3:5)
hold on
patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-3 ymax ymax -3],Colors('gray'),'edgecolor','none');alpha(alph);
a=plot(activity_xaxis,psth1,'Color',Colors('sucrose'),'linewidth',1);
b=plot(activity_xaxis,psth2,'Color',Colors('maltodextrin'),'linewidth',1);
c=plot(activity_xaxis,psth3,'Color',Colors('water'),'linewidth',1);
patch([activity_xaxis,activity_xaxis(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(alph);
patch([activity_xaxis,activity_xaxis(end:-1:1)],[up2,down2(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(alph);
patch([activity_xaxis,activity_xaxis(end:-1:1)],[up3,down3(end:-1:1)],Colors('water'),'EdgeColor','none');alpha(alph);
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
plot([0 0],[-3 ymax],'color','k','linewidth',0.75);
plot([-0.5 -0.5],[-3 ymax],'color','b','linewidth',0.75);
axis([-2 4 -3 ymax]);
title(['Reward response, all neurons (n=' num2str(sum(Sel)) ')']);
legend([a b c],'Suc trials','Mal trials','Wat trials');
text(-1.3,ymax-0.5,'PE','color','b');
xlabel('Seconds from reward delivery');


disp('Interspersed event responses')


%% categories and linear models
trialsback=10;

%neuron categories
masks=cat(2,bm_RD.mask_base',bm_RD.mask_curr',bm_RD.mask_mean');
%masks=cat(2,(bm_RD.mask_base' & bm_RD_BIC.mask_base'==0),bm_RD.mask_curr',bm_RD.mask_mean');

stacked_data = zeros(2,3);
for subset=1:3
    stacked_data(1,subset)=sum(masks(:,subset));
end

subplot(3,4,4)
hold on;
b = bar(stacked_data./sum(stacked_data,2),'stacked','edgecolor','none');
b(3).FaceColor = 'Flat';
b(3).CData = Colors('mean');
b(2).FaceColor = 'Flat';
b(2).CData = Colors('current');
b(1).FaceColor = 'Flat';
b(1).CData = Colors('rpe');
axis([0 7 0 1])
text(1.5,0.8,'Mean','color',Colors('mean'));
text(1.5,0.5,'Current','color',Colors('current'));
text(1.5,0.2,'RPE','color',Colors('rpe'));
set(gca,'xtick',[]);
yticks([0,.2,.4,.6,.8,1])
ylabel('Fraction of population')



linecolor{1}=Colors('rpe');
linecolor{2}=Colors('current');
linecolor{3}=Colors('mean');

for i=1:3
    if i==1 Sel=masks(:,1); end %base neurons
    if i==2 Sel=bm_RD.mask_curr'; end
    if i==3 Sel=bm_RD.mask_mean'; end

    %calculate model based on all reward selective neurons
    allRDHz=[];
    allpreds=[];
    selindex=find(Sel);
    for j=1:sum(Sel)
        allRDHz=cat(1,allRDHz,RewHist.RDHz{selindex(j),1});
        allpreds=cat(1,allpreds,RewHist.Predictors{selindex(j),1}(:,1:trialsback+1));
    end
    
    if i==2 allpreds(allpreds==0.75)=0.82; end %set to 0.75, which is best for base neurons, but 0.82 is best for current only
    allVPmdl=fitlm(allpreds,allRDHz);
    stats.trials.subsets{i,1}=allVPmdl.Coefficients.pValue(2:12);

    %plot coefficients
    subplot(3,4,8);
    hold on;
    errorbar(0:trialsback,allVPmdl.Coefficients.Estimate(2:trialsback+2),allVPmdl.Coefficients.SE(2:trialsback+2),'color',linecolor{i},'linewidth',1.5);
    

end
plot([-1 trialsback+1],[0 0],':','color','k','linewidth',0.75);
xlabel('Trials back');
ylabel('Coefficient weight');
title('Trial history regression, subsets');
legend('RPE','Current','Mean');
axis([-1 trialsback+1 -1.2 7]);

disp('B&G regressions')

%% RPE neurons based on different RPE values

bins=[-3 -1.2;-1.2 -0.7;-0.7 -0.3;-0.3 0;0 0.3;0.3 0.7;0.7 1.2;1.2 3];

Sel=masks(:,1);
mintrials=1;

global Dura Tm BSIZE
BSIZE=0.01; %Do not change
Dura=[-4 10]; Tm=Dura(1):BSIZE:Dura(2); %should match whatever was used to create R for simplicity
Baseline=[-11 -1]; %for PSTHs, should match whatever is in R

%smoothing filter for whole PSTH
PSTHsmoothbins=50; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',20); %std=3.98
PSTHfilterweights=pdf(halfnormal,0:PSTHsmoothbins);

%smoothing filter for individual trials
trialsmoothbins=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',10); %std=3.98
trialfilterweights=pdf(halfnormal,0:trialsmoothbins);

sessions=unique(R.Ninfo(:,1));
index=find(Sel);
RPEs=[];
for i=1:sum(Sel)
    
    %collect all RPE values for all cells of interest
    RPEs=cat(1,RPEs,os(index(i)).mod_RD.base.RPEs);
    
    session=R.Ninfo(index(i),1); %what's the name of session
    sessnum=ismember(sessions,session); %which session of all sessions
    neuronnum=ismember(R.Ninfo(strcmp(session,R.Ninfo(:,1)),2),R.Ninfo(index(i),2)); %which number neuron in the session
    RD=strmatch('RD',RAW(sessnum).Einfo(:,2),'exact'); %get index for reward delivery times
    Cue=strmatch('Cue',RAW(sessnum).Einfo(:,2),'exact'); %get index for reward delivery times
    
    %get mean baseline firing for all reward trials
    [Bcell1,B1n]=MakePSR04(RAW(sessnum).Nrast(neuronnum),RAW(sessnum).Erast{Cue},Baseline,{2});% makes trial by trial rasters for baseline
    basespk=[];
    for y= 1:B1n
        basespk(1,y)=sum(Bcell1{1,y}>Baseline(1));
    end

    FiringRate=basespk/(Baseline(1,2)-Baseline(1,1));
    Bmean=nanmean(FiringRate);
    Bstd=nanstd(FiringRate);   
    
    %make PSTHs around reward delivery
    trialtimes=RAW(sessnum).Erast{RD,1}(:,1); %reward delivery times
    for j=1:length(bins)
        RPEs=os(index(i)).mod_RD.base.RPEs;
        RPEsZ = (RPEs - mean(RPEs))/std(RPEs);
        toi=trialtimes(RPEsZ>=bins(j,1) & RPEsZ<=bins(j,2)); %which trials, binned by model-produced RPE
        if length(toi)>=mintrials
            %make PSTHs
            [PSR1,~]=MakePSR04(RAW(sessnum).Nrast(neuronnum),toi,Dura,{2});% makes collpased rasters. PSR1 is a cell(neurons)
            smoothedtrials=[];
            for trial=1:length(PSR1)
                binned=histcounts(PSR1{trial},[Tm Tm(end)+(Tm(end)-Tm(end-1))]);
                for l=1:length(Tm)
                    smoothedtrials(trial,l)=sum(binned(1,l-min([l-1 trialsmoothbins]):l).*fliplr(trialfilterweights(1:min([l trialsmoothbins+1]))))/sum(trialfilterweights(1:min([l trialsmoothbins+1])));
                end
            end
            PTH1=mean(smoothedtrials,1)/BSIZE;
            PTH1smooth=[];
            for l=1:length(Tm)
                PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 PSTHsmoothbins]):l).*fliplr(PSTHfilterweights(1:min([l PSTHsmoothbins+1]))))/sum(PSTHfilterweights(1:min([l PSTHsmoothbins+1])));
            end
            RPEPSTHs{j,1}.PSTHraw(i,:)=PTH1smooth;

            %normalize already smoothed activity
            for l=1:length(PTH1smooth)
                RPEPSTHs{j,1}.PSTHz(i,l)=(PTH1smooth(l)-Bmean)/Bstd;
            end
        else
            RPEPSTHs{j,1}.PSTHraw(i,1:length(Tm))=NaN;
            RPEPSTHs{j,1}.PSTHz(i,1:length(Tm))=NaN;
        end
    end
end

%plotting the histograms
ymax=10;
Xaxis=[-1.5 4];
Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
time1=R.Param.Tm(Ishow);

alph=0.2;

psthcolor=[0 0.9 1;...
           0 0.7 0.75;...
           0 0.5 0.5;...
           0 0.3 0.25;...
           0.25 0.3 0;...
           0.5 0.5 0;...
           0.75 0.7 0;...
           1 0.9 0];


subplot(3,4,6:7);
hold on;
patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-4 ymax ymax -4],Colors('gray'),'edgecolor','none');alpha(alph);
a=[];
for i=1:length(RPEPSTHs) %for each bin
    
    %plot activity for each bin
    psth1=nanmean(RPEPSTHs{i,1}.PSTHz(:,Ishow),1); 
    sem1=nanste(RPEPSTHs{i,1}.PSTHz(:,Ishow),1); %calculate standard error of the mean
    up1=psth1+sem1;
    down1=psth1-sem1;

    %make the plot
    a{i}=plot(time1,psth1,'Color',psthcolor(i,:),'linewidth',1.5);
    patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],psthcolor(i,:),'EdgeColor','none');alpha(alph);
    labels{i}=[num2str(bins(i,1)) ' < RPE < ' num2str(bins(i,2))];
end
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
plot([0 0],[-5 ymax],'color','k','linewidth',0.5);
axis([Xaxis(1) Xaxis(2) -4 ymax]);
ylabel('Mean firing (z-score)');
xlabel('Seconds from RD');
colormap(psthcolor);
h=colorbar;
set(h,'YTick',[0:0.125:1]);
set(h,'yticklabel',[-3;-1.2;-.7;-.3;0;0.3;0.7;1.2;3])
set(get(h,'title'),'string','RPE');
%legend([a{:}],labels,'location','northeast');

if plot_sub_figures

    %% cue responses

    masks=cat(2,bm_cue.mask_base',bm_cue.mask_mean');


    figure;


    %% bar graphs
    stacked_data = zeros(2,3);

    for subset=1:2
        stacked_data(1,subset)=sum(masks(:,subset));
    end

    subplot(2,3,1)
    hold on;
    b = bar(stacked_data./sum(stacked_data,2),'stacked','edgecolor','none');
    b(2).FaceColor = 'Flat';
    b(2).CData = Colors('mean');
    b(1).FaceColor = 'Flat';
    b(1).CData = Colors('rpe');
    axis([0 5 0 1])
    set(gca,'xtick',[]);
    yticks([0,.2,.4,.6,.8,1])
    ylabel('Fraction of population')
    legend('Value','Unmod.','location','northwest');


    %% activity plots for subsets
    ymax=9;
    %ymax=5;
    Xaxis=[-1 2];
    Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
    time1=R.Param.Tm(Ishow);

    suc_cue=strcmp('CueP1', R.Erefnames);
    mal_cue=strcmp('CueP2', R.Erefnames);
    wat_cue=strcmp('CueP3', R.Erefnames);

    titles={'Value neurons (n=';'Unmodulated neurons (n='};
    for i=1:2
        if i==1 Sel=bm_cue.mask_base'; end %base neurons
        if i==2 Sel=bm_cue.mask_mean'; end %mean neurons


        %plot cue after suc
        psth1=nanmean(R.Ev(suc_cue).PSTHz(Sel,Ishow),1); 
        sem1=nanste(R.Ev(suc_cue).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
        up1=psth1+sem1;
        down1=psth1-sem1;

        %plot cue after malt
        psth2=nanmean(R.Ev(mal_cue).PSTHz(Sel,Ishow),1); 
        sem2=nanste(R.Ev(mal_cue).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
        up2=psth2+sem2;
        down2=psth2-sem2;

        %plot cue after water
        psth3=nanmean(R.Ev(wat_cue).PSTHz(Sel,Ishow),1); 
        sem3=nanste(R.Ev(wat_cue).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
        up3=psth3+sem3;
        down3=psth3-sem3;

        %make the plot
        subplot(2,3,1+i);
        hold on;
        title([titles{i} num2str(sum(Sel)) ')'])
        patch([CueWindow(1) CueWindow(1) CueWindow(2) CueWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
        a=plot(time1,psth2,'Color',Colors('maltodextrin'),'linewidth',1);
        b=plot(time1,psth1,'Color',Colors('sucrose'),'linewidth',1);
        c=plot(time1,psth3,'Color',Colors('water'),'linewidth',1);

        patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(alph);
        patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(alph);
        patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('water'),'EdgeColor','none');alpha(alph);

        plot([-2 5],[0 0],':','color','k','linewidth',0.75);
        plot([0 0],[-2 ymax],'color','k','linewidth',0.5);
        axis([Xaxis(1) Xaxis(2) -2 ymax]);
        ylabel('Mean firing (z-score)');
        xlabel('Seconds from cue');
        if i==2 legend([b a c],'After suc','After mal','After wat','location','northeast'); end




    end
    disp('Subsets event responses')

    %% B&G for cue subsets

    linecolor{1}=Colors('rpe');
    linecolor{2}=Colors('mean');

    for i=1:2
        if i==1 Sel=bm_cue.mask_base'; end %base neurons
        if i==2 Sel=bm_cue.mask_mean'; end %base neurons

        %calculate model based on all reward selective neurons
        allRDHz=[];
        allpreds=[];
        selindex=find(Sel);
        for j=1:sum(Sel)
            allRDHz=cat(1,allRDHz,RewHist.CueHz{selindex(j),1});
            allpreds=cat(1,allpreds,RewHist.Predictors{selindex(j),1}(:,2:trialsback+1));
        end
        
        %allpreds(allpreds==0.75)=0.5;

        allVPmdl=fitlm(allpreds,allRDHz);
        stats.cue_trials.subsets{i,1}=allVPmdl.Coefficients.pValue(2:11);

        %plot coefficients
        subplot(2,3,4);
        hold on;
        errorbar(1:trialsback,allVPmdl.Coefficients.Estimate(2:trialsback+1),allVPmdl.Coefficients.SE(2:trialsback+1),'color',linecolor{i},'linewidth',1.5);


    end
    plot([-1 trialsback+1],[0 0],':','color','k','linewidth',0.75);
    xlabel('Trials back');
    ylabel('Coefficient weight');
    title('Trial history regression, subsets');
    legend('Value','Unmod.');
    axis([0 trialsback+1 -2.5 4]);

    disp('B&G regressions')

    %% Value neurons based on different values

    bins=[-3 -1.2;-1.2 -0.5;-0.5 0.5;0.5 1.2;1.2 3];

    Sel=bm_cue.mask_base';
    mintrials=1;

    global Dura Tm BSIZE
    BSIZE=0.01; %Do not change
    Dura=[-4 10]; Tm=Dura(1):BSIZE:Dura(2); %should match whatever was used to create R for simplicity
    Baseline=[-11 -1]; %for PSTHs, should match whatever is in R

    %smoothing filter for whole PSTH
    PSTHsmoothbins=50; %number of previous bins used to smooth
    halfnormal=makedist('HalfNormal','mu',0,'sigma',10); %std=3.98
    PSTHfilterweights=pdf(halfnormal,0:PSTHsmoothbins);

    %smoothing filter for individual trials
    trialsmoothbins=25; %number of previous bins used to smooth
    halfnormal=makedist('HalfNormal','mu',0,'sigma',5); %std=3.98
    trialfilterweights=pdf(halfnormal,0:trialsmoothbins);

    sessions=unique(R.Ninfo(:,1));
    index=find(Sel);
    RPEPSTHs={};
    for i=1:sum(Sel)
        session=R.Ninfo(index(i),1); %what's the name of session
        sessnum=ismember(sessions,session); %which session of all sessions
        neuronnum=ismember(R.Ninfo(strcmp(session,R.Ninfo(:,1)),2),R.Ninfo(index(i),2)); %which number neuron in the session
        RD=strmatch('RD',RAW(sessnum).Einfo(:,2),'exact'); %get index for reward delivery times
        Cue=strmatch('Cue',RAW(sessnum).Einfo(:,2),'exact'); %get index for reward delivery times

        %get mean baseline firing for all reward trials
        [Bcell1,B1n]=MakePSR04(RAW(sessnum).Nrast(neuronnum),RAW(sessnum).Erast{Cue},Baseline,{2});% makes trial by trial rasters for baseline
        basespk=[];
        for y= 1:B1n
            basespk(1,y)=sum(Bcell1{1,y}>Baseline(1));
        end

        FiringRate=basespk/(Baseline(1,2)-Baseline(1,1));
        Bmean=nanmean(FiringRate);
        Bstd=nanstd(FiringRate);   

        %make PSTHs around cue onset
        cuetimes=RAW(sessnum).Erast{Cue,1}(:,1); %cue times
        RDtimes=RAW(sessnum).Erast{RD,1}(:,1); %reward delivery times
        trialtimes=[];
        for trial = 1:length(RDtimes)
            trialtimes(trial,1) = max(cuetimes(cuetimes<RDtimes(trial)));
        end

        for j=1:length(bins)
            RPEs=os(index(i)).mod_cue.base.V;
            RPEsZ = (RPEs - mean(RPEs))/std(RPEs);
            toi=trialtimes(RPEsZ>=bins(j,1) & RPEsZ<=bins(j,2)); %which trials, binned by model-produced RPE
            if length(toi)>=mintrials
                %make PSTHs
                [PSR1,~]=MakePSR04(RAW(sessnum).Nrast(neuronnum),toi,Dura,{2});% makes collpased rasters. PSR1 is a cell(neurons)
                smoothedtrials=[];
                for trial=1:length(PSR1)
                    binned=histcounts(PSR1{trial},[Tm Tm(end)+(Tm(end)-Tm(end-1))]);
                    for l=1:length(Tm)
                        smoothedtrials(trial,l)=sum(binned(1,l-min([l-1 trialsmoothbins]):l).*fliplr(trialfilterweights(1:min([l trialsmoothbins+1]))))/sum(trialfilterweights(1:min([l trialsmoothbins+1])));
                    end
                end
                PTH1=mean(smoothedtrials,1)/BSIZE;
                PTH1smooth=[];
                for l=1:length(Tm)
                    PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 PSTHsmoothbins]):l).*fliplr(PSTHfilterweights(1:min([l PSTHsmoothbins+1]))))/sum(PSTHfilterweights(1:min([l PSTHsmoothbins+1])));
                end
                RPEPSTHs{j,1}.PSTHraw(i,:)=PTH1smooth;

                %normalize already smoothed activity
                for l=1:length(PTH1smooth)
                    RPEPSTHs{j,1}.PSTHz(i,l)=(PTH1smooth(l)-Bmean)/Bstd;
                end
            else
                RPEPSTHs{j,1}.PSTHraw(i,1:length(Tm))=NaN;
                RPEPSTHs{j,1}.PSTHz(i,1:length(Tm))=NaN;
            end
        end
    end

    %plotting the histograms
    ymax=9;
    Xaxis=[-1 2];
    Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
    time1=R.Param.Tm(Ishow);


    alph=0.2;

    psthcolor=[0 0.1 0;...
               0 0.3 0.2;...
               0 0.5 0.3;...
               0 0.75 0.5;...
               0 1 0.7];

    subplot(2,2,4);
    hold on;
    patch([CueWindow(1) CueWindow(1) CueWindow(2) CueWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
    a={};
    for i=1:length(RPEPSTHs) %for each bin

        %plot activity for each bin
        psth1=nanmean(RPEPSTHs{i,1}.PSTHz(:,Ishow),1); 
        sem1=nanste(RPEPSTHs{i,1}.PSTHz(:,Ishow),1); %calculate standard error of the mean
        up1=psth1+sem1;
        down1=psth1-sem1;

        %make the plot
        a{i}=plot(time1,psth1,'Color',psthcolor(i,:),'linewidth',1.5);
        patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],psthcolor(i,:),'EdgeColor','none');alpha(alph);
        labels{i}=[num2str(bins(i,1)) ' < RPE < ' num2str(bins(i,2))];
    end
    plot([-2 5],[0 0],':','color','k','linewidth',0.75);
    plot([0 0],[-3 ymax],'color','k','linewidth',0.5);
    axis([Xaxis(1) Xaxis(2) -2 ymax]);
    ylabel('Mean firing (z-score)');
    xlabel('Seconds from cue');
    colormap(psthcolor);
    h=colorbar;
    set(h,'YTick',[0:0.2:1]);
    set(h,'yticklabel',[-3;-1.2;-.5;0.5;1.2;3])
    set(get(h,'title'),'string','Value');
    %legend([a{:}],labels,'location','northeast');
end