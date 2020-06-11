function stats=plot_int(Colors, R, RewHist, bm_RD, RAW, os, bm_cue, plot_value, bm_RD_BIC)

color = {Colors('sucrose'),Colors('maltodextrin'),Colors('total'),Colors('blockrose'),Colors('blockrodextrin'),Colors('total')};

xaxis=[-2 5];
time_indicies=find(R.Param.Tm>=xaxis(1) & R.Param.Tm<=xaxis(end)); %returns the indicies where these are both true
activity_xaxis=R.Param.Tm(time_indicies);
suc_delivery=strcmp('RD1', R.Erefnames);
malt_delivery=strcmp('RD2', R.Erefnames);
suc_cue=strcmp('CueP1', R.Erefnames);
mal_cue=strcmp('CueP2', R.Erefnames);
RD1P1=strcmp('RD1P1', R.Erefnames);
RD1P2=strcmp('RD1P2', R.Erefnames);
RD2P1=strcmp('RD2P1', R.Erefnames);
RD2P2=strcmp('RD2P2', R.Erefnames);
%masks=cat(2,(bm_RD.mask_base' & bm_RD_BIC.mask_base'==0),bm_RD.mask_curr',bm_RD.mask_mean');
masks=cat(2,bm_RD.mask_base',bm_RD.mask_curr',bm_RD.mask_mean');
%masks=cat(2,bm_RD_BIC.mask_base',bm_RD_BIC.mask_curr',bm_RD_BIC.mask_mean');

region = strcmp(R.Region,'VP');
task = R.Blocks==0;

HistWindow=RewHist.RDWindow;
CueWindow=RewHist.CueWindow;

figure;

%% activity plots
ymax=4;
alph=0.3;
Sel=region&task;

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

subplot(3,2,1)
hold on
plot(activity_xaxis,psth1,'Color',Colors('sucrose'),'linewidth',1);
plot(activity_xaxis,psth2,'Color',Colors('maltodextrin'),'linewidth',1);
patch([activity_xaxis,activity_xaxis(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(0.5);
patch([activity_xaxis,activity_xaxis(end:-1:1)],[up2,down2(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(0.5);
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
plot([0 0],[-2 ymax],'color','k','linewidth',0.5);
plot([-0.5 -0.5],[-2 ymax],'color','b','linewidth',0.5);
axis([-2 4 -1 ymax]);
title(['Reward response, all neurons (n=' num2str(sum(Sel)) ')']);
legend('Suc trials','Mal trials');
text(-1.3,ymax-0.5,'PE','color','b');
xlabel('Seconds from reward delivery');

%broken down by previous trial
ymax=5;
Xaxis=[0 2.5];
Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
time1=R.Param.Tm(Ishow);

%plot suc after suc
psth1=nanmean(R.Ev(RD1P1).PSTHz(Sel,Ishow),1);
sem1=nanste(R.Ev(RD1P1).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up1=psth1+sem1;
down1=psth1-sem1;

%plot suc after malt
psth2=nanmean(R.Ev(RD1P2).PSTHz(Sel,Ishow),1);
sem2=nanste(R.Ev(RD1P2).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up2=psth2+sem2;
down2=psth2-sem2;

%plot malt after suc
psth3=nanmean(R.Ev(RD2P1).PSTHz(Sel,Ishow),1);
sem3=nanste(R.Ev(RD2P1).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up3=psth3+sem3;
down3=psth3-sem3;

%plot malt after malt
psth4=nanmean(R.Ev(RD2P2).PSTHz(Sel,Ishow),1);
sem4=nanste(R.Ev(RD2P2).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up4=psth4+sem4;
down4=psth4-sem4;

%make the plot
subplot(3,3,3);
hold on;
title(['Current/previous reward'])
patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-5 10 10 -5],Colors('gray'),'EdgeColor','none');
a=plot(time1,psth2,'Color',Colors('sucrosemnp'),'linewidth',1);
b= plot(time1,psth1,'Color',Colors('sucrosesnp'),'linewidth',1);
c=plot(time1,psth4,'Color',Colors('maltodextrinmnp'),'linewidth',1);
d=plot(time1,psth3,'Color',Colors('maltodextrinsnp'),'linewidth',1);
patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('sucrosemnp'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrosesnp'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('maltodextrinsnp'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up4,down4(end:-1:1)],Colors('maltodextrinmnp'),'EdgeColor','none');alpha(alph);
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
%plot([0 0],[-2 ymax],'color','k','linewidth',0.5);
% plot([HistWindow(1) HistWindow(1)],[-5 10],'--','color','k','linewidth',0.85);
% plot([HistWindow(2) HistWindow(2)],[-5 10],'--','color','k','linewidth',0.85);
%text(0.1,ymax-0.5,'RD','color','k');
axis([Xaxis(1) Xaxis(2) -1 ymax]);
ylabel('Mean firing (z-score)');
xlabel('Seconds from RD');
legend([a b c d],'Suc after mal','Suc after suc','Mal after mal','Mal after suc','location','northeast');


disp('Interspersed event responses')


%% fixed window linear models
trialsback=10;

%all neurons
Sel=region>=0; %allneurons
Sel=Sel(region&task); %only the neurons from VP interspersed

%calculate model based on all reward selective neurons
allRDHz=[];
allpreds=[];
selindex=find(Sel);
for i=1:sum(Sel)
    allRDHz=cat(1,allRDHz,RewHist.RDHz{selindex(i),1});
    allpreds=cat(1,allpreds,RewHist.Predictors{selindex(i),1}(:,1:trialsback+1));
end


allVPmdl=fitlm(allpreds,allRDHz);
stats.trials.all=allVPmdl.Coefficients.pValue(2:12);

%plot coefficients
subplot(3,3,4);
hold on;
errorbar(0:trialsback,allVPmdl.Coefficients.Estimate(2:trialsback+2),allVPmdl.Coefficients.SE(2:trialsback+2),'linewidth',1.5);
axis([-1 trialsback+1 -1 2]);
plot([-1 trialsback+1],[0 0],':','color','k','linewidth',0.75);
xlabel('Trials back');
ylabel('Coefficient weight');
title('Trial history regression, all neurons');

linecolor{1}=Colors('rpe');
linecolor{2}=Colors('current');
linecolor{3}=Colors('mean');

for i=1:3
    Sel=masks(:,i);
    Sel=Sel(region&task); %only the neurons from VP interspersed
    
    %calculate model based on all reward selective neurons
    allRDHz=[];
    allpreds=[];
    selindex=find(Sel);
    for j=1:sum(Sel)
        allRDHz=cat(1,allRDHz,RewHist.RDHz{selindex(j),1});
        allpreds=cat(1,allpreds,RewHist.Predictors{selindex(j),1}(:,1:trialsback+1));
    end
    
    allVPmdl=fitlm(allpreds,allRDHz);
    stats.trials.subsets{i,1}=allVPmdl.Coefficients.pValue(2:12);
    
    %plot coefficients
    subplot(3,3,6);
    hold on;
    errorbar(0:trialsback,allVPmdl.Coefficients.Estimate(2:trialsback+2),allVPmdl.Coefficients.SE(2:trialsback+2),'color',linecolor{i},'linewidth',1.5);
    
    
end
plot([-1 trialsback+1],[0 0],':','color','k','linewidth',0.75);
xlabel('Trials back');
ylabel('Coefficient weight');
title('Trial history regression, subsets');
legend('RPE','Current','Unmod.');
axis([-1 trialsback+1 -3 5]); %AIC
%axis([-1 trialsback+1 -3 6]); %BIC

disp('B&G regressions')


%% bar graphs
stacked_data = zeros(2,3);

for subset=1:3
    stacked_data(1,subset)=sum(region&masks(:,subset)&task);
end

subplot(3,3,5)
hold on;
b = bar(stacked_data./sum(stacked_data,2),'stacked','edgecolor','none');
b(3).FaceColor = 'Flat';
b(3).CData = Colors('mean');
b(2).FaceColor = 'Flat';
b(2).CData = Colors('current');
b(1).FaceColor = 'Flat';
b(1).CData = Colors('rpe');
axis([0 5 0 1])
% text(1.5,0.7,'Mean','color',[0.6 0.6 0.6]);
% text(1.5,0.3,'Current only','color',[0.6 0 0.6]);
% text(1.5,0.1,'RPE','color',[0 0.8 0.6]);
set(gca,'xtick',[]);
yticks([0,.2,.4,.6,.8,1])
ylabel('Fraction of population')
legend('RPE','Current','Unmod.','location','northwest');


%% activity plots for subsets
ymax=11;
%ymax=5;
Xaxis=[-0.5 3];
Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
time1=R.Param.Tm(Ishow);

titles={'RPE neurons (n=';'Current outcome neurons (n=';'Unmodulated neurons (n='};
for i=1:3
    if i==1 Sel=masks(:,1); end %base neurons
    if i==2 Sel=masks(:,2); end %current outcome neurons
    if i==3 Sel=masks(:,3); end %mean neurons
    Sel=Sel&region&task; %only the neurons from VP interspersed
    
    
    
    %plot suc after suc
    psth1=nanmean(R.Ev(RD1P1).PSTHz(Sel,Ishow),1);
    sem1=nanste(R.Ev(RD1P1).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
    up1=psth1+sem1;
    down1=psth1-sem1;
    
    %plot suc after malt
    psth2=nanmean(R.Ev(RD1P2).PSTHz(Sel,Ishow),1);
    sem2=nanste(R.Ev(RD1P2).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
    up2=psth2+sem2;
    down2=psth2-sem2;
    
    %plot malt after suc
    psth3=nanmean(R.Ev(RD2P1).PSTHz(Sel,Ishow),1);
    sem3=nanste(R.Ev(RD2P1).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
    up3=psth3+sem3;
    down3=psth3-sem3;
    
    %plot malt after malt
    psth4=nanmean(R.Ev(RD2P2).PSTHz(Sel,Ishow),1);
    sem4=nanste(R.Ev(RD2P2).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
    up4=psth4+sem4;
    down4=psth4-sem4;
    
    %make the plot
    subplot(3,3,6+i);
    hold on;
    title([titles{i} num2str(sum(Sel)) ')'])
    patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
    a=plot(time1,psth2,'Color',Colors('sucrosemnp'),'linewidth',1);
    b=plot(time1,psth1,'Color',Colors('sucrosesnp'),'linewidth',1);
    c=plot(time1,psth4,'Color',Colors('maltodextrinmnp'),'linewidth',1);
    d=plot(time1,psth3,'Color',Colors('maltodextrinsnp'),'linewidth',1);
    patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('sucrosemnp'),'EdgeColor','none');alpha(alph);
    patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrosesnp'),'EdgeColor','none');alpha(alph);
    patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('maltodextrinsnp'),'EdgeColor','none');alpha(alph);
    patch([time1,time1(end:-1:1)],[up4,down4(end:-1:1)],Colors('maltodextrinmnp'),'EdgeColor','none');alpha(alph);
    plot([-2 5],[0 0],':','color','k','linewidth',0.75);
    plot([0 0],[-2 ymax],'color','k','linewidth',0.5);
    axis([Xaxis(1) Xaxis(2) -2 ymax]);
    ylabel('Mean firing (z-score)');
    xlabel('Seconds from RD');
    if i==3 legend([a b c d],'Suc after mal','Suc after suc','Mal after mal','Mal after suc','location','northeast'); end
    
    
    
    
end
disp('Subsets event responses')

%% evaluation of RPE models
figure;
alph=0.3;

%raster
window=[0 3];
sample_neuron_number=795; %794-798, 829, 1101

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
xlabel('Seconds from reward delivery');
xlim(window);

disp('Raster')

%% RPE neurons based on different RPE values

plot_RPEs=1;
if plot_RPEs
    
    bins=[-3 -1.2;-1.2 -0.7;-0.7 -0.3;-0.3 0;0 0.3;0.3 0.7;0.7 1.2;1.2 3];
    
    Sel=masks(:,1) & region & task;
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
    ymax=10; %AIC
    %ymax=12; %BIC
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
    
    subplot(2,4,2:3);
    hold on;
    patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
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
    xlabel('Seconds from RD');
    colormap(psthcolor);
    h=colorbar;
    set(h,'YTick',[0:0.125:1]);
    set(h,'yticklabel',[-3;-1.2;-.7;-.3;0;0.3;0.7;1.2;3])
    set(get(h,'title'),'string','RPE');
    %legend([a{:}],labels,'location','northeast');
end

%% bar graphs comparing regions
stacked_data = zeros(2,3);

for reg=1:2
    if reg==1 regsel = strcmp(R.Region,'VP'); end
    if reg==2 regsel = strcmp(R.Region,'NA'); end
    
    for subset=1:3
        stacked_data(reg,subset)=sum(regsel&masks(:,subset)&task);
    end
end

subplot(2,4,7)
hold on;
b = bar(stacked_data./sum(stacked_data,2),'stacked','edgecolor','none');
b(3).FaceColor = 'Flat';
b(3).CData = Colors('mean');
b(2).FaceColor = 'Flat';
b(2).CData = Colors('current');
b(1).FaceColor = 'Flat';
b(1).CData = Colors('rpe');
axis([0 4 0 1])
% text(1.5,0.7,'Mean','color',[0.6 0.6 0.6]);
% text(1.5,0.3,'Current only','color',[0.6 0 0.6]);
% text(1.5,0.1,'RPE','color',[0 0.8 0.6]);
set(gca,'xtick',[]);
yticks([0,.2,.4,.6,.8,1])
ylabel('Fraction of population')
legend('RPE','Current','Unmod.','location','northwest');

if plot_value
    
    %% cue responses
    
    masks=cat(2,bm_cue.mask_base',bm_cue.mask_mean');
    
    
    figure;
    
    
    %% bar graphs
    stacked_data = zeros(2,3);
    
    for subset=1:2
        stacked_data(1,subset)=sum(region&masks(:,subset)&task);
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
    ymax=10;
    %ymax=5;
    Xaxis=[-1 2];
    Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
    time1=R.Param.Tm(Ishow);
    
    titles={'Value neurons (n=';'Unmodulated neurons (n='};
    for i=1:2
        if i==1 Sel=bm_cue.mask_base'; end %base neurons
        if i==2 Sel=bm_cue.mask_mean'; end %mean neurons
        Sel=Sel&region&task; %only the neurons from VP interspersed
        
        
        
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
        
        %make the plot
        subplot(2,3,1+i);
        hold on;
        title([titles{i} num2str(sum(Sel)) ')'])
        patch([CueWindow(1) CueWindow(1) CueWindow(2) CueWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
        a=plot(time1,psth2,'Color',Colors('maltodextrin'),'linewidth',1);
        b=plot(time1,psth1,'Color',Colors('sucrose'),'linewidth',1);
        patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(alph);
        patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(alph);
        plot([-2 5],[0 0],':','color','k','linewidth',0.75);
        plot([0 0],[-2 ymax],'color','k','linewidth',0.5);
        axis([Xaxis(1) Xaxis(2) -2 ymax]);
        ylabel('Mean firing (z-score)');
        xlabel('Seconds from cue');
        if i==2 legend([b a],'After suc','After mal','location','northeast'); end
        
        
        
        
    end
    disp('Subsets event responses')
    
    %% B&G for cue subsets
    
    linecolor{1}=Colors('rpe');
    linecolor{2}=Colors('mean');
    
    for i=1:2
        if i==1 Sel=bm_cue.mask_base'; end %base neurons
        if i==2 Sel=bm_cue.mask_mean'; end %base neurons
        Sel=Sel(region&task); %only the neurons from VP interspersed
        
        %calculate model based on all reward selective neurons
        allRDHz=[];
        allpreds=[];
        selindex=find(Sel);
        for j=1:sum(Sel)
            allRDHz=cat(1,allRDHz,RewHist.CueHz{selindex(j),1});
            allpreds=cat(1,allpreds,RewHist.Predictors{selindex(j),1}(:,2:trialsback+1));
        end
        
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
    axis([0 trialsback+1 -1 3]);
    
    disp('B&G regressions')
    
    %% Value neurons based on different values
    
    bins=[-3 -1.2;-1.2 -0.5;-0.5 0.5;0.5 1.2;1.2 3];
    
    Sel=bm_cue.mask_base' & region & task;
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