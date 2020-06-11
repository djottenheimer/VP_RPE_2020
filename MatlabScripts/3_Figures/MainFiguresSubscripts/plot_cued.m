function Stats=plot_cued(Colors, R, RewHist, bm_RD, DOI, ROI, os, bm_cue, RAW)

%which neurons are plotted
IncludedNeurons=ChooseNs(R,DOI,ROI);

%subselect sessions of interest
for i=1:length(RAW)
    name=char(RAW(i).Ninfo(1,1));
    sessrat(i,1)=str2num(name(2:3));
    sessday(i,1)=str2num(name(8:9));
end

%decide which sessions to include
inc_sessions = (ismember(sessrat,ROI) & ismember(sessday,DOI));
RAW=RAW(inc_sessions);


xaxis=[-2 5];

%event names
RD1NP=strcmp('Cue3RD1', R.Erefnames);
RD1PR=strcmp('Cue1RD', R.Erefnames);
RD2NP=strcmp('Cue3RD2', R.Erefnames);
RD2PR=strcmp('Cue2RD', R.Erefnames);
Cue1=strcmp('Cue1', R.Erefnames);
Cue2=strcmp('Cue2', R.Erefnames);
Cue31=strcmp('Cue31', R.Erefnames);
Cue32=strcmp('Cue32', R.Erefnames);

%when the model data is only on subselected neurons, make the logicals in terms of all neurons           
NNumber=[1:length(R.Ninfo)]';
index=NNumber(IncludedNeurons);
histmasks=zeros(length(R.Ninfo),3);
histmasks(index,1)=bm_RD.mask_base'|bm_RD.mask_base_cue';
histmasks(index,2)=bm_RD.mask_curr'|bm_RD.mask_curr_cue';
histmasks(index,3)=bm_RD.mask_mean'|bm_RD.mask_mean_cue';

cuemasks=zeros(length(R.Ninfo),2);
cuemasks(index,1)=bm_RD.mask_base_cue'|bm_RD.mask_curr_cue'|bm_RD.mask_mean_cue';
cuemasks(index,2)=bm_RD.mask_base'|bm_RD.mask_curr'|bm_RD.mask_mean';

%analysis windows
HistWindow = [0.75 1.95];
CueWindow = [0 0.75];

figure;

%% bar graph for history effects
stacked_data = zeros(2,3);


for subset=1:3
    stacked_data(1,subset)=sum(IncludedNeurons&histmasks(:,subset)&cuemasks(:,2));
    stacked_data(2,subset)=sum(IncludedNeurons&histmasks(:,subset)&cuemasks(:,1));
end

subplot(3,6,1)
hold on;
b = bar([1 2],stacked_data./sum(IncludedNeurons),'stacked','edgecolor','none');
b(3).FaceColor = 'Flat';
b(3).CData = Colors('mean');
b(2).FaceColor = 'Flat';
b(2).CData = Colors('current');
b(1).FaceColor = 'Flat';
b(1).CData = Colors('rpe');
axis([0 7 0 1])
legend('RPE','Current','Unmod.');
set(gca,'xtick',[]);
yticks([0,.2,.4,.6,.8,1])
ylabel('Fraction of population')




%chi2 for distribution of cue effects among history effects
for i = 1:length(histmasks)
    histcategory(i,1)=sum(find(histmasks(i,:)));
    cuecategory(i,1)=sum(find(cuemasks(i,:)));
end

[table,Stats.HistCueDistX2,Stats.HistCueDistpval]=crosstab(histcategory(IncludedNeurons),cuecategory(IncludedNeurons));
%% mean activity
Sel=IncludedNeurons;
ymax=5;
alph=0.3;

Xaxis=[-2 4];
Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
time1=R.Param.Tm(Ishow);

%plot suc after non-specific cue
psth1=nanmean(R.Ev(RD1NP).PSTHz(Sel,Ishow),1); 
sem1=nanste(R.Ev(RD1NP).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up1=psth1+sem1;
down1=psth1-sem1;

%plot suc after sucrose cue
psth2=nanmean(R.Ev(RD1PR).PSTHz(Sel,Ishow),1); 
sem2=nanste(R.Ev(RD1PR).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up2=psth2+sem2;
down2=psth2-sem2;

%plot malt after non-specific cue
psth3=nanmean(R.Ev(RD2NP).PSTHz(Sel,Ishow),1); 
sem3=nanste(R.Ev(RD2NP).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up3=psth3+sem3;
down3=psth3-sem3;

%plot malt after maltodextrin cue
psth4=nanmean(R.Ev(RD2PR).PSTHz(Sel,Ishow),1); 
sem4=nanste(R.Ev(RD2PR).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up4=psth4+sem4;
down4=psth4-sem4;

%make the plot
subplot(3,3,4);
hold on;
title(['All neurons (n=' num2str(sum(Sel)) ')']);
patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
a=plot(time1,psth1,'Color',Colors('sucrose'),'linewidth',1);
b=plot(time1,psth2,'Color',Colors('blockrose'),'linewidth',1);
c=plot(time1,psth3,'Color',Colors('maltodextrin'),'linewidth',1);
d=plot(time1,psth4,'Color',Colors('blockrodextrin'),'linewidth',1);
patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('blockrose'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up4,down4(end:-1:1)],Colors('blockrodextrin'),'EdgeColor','none');alpha(alph);
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
plot([0 0],[-3 ymax],'color','k','linewidth',0.5);
text(0.1,ymax-0.5,'RD','color','k');
axis([Xaxis(1) Xaxis(2) -0.5 ymax]);
ylabel('Mean firing (z-score)');
xlabel('Seconds from RD');
legend([a b c d],'Non-pr suc','Pred suc','Non-pr mal','Pred mal','location','northeast');

disp('Mean activity');

%% reward history plots

%broken down by previous trial
ymax=10;
alph=0.3;
trialsback=10;
linecolor{1}=Colors('rpe');
linecolor{2}=Colors('current');
linecolor{3}=Colors('mean');

titles={'RPE neurons (n=';'Current outcome neurons (n=';'Mean neurons (n='};
for subset=1:3
    Sel=IncludedNeurons&histmasks(:,subset)&cuemasks(:,2);
    
    %calculate model based on all reward selective neurons
    allRDHz=[];
    allpreds=[];
    selindex=find(Sel);
    for i=1:sum(Sel)
        allRDHz=cat(1,allRDHz,RewHist.RDHz{selindex(i),1});
        allpreds=cat(1,allpreds,RewHist.Predictors{selindex(i),1}(:,1:trialsback+1));
    end

    allVPmdl=fitlm(allpreds,allRDHz);
    Stats.trials.subsets{subset,1}=allVPmdl.Coefficients.pValue(2:12);

    %plot coefficients
    subplot(3,6,2);
    hold on;
    errorbar(0:trialsback,allVPmdl.Coefficients.Estimate(2:trialsback+2),allVPmdl.Coefficients.SE(2:trialsback+2),'color',linecolor{subset},'linewidth',1.5);

    
end

axis([-1 trialsback+1 -1.5 4]);
plot([-1 trialsback+1],[0 0],':','color','k','linewidth',0.75);
xlabel('Trials back');
ylabel('Coefficient weight');
title('Trial history regression, subsets');
legend('RPE','Current','Mean');

disp('Rew History Responses')


%% scatterplot of cue values for effects on RD

subplot(3,6,5);
sz=25;
hold on;
sucvalues=nan(sum(IncludedNeurons),1);
malvalues=nan(sum(IncludedNeurons),1);

%mean
Sel=bm_RD.mask_mean_cue;
index=find(Sel);
vsuc = getParameters_ott(os, 'RD', 'mean_cue', 'vsuc');
vmal = getParameters_ott(os, 'RD', 'mean_cue', 'vmal');
scatter(-vsuc(Sel),-vmal(Sel),sz,Colors('mean'));
sucvalues(index)=-vsuc(Sel)';
malvalues(index)=-vmal(Sel)';

%current
Sel=bm_RD.mask_curr_cue;
index=find(Sel);
vsuc = getParameters_ott(os, 'RD', 'curr_cue', 'vsuc');
vmal = getParameters_ott(os, 'RD', 'curr_cue', 'vmal');
scatter(-vsuc(Sel),-vmal(Sel),sz,Colors('current'));
sucvalues(index)=-vsuc(Sel)';
malvalues(index)=-vmal(Sel)';

%base
Sel=bm_RD.mask_base_cue;
index=find(Sel);
vsuc = getParameters_ott(os, 'RD', 'base_cue', 'vsuc');
vmal = getParameters_ott(os, 'RD', 'base_cue', 'vmal');
scatter(-vsuc(Sel),-vmal(Sel),sz,Colors('rpe'));
sucvalues(index)=-vsuc(Sel)';
malvalues(index)=-vmal(Sel)';

axis([-1 1 -1 1]);
xlabel('Sucrose value');
ylabel('Maltodextrin value');
plot([0 0],[-1 1],'color','k');
plot([-1 1],[0 0],'color','k');


quadrant=[];
for i=1:length(sucvalues)
    if sucvalues(i)<0 & malvalues(i)>0
        quadrant(i)=1;
    elseif sucvalues(i)>0 & malvalues(i)>0
        quadrant(i)=2;
    elseif sucvalues(i)<0 & malvalues(i)<0
        quadrant(i)=3;
    elseif sucvalues(i)>0 & malvalues(i)<0
        quadrant(i)=4;
    end
end
quadrant=quadrant(find(quadrant));
text(-0.9,0.8,[num2str(round(100*sum(quadrant==1)/length(quadrant),1)) '%']);
text(0.6,0.8,[num2str(round(100*sum(quadrant==2)/length(quadrant),1)) '%']);
text(-0.9,-0.8,[num2str(round(100*sum(quadrant==3)/length(quadrant),1)) '%']);
text(0.6,-0.8,[num2str(round(100*sum(quadrant==4)/length(quadrant),1)) '%']);

[h,Stats.RDpval,Stats.RDX2] = chi2gof(quadrant,'Ctrs',[1 2 3 4],'expected',[length(quadrant)/4 length(quadrant)/4 length(quadrant)/4 length(quadrant)/4]);
if h
    text(0.5,-0.8,'*');
end

%% effect of cue plots

AllSel=zeros(length(IncludedNeurons),1);
index=find(IncludedNeurons);
Sel=sucvalues>0 & malvalues<0;
AllSel(index)=Sel;
Sel=logical(AllSel);

ymax=7;
alph=0.3;
trialsback=10;

Xaxis=[-0.5 2.5];
Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
time1=R.Param.Tm(Ishow);

%plot sucrose after non-specific cue
psth1=nanmean(R.Ev(RD1NP).PSTHz(Sel,Ishow),1); 
sem1=nanste(R.Ev(RD1NP).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up1=psth1+sem1;
down1=psth1-sem1;

%plot sucrose after sucrose cue
psth2=nanmean(R.Ev(RD1PR).PSTHz(Sel,Ishow),1); 
sem2=nanste(R.Ev(RD1PR).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up2=psth2+sem2;
down2=psth2-sem2;

%plot malt after non-specific cue
psth3=nanmean(R.Ev(RD2NP).PSTHz(Sel,Ishow),1); 
sem3=nanste(R.Ev(RD2NP).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up3=psth3+sem3;
down3=psth3-sem3;

%plot malt after maltodextrin cue
psth4=nanmean(R.Ev(RD2PR).PSTHz(Sel,Ishow),1); 
sem4=nanste(R.Ev(RD2PR).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up4=psth4+sem4;
down4=psth4-sem4;

%make the plot
subplot(3,6,6);
hold on;
title(['Suc>0, Mal<0 (n=' num2str(sum(Sel)) ')']);
patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');

a=plot(time1,psth1,'Color',Colors('sucrose'),'linewidth',1);
b=plot(time1,psth2,'Color',Colors('blockrose'),'linewidth',1);
c=plot(time1,psth3,'Color',Colors('maltodextrin'),'linewidth',1);
d=plot(time1,psth4,'Color',Colors('blockrodextrin'),'linewidth',1);
patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('blockrose'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up4,down4(end:-1:1)],Colors('blockrodextrin'),'EdgeColor','none');alpha(alph);
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
plot([0 0],[-3 ymax],'color','k','linewidth',0.5);
text(0.1,ymax-0.5,'RD','color','k');
axis([Xaxis(1) Xaxis(2) -2 ymax]);
ylabel('Mean firing (z-score)');
xlabel('Seconds from RD');
legend([a b c d],'Non-pr suc','Pred suc','Non-pr mal','Pred mal','location','northeast');

disp('RD cue impact');
%% bar graph for cue effects on cue

%when the model data is only on subselected neurons           
NNumber=[1:length(R.Ninfo)]';
index=NNumber(IncludedNeurons);
cuemasksc=zeros(length(R.Ninfo),2);
cuemasksc(index,1)=bm_cue.mask_base_cue'|bm_cue.mask_mean_cue';
cuemasksc(index,2)=bm_cue.mask_base'|bm_cue.mask_mean';

stacked_data = zeros(2,2);


for subset=1:2
    stacked_data(1,subset)=sum(IncludedNeurons&cuemasksc(:,subset));
end

subplot(3,6,10)
hold on;
b = bar(stacked_data./sum(stacked_data,2),'stacked','edgecolor','none');
b(2).FaceColor = 'Flat';
b(2).CData = [0 0 0];
b(1).FaceColor = 'Flat';
b(1).CData = Colors('cueeffect');
axis([0 7 0 1])
text(1.5,0.7,'No cue effect','color',[0 0 0]);
%text(1.5,0.15,'Info','color',[0.3 0 0.8]);
text(1.5,0.2,'Cue effect','color',Colors('cueeffect'));
set(gca,'xtick',[]);
yticks([0,.2,.4,.6,.8,1])
ylabel('Fraction of population')

%% scatterplot of cue responses

subplot(3,6,11);
sz=25;
hold on;
sucvalues=nan(sum(IncludedNeurons),1);
malvalues=nan(sum(IncludedNeurons),1);

%mean
Sel=bm_cue.mask_mean_cue;
index=find(Sel);
vsuc = getParameters_ott(os, 'cue', 'mean_cue', 'vsuc');
vmal = getParameters_ott(os, 'cue', 'mean_cue', 'vmal');
scatter(vsuc(Sel),vmal(Sel),sz,Colors('cueeffect'));
sucvalues(index)=vsuc(Sel)';
malvalues(index)=vmal(Sel)';

%base
Sel=bm_cue.mask_base_cue;
index=find(Sel);
vsuc = getParameters_ott(os, 'cue', 'base_cue', 'vsuc');
vmal = getParameters_ott(os, 'cue', 'base_cue', 'vmal');
scatter(vsuc(Sel),vmal(Sel),sz,Colors('cueeffect'));
sucvalues(index)=vsuc(Sel)';
malvalues(index)=vmal(Sel)';

axis([-1 1 -1 1]);
xlabel('Sucrose value');
ylabel('Maltodextrin value');
plot([0 0],[-1 1],'color','k');
plot([-1 1],[0 0],'color','k');

quadrant=[];
for i=1:length(sucvalues)
    if sucvalues(i)<0 & malvalues(i)>0
        quadrant(i)=1;
    elseif sucvalues(i)>0 & malvalues(i)>0
        quadrant(i)=2;
    elseif sucvalues(i)<0 & malvalues(i)<0
        quadrant(i)=3;
    elseif sucvalues(i)>0 & malvalues(i)<0
        quadrant(i)=4;
    end
end
quadrant=quadrant(find(quadrant));
text(-0.9,0.8,[num2str(round(100*sum(quadrant==1)/length(quadrant),1)) '%']);
text(0.6,0.8,[num2str(round(100*sum(quadrant==2)/length(quadrant),1)) '%']);
text(-0.9,-0.8,[num2str(round(100*sum(quadrant==3)/length(quadrant),1)) '%']);
text(0.6,-0.8,[num2str(round(100*sum(quadrant==4)/length(quadrant),1)) '%']);

[h,Stats.Cuepval,Stats.CueX2] = chi2gof(quadrant,'Ctrs',[1 2 3 4],'expected',[length(quadrant)/4 length(quadrant)/4 length(quadrant)/4 length(quadrant)/4]);
if h
    text(0.5,-0.8,'*');
end

%% effect of cue plots

AllSel=zeros(length(IncludedNeurons),1);
index=find(IncludedNeurons);
Sel=sucvalues>0 & malvalues<0;
AllSel(index)=Sel;
Sel=logical(AllSel);

ymax=8;
alph=0.3;
trialsback=10;

Xaxis=[-0.5 1.5];
Ishow=find(R.Param.Tm>=Xaxis(1) & R.Param.Tm<=Xaxis(2));
time1=R.Param.Tm(Ishow);

%plot non-specific cue, sucrose trials
psth1=nanmean(R.Ev(Cue31).PSTHz(Sel,Ishow),1); 
sem1=nanste(R.Ev(Cue31).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up1=psth1+sem1;
down1=psth1-sem1;

%plot sucrose cue
psth2=nanmean(R.Ev(Cue1).PSTHz(Sel,Ishow),1); 
sem2=nanste(R.Ev(Cue1).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up2=psth2+sem2;
down2=psth2-sem2;

%plot non-specific cue, maltodextrin trials
psth3=nanmean(R.Ev(Cue32).PSTHz(Sel,Ishow),1); 
sem3=nanste(R.Ev(Cue32).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up3=psth3+sem3;
down3=psth3-sem3;

%plot maltodextrin cue
psth4=nanmean(R.Ev(Cue2).PSTHz(Sel,Ishow),1); 
sem4=nanste(R.Ev(Cue2).PSTHz(Sel,Ishow),1); %calculate standard error of the mean
up4=psth4+sem4;
down4=psth4-sem4;

%make the plot
subplot(3,6,12);
hold on;
title(['Suc>0, Mal<0 (n=' num2str(sum(Sel)) ')']);
patch([CueWindow(1) CueWindow(1) CueWindow(2) CueWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');

a=plot(time1,psth1,'Color',Colors('sucrose'),'linewidth',1);
b=plot(time1,psth2,'Color',Colors('blockrose'),'linewidth',1);
c=plot(time1,psth3,'Color',Colors('maltodextrin'),'linewidth',1);
d=plot(time1,psth4,'Color',Colors('blockrodextrin'),'linewidth',1);
patch([time1,time1(end:-1:1)],[up1,down1(end:-1:1)],Colors('sucrose'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up2,down2(end:-1:1)],Colors('blockrose'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up3,down3(end:-1:1)],Colors('maltodextrin'),'EdgeColor','none');alpha(alph);
patch([time1,time1(end:-1:1)],[up4,down4(end:-1:1)],Colors('blockrodextrin'),'EdgeColor','none');alpha(alph);
plot([-2 5],[0 0],':','color','k','linewidth',0.75);
plot([0 0],[-3 ymax],'color','r','linewidth',0.5);
text(Xaxis(1)+0.1,ymax-0.5,'Cue','color','r');
axis([Xaxis(1) Xaxis(2) -2 ymax]);
ylabel('Mean firing (z-score)');
xlabel('Seconds from cue');
legend([a b c d],'Non-pr suc','Pred suc','Non-pr mal','Pred mal','location','northeast');
disp('Cue cue impact')



%% RPE neurons based on different RPE values

bins=[-3 -1.2;-1.2 -0.7;-0.7 -0.4;-0.4 0;0 0.4;0.4 0.7;0.7 1.2;1.2 3];

Sel=bm_RD.mask_base';
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
    
    neuroninfo = R.Ninfo(IncludedNeurons,:);
    session=neuroninfo(index(i),1); %what's the name of session
    sessnum=ismember(sessions,session); %which session of all sessions
    neuronnum=ismember(R.Ninfo(strcmp(session,R.Ninfo(:,1)),2),neuroninfo(index(i),2)); %which number neuron in the session
    RD=strmatch('RD',RAW(sessnum).Einfo(:,2),'exact'); %get index for reward delivery times
    Cue=strmatch('Cues',RAW(sessnum).Einfo(:,2),'exact'); %get index for reward delivery times
    
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
ymax=8;
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


subplot(3,3,7);
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
axis([Xaxis(1) Xaxis(2) -1 ymax]);
ylabel('Mean firing (z-score)');
xlabel('Seconds from RD');
colormap(psthcolor);
h=colorbar;
set(h,'YTick',[0:0.125:1]);
set(h,'yticklabel',[-3;-1.2;-.7;-.4;0;0.4;0.7;1.2;3])
set(get(h,'title'),'string','RPE');
%legend([a{:}],labels,'location','northeast');

disp('RPE responses')
