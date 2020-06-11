clear all;

files={'42Sess2foreverDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '42Sess3DeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '42Sess4DeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '42Sess4redoDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '42Sess5DeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '51Sess1againDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '51Sess1foreverDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '51Sess2redoDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '53Sess1againDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '53Sess1redoDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '54Sess1foreverDeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '54Sess1DeepCut_resnet50_ThreeAug21shuffle1_1030000.h5';...
        '54Sess2DeepCut_resnet50_ThreeAug21shuffle1_1030000.h5'};
    
portcoord=[50 65;50 65;50 65;50 65;50 65;60 70;60 70;60 70;60 70;60 70;60 70;60 70;60 70];

load('RAWintBlocks');
load('ModData_intBlocks.mat');
load ('intBlocks_MLEfits.mat');
bm_RD=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
bm_cue=select_RPEmods(os, 'cue', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
    
for i=1:length(RAWblocks)
    name=char(RAWblocks(i).Region);
    region{i,1}=name(1:3);
    name=char(RAWblocks(i).Blocks);
    type{i,1}=name(1);
end
RAW=RAWblocks((strcmp('VP2',region) | strcmp('VP3',region) | strcmp('VP4',region) | strcmp('VP5',region)) & strcmp('I',type));

included=(strcmp('VP2',CS.Rat) | strcmp('VP3',CS.Rat) | strcmp('VP4',CS.Rat) | strcmp('VP5',CS.Rat)) & CS.Blocks==0;
%included=(strcmp('VP2',CS.Rat) | strcmp('VP3',CS.Rat) | strcmp('VP4',CS.Rat)) & CS.Blocks==0;
Predictors=CS.Predictors(included);
RDHz=CS.RDHz(included);
CueHz=CS.CueHzAll(included);

Colors = load_colors();
[magma,inferno,plasma,viridis]=colormaps;
masks=cat(2,bm_RD.mask_base(included)',bm_RD.mask_curr(included)',bm_RD.mask_mean(included)');
cue_masks=cat(2,bm_cue.mask_base(included)',bm_cue.mask_mean(included)');


%%
NN=0;
trlactivity={[];[];[]};
trldistance={[];[];[]};
all_distance=[];
all_preds=[];

for session=1:length(RAW)
    %get coordinates from deeplabcut analysis
    [timestamps,xcoordinates,ycoordinates]=AnalyzeDeepLabCut(files{session});
    
    %get event times
    cue = strcmp('Cue',RAW(session).Einfo(:,2));
    cuetimes = RAW(session).Erast{cue};
    pe = strcmp('PE',RAW(session).Einfo(:,2));
    petimes = RAW(session).Erast{pe};
    rd = strcmp('RD',RAW(session).Einfo(:,2));
    rdtimes = RAW(session).Erast{rd};
    lick = strcmp('Licks',RAW(session).Einfo(:,2));
    licks = RAW(session).Erast{lick};
    
    %find included trials
    included_RDs=rdtimes<cuetimes(end);
    included_cue_trials=[];
    for trial=1:length(rdtimes)
        if sum(cuetimes>rdtimes(trial))>0
            included_cue_trials(trial,1)=find(cuetimes>rdtimes(trial),1,'first');
        end
    end
    
    %for each cue onset
    xypoints={};
    xtrial=[];
    ytrial=[];
    ITI_dfp=[];
    
    for trial=1:length(cuetimes)
        
        %coordinates at cue onset
        if sum(timestamps>(cuetimes(trial)-0.15) & timestamps<(cuetimes(trial)+0.15))>0
            xtrial(trial)=mean(xcoordinates((timestamps>cuetimes(trial)-0.15) & (timestamps<cuetimes(trial)+0.15)));
            ytrial(trial)=mean(ycoordinates((timestamps>cuetimes(trial)-0.15) & (timestamps<cuetimes(trial)+0.15)));      
        else
            xtrial(trial)=NaN;
            ytrial(trial)=NaN;
        end
        
        %total distance traveled during ITI, and mean distance from port
        binsize=0.2; %seconds

        if sum(rdtimes<cuetimes(trial))>0 & sum(included_cue_trials==trial)>0
            rd_time=max(rdtimes(rdtimes<cuetimes(trial)));
            start_time=max(licks(licks>rd_time & licks<(rd_time+15)));
            %bins=start_time:binsize:cuetimes(trial);
            bins=start_time:binsize:cuetimes(trial)+binsize/2;
            bn=0;
            bin_t=[];
            bin_x=[];
            bin_y=[];
            bin_dfp=[];
            
            %location for each bin
            for bin=1:length(bins)-1
                if sum(timestamps>bins(bin) & timestamps<bins(bin+1))>0
                    bn=bn+1;
                    bin_t(bn,1)=bins(bin)+binsize/2;
                    bin_x(bn,1)=mean(xcoordinates(timestamps>bins(bin) & timestamps<bins(bin+1)));
                    bin_y(bn,1)=mean(ycoordinates(timestamps>bins(bin) & timestamps<bins(bin+1)));
                    bin_dfp(bn,1)=sqrt((bin_x(bn,1)-portcoord(session,1))^2+(bin_y(bn,1)-portcoord(session,2))^2);
                end
            end
            
            norm_time=(bin_t-bin_t(1))/(bin_t(end)-bin_t(1));
            
            xypoints{trial,1}=cat(2,bin_x,bin_y);
            
            ITI_dfp(trial)=trapz(bin_t,bin_dfp)/(bin_t(end)-bin_t(1)); %dfp = distance from port
            
           
        else
            ITI_dfp(trial)=NaN;
            xypoints{trial,1}=NaN;
        end
    end
    
    %only look at included trials
    distance_inc=ITI_dfp(included_cue_trials)';    
   
    %relate to neural activity
    noi=[1:length(RAW(session).Nrast)]+NN;
    noi_activity=RDHz(noi);
    noi_activity_cue=CueHz(noi);

    sessions=[5];
    for neuron=1:length(noi)
        NN=NN+1;
        
        %plot example traces
        if neuron==1 & sum(sessions==session)>0
             
            %scatterplot of rat locations
            xypoints=xypoints(included_cue_trials,:);
            xtrial_inc=xtrial(included_cue_trials);
            ytrial_inc=ytrial(included_cue_trials);
            sucrosexy=cat(1,xypoints{Predictors{NN,1}(1:sum(included_RDs),1)==1});
            maltodextrinxy=cat(1,xypoints{Predictors{NN,1}(1:sum(included_RDs),1)==0});
            opacity=0.1;
            dotsize=24;
            
            sucrosetrls=xypoints(Predictors{NN,1}(1:sum(included_RDs),1)==1);
            maltodextrintrls=xypoints(Predictors{NN,1}(1:sum(included_RDs),1)==0);
            
            %sucrose traces
            figure;
            subplot(2,3,1);
            hold on;
            s1=scatter(sucrosexy(:,1),sucrosexy(:,2),dotsize,Colors('sucrose'),'filled');
            s1.MarkerFaceAlpha = opacity;
            sc=scatter(xtrial_inc(Predictors{NN,1}(1:sum(included_RDs),1)==1),ytrial_inc(Predictors{NN,1}(1:sum(included_RDs),1)==1),'k','x');

            
            axis([0 450 0 350]);
            set(gca,'Ydir','reverse')
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            if session>5 plot([40 60],[25 70],'color','k'); end
            if session<=5 plot([40 50],[25 65],'color','k'); end
            legend(sc,'Location at cue onset');
            text(5,10,'Reward Port');
            title('Location post-sucrose');

            
            %maltodextrin traces
            subplot(2,3,2);
            hold on;
            s2=scatter(maltodextrinxy(:,1),maltodextrinxy(:,2),dotsize,Colors('maltodextrin'),'filled');
            s2.MarkerFaceAlpha = opacity;
            scatter(xtrial_inc(Predictors{NN,1}(1:sum(included_RDs),1)==0),ytrial_inc(Predictors{NN,1}(1:sum(included_RDs),1)==0),[],'k','x')
            
            %plot individual trial traces
%             for i=1:length(maltodextrintrls)
%                 if maltodextrintrls{i,1}
%                     plot(maltodextrintrls{i,1}(:,1),maltodextrintrls{i,1}(:,2),'color',Colors('maltodextrin'));
%                 end
%             end
            
            set(gca,'Ydir','reverse')
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            axis([0 450 0 350]);
            if session>5 plot([40 60],[25 70],'color','k'); end
            if session<=5 plot([40 50],[25 65],'color','k'); end
            text(5,10,'Reward Port');
            title('Location post-maltodextrin');

          end


        [distance_corr_cue(NN,1),distance_corr_cue(NN,2)]=corr(noi_activity_cue{neuron,1}(included_cue_trials,1),distance_inc,'rows','complete','type','spearman');
        [distance_corr(NN,1),distance_corr(NN,2)]=corr(noi_activity{neuron,1}(included_RDs,1),distance_inc,'rows','complete','type','spearman');

        
        %shuffled controls
        for i=1:1000
            distance_shuff=distance_inc(randperm(length(distance_inc)));
            dist_shuff_corr(NN,i)=corr(noi_activity{neuron,1}(included_RDs,1),distance_shuff,'rows','complete','type','spearman');
            dist_shuff_corr_cue(NN,i)=corr(noi_activity_cue{neuron,1}(included_cue_trials,1),distance_shuff,'rows','complete','type','spearman');

        end
        
        if neuron==1 %only do this once per session
        
            suc_distance(session,1)=nanmean(distance_inc(Predictors{NN,1}(1:sum(included_RDs),1)==1));
            mal_distance(session,1)=nanmean(distance_inc(Predictors{NN,1}(1:sum(included_RDs),1)==0));
            
            distance_inc_norm = (distance_inc - nanmean(distance_inc))/nanstd(distance_inc);
            
            all_distance = cat(1,all_distance,distance_inc_norm);
            all_preds = cat(1,all_preds,Predictors{NN,1}(1:sum(included_RDs),:));
            
        end
        
 

        
        
    end
    
    disp(['Session #' num2str(session)]);
    
end

%% plotting distance from port graphs

subplot(2,6,6);
hold on;
plot([1 2],[suc_distance mal_distance],'color',[0.6 0.6 0.6]);
errorbar(1,nanmean(suc_distance),nanste(suc_distance,1),'color',Colors('sucrose'),'marker','o','linewidth',1.5);
errorbar(2,nanmean(mal_distance),nanste(mal_distance,1),'color',Colors('maltodextrin'),'marker','o','linewidth',1.5);
ylabel('Distance from port during ITI (pixels)');
xticks([1 2]);
xticklabels({'Post-suc','Post-mal'});
xtickangle(45);
axis([0.5 2.5 0 200]);

%stats
%signrank(suc_distance,mal_distance);

colors{1,1}=Colors('rpe');
colors{2,1}=Colors('current');
colors{3,1}=Colors('mean');

subplot(2,3,5);
hold on;
plots={};
for mask=1:length(masks(1,:))
    [cdf,x] = ecdf(distance_corr(masks(:,mask)));
    plots{mask} = plot(x,cdf,'linewidth',1.5,'color',colors{mask,1});
    plot([mean(distance_corr(masks(:,mask))) mean(distance_corr(masks(:,mask)))],[0 1],'color',colors{mask,1},'linewidth',1);
end
%histogram(distance_corr(masks(:,1)),-0.5:0.05:0.5,'normalization','probability','edgecolor','none','facecolor',Colors('rpe'));
plot([0 0],[0 1],'color','k');
plot([-0.5 0.5],[0.5 0.5],'color','k');
axis([-.5 .5 0 1]);
legend([plots{:}],'RPE','Current','Unmod.','location','southeast')
xlabel('Spearman''s rho');
ylabel('Cumulative fraction of neurons');
title('Correlation between VP activity and ITI distance');
text(-0.2,1,'*','color','k','fontsize',24);

%stats tests
%ranksum(distance_corr(masks(:,1)),distance_corr(masks(:,3)));
%signrank(distance_corr(masks(:,2)),dist_shuff_corr(masks(:,2)));

%% cue value
colors{1,1}=Colors('rpe');
colors{2,1}=Colors('mean');

figure;

subplot(1,1,1);
hold on;
plots={};
for mask=1:length(cue_masks(1,:))
    [cdf,x] = ecdf(distance_corr_cue(cue_masks(:,mask)));
    plots{mask} = plot(x,cdf,'linewidth',1.5,'color',colors{mask,1});
    plot([mean(distance_corr_cue(cue_masks(:,mask))) mean(distance_corr_cue(cue_masks(:,mask)))],[0 1],'color',colors{mask,1},'linewidth',1);
end
%histogram(distance_corr(masks(:,1)),-0.5:0.05:0.5,'normalization','probability','edgecolor','none','facecolor',Colors('rpe'));
plot([0 0],[0 1],'color','k');
plot([-0.5 0.5],[0.5 0.5],'color','k');
axis([-.5 .5 0 1]);
legend([plots{:}],'Value','Unmod.','location','northwest')
xlabel('Spearman''s rho');
ylabel('Cumulative fraction of neurons');
title('Correlation between VP activity and ITI distance');
text(-0.2,1,'*','color','k','fontsize',24);

%stats tests
%ranksum(distance_corr_cue(cue_masks(:,1)),distance_corr(cue_masks(:,2)));
%signrank(distance_corr_cue(masks(:,1)),dist_shuff_corr_cue(masks(:,1)));
