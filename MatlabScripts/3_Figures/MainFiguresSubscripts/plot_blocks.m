function R=plot_blocks(Colors, R, events, trial_info, plot_sub_figures,bm_RD)
%PLOT_PAPER_F3 Plots paper figure 3 for the blocks project. This is the
%evolution through the trial figure

xaxis=[-2 4];

num_neurons = size(events{1,2},1);

interspersed = [trial_info.blocks]' == 0;
blocks1 = [trial_info.blocks]' == 1;
blocks2 = [trial_info.blocks]' == 2;
area_vectors = [strcmp({trial_info.area},'NA')',strcmp({trial_info.area},'VP')'];

suc_delivery=strcmp('RD1', R.Erefnames);
malt_delivery=strcmp('RD2', R.Erefnames);
suc_cue=strcmp('CueR1', R.Erefnames);
malt_cue=strcmp('CueR2', R.Erefnames);
time_indicies=find(R.Param.Tm>=xaxis(1) & R.Param.Tm<=xaxis(end)); %returns the indicies where these are both true
activity_xaxis=R.Param.Tm(time_indicies);

masks=cat(2,bm_RD.mask_base',bm_RD.mask_curr',bm_RD.mask_mean');
region = strcmp(R.Region,'VP');

HistWindow=[0.75 1.95];

color = {Colors('sucrose'),Colors('maltodextrin'),Colors('total'),Colors('blockrose'),Colors('blockrodextrin'),Colors('total')};
figure;

%% bar graphs
stacked_data = zeros(2,3);

for session_type=1:2
    for subset=1:3
        stacked_data(session_type,subset)=sum(region&masks(:,subset)&R.Blocks==(session_type-1));
    end
end

subplot(5,4,1)
hold on;
b = bar(stacked_data./sum(stacked_data,2),'stacked','edgecolor','none');
b(3).FaceColor = 'Flat';
b(3).CData = Colors('mean');
b(2).FaceColor = 'Flat';
b(2).CData = Colors('current');
b(1).FaceColor = 'Flat';
b(1).CData = Colors('rpe');
axis([0 5 0 1])

xticks([1 2]);
xticklabels({'Interspersed','Blocked'});
yticks([0,.2,.4,.6,.8,1])
ylabel('Fraction of population')
view(-90,90);
legend('RPE','Current','Mean','location','northwest');
%chisquared
%all categories
category=sum(masks.*[1:3],2);
category=category(region);
task=R.Blocks(region);
[~,x2,pval]=crosstab(category,task);
R.Stats.MLEX2(1).chi2=x2;
R.Stats.MLEX2(1).pval=pval;

%RPE cells
category=masks(region,1);
task=R.Blocks(region);
[~,x2,pval]=crosstab(category,task);
R.Stats.MLEX2(2).chi2=x2;
R.Stats.MLEX2(2).pval=pval;
%% activity plots
ymax=11;

for session_type=1:2
    if session_type==1 Sel=region & R.Blocks==0 & bm_RD.mask_base'; end
    if session_type==2 Sel=region & R.Blocks & bm_RD.mask_base'; end
    
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

    subplot(5,8,[5:8]-8*(session_type-2))
    hold on
    patch([HistWindow(1) HistWindow(1) HistWindow(2) HistWindow(2)],[-5 12 12 -5],Colors('gray'),'EdgeColor','none');
    a=plot(activity_xaxis,psth1,'Color',color{1+3*(session_type-1)},'linewidth',1);
    b=plot(activity_xaxis,psth2,'Color',color{2+3*(session_type-1)},'linewidth',1);
    patch([activity_xaxis,activity_xaxis(end:-1:1)],[up1,down1(end:-1:1)],color{1+3*(session_type-1)},'EdgeColor','none');alpha(0.5);
    patch([activity_xaxis,activity_xaxis(end:-1:1)],[up2,down2(end:-1:1)],color{2+3*(session_type-1)},'EdgeColor','none');alpha(0.5);
    plot([-2 5],[0 0],':','color','k','linewidth',0.75);
    plot([0 0],[-2 ymax],'color','k','linewidth',0.5);
    plot([-0.5 -0.5],[-2 ymax],'color','b','linewidth',0.5);
    axis([-2 4 -1 ymax]);
    text(-1.75,ymax+0.5,['RPE (n=' num2str(sum(Sel)) ' of ' num2str(sum(region&R.Blocks==(session_type-1))) ')'],'FontSize',8)
    legend([a b],'Suc trials','Mal trials');
    text(-1.3,ymax-0.5,'PE','color','b');
    text(0.1,ymax-0.5,'RD','color','k');
    xlabel('Seconds from RD');
    set(gca,'ytick',[])

    %plot sucrose preferring neurons to sucrose cue
    psth1=nanmean(R.Ev(suc_cue).PSTHz(Sel,time_indicies),1); 
    sem1=nanste(R.Ev(suc_cue).PSTHz(Sel,time_indicies),1); %calculate standard error of the mean
    up1=psth1+sem1;
    down1=psth1-sem1;

    %plot sucrose preferring neurons to malt cue
    psth2=nanmean(R.Ev(malt_cue).PSTHz(Sel,time_indicies),1); 
    sem2=nanste(R.Ev(malt_cue).PSTHz(Sel,time_indicies),1); %calculate standard error of the mean
    up2=psth2+sem2;
    down2=psth2-sem2;

    subplot(5,8,[3:4]-8*(session_type-2))
    hold on
    plot(activity_xaxis,psth1,'Color',color{1+3*(session_type-1)},'linewidth',1);
    plot(activity_xaxis,psth2,'Color',color{2+3*(session_type-1)},'linewidth',1);
    patch([activity_xaxis,activity_xaxis(end:-1:1)],[up1,down1(end:-1:1)],color{1+3*(session_type-1)},'EdgeColor','none');alpha(0.5);
    patch([activity_xaxis,activity_xaxis(end:-1:1)],[up2,down2(end:-1:1)],color{2+3*(session_type-1)},'EdgeColor','none');alpha(0.5);
    plot([-2 5],[0 0],':','color','k','linewidth',0.75);
    axis([-1 2 -1 ymax]);
    ylabel('Mean firing (z-score)');
    xlabel('Seconds from cue onset')
    plot([0 0],[-1 ymax],'color','r','linewidth',0.5);
    text(-0.6,ymax-0.5,'Cue','color','r');
    title('Event-evoked firing');

end
disp('Event firing')

%% firing across session
session_types = {interspersed,blocks1,blocks2};
event=3; %RD
area=2; %VP
SessType={'RPE','Current only','Mean'};

% Visualize (equal axis)
for subset=1:3
    for i = 1:3 % Each session type (interspersed, blocks1, blocks2)
        subplot(5,3,6+i+3*(subset-1)); 
        hold on


        alpha(0.15)

        data_suc = events{event,5}(session_types{i}&area_vectors(:,area)&masks(:,subset),:);
        error_suc = std(data_suc)/sqrt(size(data_suc,1));

        data_means = mean(data_suc);            
        error_bars = error_suc;
        if i == 1
            segments = [1 3 5 7 9;2 4 6 8 10]; % Change to [1:6] to make all connected
        elseif i == 2
            segments = [1:5;6:10];
        else
            segments = [6:10;1:5];
        end

        segcolors={Colors('sucrose'),Colors('maltodextrin'),Colors('blockrose'),Colors('blockrodextrin'),Colors('blockrose'),Colors('blockrodextrin')};
        segnum=0;                      

        for segment = segments'
            segnum=segnum+1;
            if i==1
                l{segnum}=errorbar(1:5,data_means(1,segment),error_bars(1,segment),'Color',segcolors{2*(i-1)+segnum},'linewidth',1.5);
                axis([0 6 -1 11])

            else
                l{segnum}=errorbar(segment,data_means(1,segment),error_bars(1,segment),'Color',segcolors{2*(i-1)+segnum},'linewidth',1.5);
                axis([0 11 -1 11])
            end
        end

%             title([session_titles{i} ': ' events{event,1}])
        if i < 3
            ylabel('Mean firing (z-score)')
            legend([l{1} l{2}],'Suc trials','Mal trials','AutoUpdate','off');
        end

        plot([0 11],[0 0],':','color','k');
        set(gca,'xtick',[]);
        xlabel('Session progress');
        title({[SessType{subset} ' (n=' num2str(length(data_suc)) ' of ' num2str(sum(session_types{i}&area_vectors(:,area))) ')']}');
        
    end
    
end
disp('Evolution')

if plot_sub_figures
    
    % Make the Linear Model
    linear_models = cell(num_neurons,size(events,1));
    linear_models_shuffled = cell(num_neurons,size(events,1));
    model_output = cell(3,1);
    model_output_shuffled = cell(3,1);

    for task = 1:3
        for class=1:3
            Sel=session_types{task} & masks(:,class) & region;
            SelIndex=find(Sel);
            data_table_all{class,task}=table;
            for event = 3 
                for neuron = 1:length(SelIndex)
                    firing_rate = events{event, 2}(SelIndex(neuron),:)';
                    firing_rate = firing_rate(~isnan(firing_rate)); 

                    if event == 3
                        reward = trial_info(SelIndex(neuron)).current_reward; 
                        progress = (1:length(firing_rate))'./length(firing_rate);
                        order = ones(length(reward),1)*trial_info(SelIndex(neuron)).blocks;
                        neuronnum = categorical(ones(length(reward),1)*SelIndex(neuron));

                    end

                    data_table = table(progress,reward,order,neuronnum,firing_rate);
                    data_table = data_table(~any(ismissing(data_table),2),:); % Removes any NaN values from Cues played at the beginning without any previous reward (should be very few, since I already removed the very first trial)

                    data_table_all{class,task}=cat(1,data_table_all{class,task},data_table);
                end  
            end
            
            R.sucrose_mdl{class,task} = fitlm(data_table_all{class,task}.progress(data_table_all{class,task}.reward==1),data_table_all{class,task}.firing_rate(data_table_all{class,task}.reward==1));
            R.maltodextrin_mdl{class,task} = fitlm(data_table_all{class,task}.progress(data_table_all{class,task}.reward==0),data_table_all{class,task}.firing_rate(data_table_all{class,task}.reward==0));
            R.lmdl{class,task}=fitlm(data_table_all{class,task},'firing_rate~progress+reward+progress:reward');
        end
    end

    R.lmdl_blocksRPE=fitlme(cat(1,data_table_all{1,2},data_table_all{1,3}),'firing_rate~progress+reward+progress:reward+progress:reward:order+(1|neuronnum)');
    R.lmdl_blocksCurrent=fitlme(cat(1,data_table_all{2,2},data_table_all{2,3}),'firing_rate~progress+reward+progress:reward+progress:reward:order+(1|neuronnum)');

    disp('Linear Models')

end

