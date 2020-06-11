clear all;
analyzedlcfiles=0;
%% get data from DeepLabCut files

if analyzedlcfiles

    ffmpeg_path='C:\FFmpeg\bin';
    address=['']; %where are the running files located
    files=dir([address,'\\*_*']);

    videoaddress=[''];
    videofiles=dir([videoaddress,'\\*Ch*']);

    for k=1:length(files)
        filename=fullfile(files(k).folder,files(k).name);
        rat=str2num(files(k).name(9:10)); %rat #
        videofile=[videofiles(k).folder '\' videofiles(k).name];

        %get timestamp data
        [timestamps,xcoordinates,ycoordinates]=AnalyzeDeepLabCutPP(filename,videofile,ffmpeg_path);

        VideoData{rat,1}=timestamps;
        VideoData{rat,2}=xcoordinates;
        VideoData{rat,3}=ycoordinates;

    end

    save('CueVideoData.mat','VideoData');

else
    load('CueVideoData.mat');

end
%%

spreadsheet=importdata('CueVideoKey.xlsx');
spreadsheetrat=spreadsheet.textdata(2:end,1);
for i=1:length(spreadsheetrat)
    sprat(i,1)=str2num(spreadsheetrat{i,1}(3:4));
end

lighton(sprat,1)=spreadsheet.data(:,1);
portcoord(sprat,1)=spreadsheet.data(:,2);
portcoord(sprat,2)=spreadsheet.data(:,3);
%cutoff(sprat,1)=spreadsheet.data(:,5);

Colors = load_colors();
trlactivity={[];[];[]};
trldistance={[];[];[]};

address=['/Data/ChR2CueBehavior'];
%files=dir([address,'\\*!*']);
files=dir([address,'//*!*']); %mac

for k=1:length(files)
    filename=fullfile(files(k).folder,files(k).name);
    file=fopen(filename);
    L=textscan(file,'%s','Delimiter',':');
    fclose('all');
    rat=str2num(files(k).name(end-1:end));

    if sum(rat==3 | rat==15)==0 %don't analyze rats without videos
        %find start of A, start of C, and end of C
        Jstrt=find(strcmp('J',L{1,1})); %Cue
        Kstrt=find(strcmp('K',L{1,1})); %Sucrose Delivery
        Ostrt=find(strcmp('O',L{1,1})); %Port entry times
        Pstrt=find(strcmp('P',L{1,1})); %Port durations
        Qstrt=find(strcmp('Q',L{1,1})); %Laser stimulations
        Xstrt=find(strcmp('X',L{1,1})); %List (marks end of laser stims)

        %trial times
        choicet=(L{1,1}(Jstrt+2:2:Kstrt-1));
        choicetime=[];
        for i=1:length(choicet)
            medtext=textscan(char(choicet(i,1)),'%f','Delimiter',' ');
            ctime=medtext{1,1}(:,1);
            choicetime=(cat(1,choicetime,ctime(isnan(ctime)==0)));
        end
        choicetime(choicetime==0)=[]; %delete 0s


        %Sucrose delivery times
        sucrosedel=(L{1,1}(Kstrt+2:2:Ostrt-1));
        sucrosedeltime=[];
        for i=1:length(sucrosedel)
            medtext=textscan(char(sucrosedel(i,1)),'%f','Delimiter',' ');
            sdtime=medtext{1,1}(:,1);
            sucrosedeltime=(cat(1,sucrosedeltime,sdtime(isnan(sdtime)==0)));
        end
        sucrosedeltime(sucrosedeltime==0)=[]; %delete 0s

        %Port entry times
        pet=(L{1,1}(Ostrt+2:2:Pstrt-1));
        petime=[];
        for i=1:length(pet)
            medtext=textscan(char(pet(i,1)),'%f','Delimiter',' ');
            portentry=medtext{1,1}(:,1);
            petime=(cat(1,petime,portentry(isnan(portentry)==0)));
        end
        removes=petime==0;
        petime(petime==0)=[]; %delete 0s

        %Port entry durations
        durationt=(L{1,1}(Pstrt+2:2:Qstrt-1));
        durationtime=[];
        for i=1:length(durationt)
            medtext=textscan(char(durationt(i,1)),'%f','Delimiter',' ');
            duration=medtext{1,1}(:,1);
            durationtime=(cat(1,durationtime,duration(isnan(duration)==0)));
        end
        durationtime(removes)=[]; %delete 0s

        pexits=petime(1:length(durationtime))+durationtime;

        %Laser stims
        ST=(L{1,1}(Qstrt+2:2:Xstrt-1));
        STtime=[];
        for i=1:length(ST)
            medtext=textscan(char(ST(i,1)),'%f','Delimiter',' ');
            stime=medtext{1,1}(:,1);
            STtime=(cat(1,STtime,stime(isnan(stime)==0)));
        end

        STtime(STtime==0)=[]; %delete 0s    


        %Sucrose or water?
        latency=[];
        for i=1:length(choicetime)
            sucrose=find(sucrosedeltime(:,1) > choicetime(i,1),1);
            sucrosetime=sucrosedeltime(sucrose);
            if sucrosetime
                latency(i,1)=sucrosetime-choicetime(i,1)-0.5; %delay is 0.5s
            else
                latency(i,1)=NaN;
            end
        end
        latency(latency>10)=NaN;

        completed=latency>0;

        %stim trials
        post_stimmed=zeros(length(choicetime),1);
        for stim=1:length(STtime)
            trialtime=min(choicetime(choicetime>STtime(stim)));
            if ~isempty(trialtime)
                post_stimmed(choicetime==trialtime)=1;
            end
        end

        %no stim trials
        postsuc=zeros(length(choicetime),1);
        for suctrial=1:length(sucrosedeltime)
            trialtime=min(choicetime(choicetime>sucrosedeltime(suctrial)));
            if ~isempty(trialtime)
                postsuc(choicetime==trialtime)=1;
            end
        end
        nostim=postsuc==1 & post_stimmed==0;    


        rdtimes=sucrosedeltime;
        cuetimes=choicetime;
        petimes=petime;

        included_cue_trials=[];
        for trial=1:length(rdtimes)
            if sum(cuetimes>rdtimes(trial))>0
                included_cue_trials(trial,1)=find(cuetimes>rdtimes(trial),1,'first');
            end
        end
        post_stim=post_stimmed(included_cue_trials);

        %get coordinates from deeplabcut analysis
        timestamps=VideoData{rat,1}-lighton(rat,1); %adjust for when session started in video
        xcoordinates=VideoData{rat,2};
        ycoordinates=VideoData{rat,3};


        %for each cue onset
        xypoints={};
        xtrial=[];
        ytrial=[];
        ITI_distance=[];
        ITI_dfp=[];
        licks_trl=[];


        for trial=1:length(cuetimes)
            if sum(timestamps>(cuetimes(trial)-0.25) & timestamps<(cuetimes(trial)+0.25))>0
                xtrial(trial)=mean(xcoordinates((timestamps>cuetimes(trial)-0.25) & (timestamps<cuetimes(trial)+0.25)));
                ytrial(trial)=mean(ycoordinates((timestamps>cuetimes(trial)-0.25) & (timestamps<cuetimes(trial)+0.25)));            



            elseif sum(timestamps>(cuetimes(trial)-1) & timestamps<(cuetimes(trial)+0.3))>0
                xtrial(trial)=mean(xcoordinates((timestamps>cuetimes(trial)-1) & (timestamps<cuetimes(trial)+0.3)));
                ytrial(trial)=mean(ycoordinates((timestamps>cuetimes(trial)-1) & (timestamps<cuetimes(trial)+0.3)));


            else
                xtrial(trial)=NaN;
                ytrial(trial)=NaN;
            end

            %total distance traveled during ITI, and mean distance from port
            binsize=0.2; %seconds

            if sum(rdtimes<cuetimes(trial))>0 & postsuc(trial)==1

                rd_time=max(rdtimes(rdtimes<cuetimes(trial)));
                pe_time=max(petimes(petimes<(rd_time+10)));
                start_time=min([pexits(petimes==pe_time) rd_time+10]);
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
                        bin_dfp(bn,1)=sqrt((bin_x(bn,1)-portcoord(rat,1))^2+(bin_y(bn,1)-portcoord(rat,2))^2);
                    end
                end
                
                xypoints{trial,1}=cat(2,bin_x,bin_y);
%                 subplot(5,5,rat);
%                 hold on;
%                 plot(bin_dfp);
%                 xypoints{trial,1}=cat(2,bin_x,bin_y);
%                 axis([0 25 0 500]);

                ITI_dfp(trial)=trapz(bin_t,bin_dfp)/(bin_t(end)-bin_t(1)); %dfp = distance from port


            else
                ITI_dfp(trial)=NaN;
                xypoints{trial,1}=NaN;
            end
        end

        

        %plot example traces
        if rat==17

            figure;

            %scatterplot of rat locations
            xypoints=xypoints(included_cue_trials,:);
            xtrial_inc=xtrial(included_cue_trials);
            ytrial_inc=ytrial(included_cue_trials);
            sucrosexy=cat(1,xypoints{post_stim==0});
            maltodextrinxy=cat(1,xypoints{post_stim==1});
            opacity=0.1;
            dotsize=24;

            sucrosetrls=xypoints(post_stim==0);
            maltodextrintrls=xypoints(post_stim==1);

            %sucrose traces
            subplot(1,4,1);
            hold on;
            s1=scatter(sucrosexy(:,1),sucrosexy(:,2),dotsize,[.4 .4 .4],'filled');
            s1.MarkerFaceAlpha = opacity;
            sc=scatter(xtrial_inc(post_stim==0),ytrial_inc(post_stim==0),'k','x');


            axis([0 950 75 750]);
            set(gca,'Ydir','reverse')
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            plot([950 portcoord(rat,1)],[450 portcoord(rat,2)],'color','k');
            legend(sc,'Location at cue onset');
            text(955,450,'Reward Port');
            title('No laser');


            %maltodextrin traces
            subplot(1,4,2);
            hold on;
            s2=scatter(maltodextrinxy(:,1),maltodextrinxy(:,2),dotsize,[0 .3 1],'filled');
            s2.MarkerFaceAlpha = opacity;
            scatter(xtrial_inc(post_stim==1),ytrial_inc(post_stim==1),[],'k','x')

            %plot individual trial traces
            %             for i=1:length(maltodextrintrls)
            %                 if maltodextrintrls{i,1}
            %                     plot(maltodextrintrls{i,1}(:,1),maltodextrintrls{i,1}(:,2),'color',Colors('maltodextrin'));
            %                 end
            %             end

            set(gca,'Ydir','reverse')
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            axis([0 950 75 750]);
            plot([950 portcoord(rat,1)],[450 portcoord(rat,2)],'color','k');
            legend(sc,'Location at cue onset');
            text(955,450,'Reward Port');
            title('Laser');

        end




        ITI_dfp_inc=ITI_dfp(included_cue_trials)';



        mean_distance(rat,1)=nanmean(ITI_dfp_inc(post_stim==0));
        mean_distance(rat,2)=nanmean(ITI_dfp_inc(post_stim==1));



        completed_trials(rat,1)=length(sucrosedeltime);
    end

    disp(['File #' num2str(k)]);
end


%% plotting measures

%included rats
yfp=[4 6 10 11 13 14 21];
chr=[1 2 8 9 16 17 18 22 23 24];

colors{2}=[0 0.3 1];
colors{1}=[0.4 0.4 0.4];
groups{1,1}=yfp;
groups{2,1}=chr;


subplot(1,8,5);
hold on;
plot([1 2],mean_distance(yfp,:),'color',[0.4 0.4 0.4]);
errorbar(1,nanmean(mean_distance(yfp,1)),nanste(mean_distance(yfp,1),1),'color',[0.4 0.4 0.4],'marker','o','linewidth',1.5);
errorbar(2,nanmean(mean_distance(yfp,2)),nanste(mean_distance(yfp,2),1),'color',[0 0.3 1],'marker','o','linewidth',1.5);
ylabel('Distance from port during ITI (pixels)');
xticks([1 2]);
xticklabels({'No laser','Laser'});
xtickangle(45);
axis([0.5 2.5 0 400]);
title('YFP');

pval=signrank(mean_distance(yfp,1),mean_distance(yfp,2));
text(1.2,375,['p = ' num2str(round(pval,2,'significant'))]);


subplot(1,8,6);
hold on;
plot([1 2],mean_distance(chr,:),'color',[0.4 0.4 0.4]);
errorbar(1,nanmean(mean_distance(chr,1)),nanste(mean_distance(chr,1),1),'color',[0.4 0.4 0.4],'marker','o','markerfacecolor',[0.4 0.4 0.4],'linewidth',1.5);
errorbar(2,nanmean(mean_distance(chr,2)),nanste(mean_distance(chr,2),1),'color',[0 0.3 1],'marker','o','markerfacecolor',[0 0.3 1],'linewidth',1.5);
ylabel('Distance from port during ITI (pixels)');
xticks([1 2]);
xticklabels({'No laser','Laser'});
xtickangle(45);
axis([0.5 2.5 0 400]);
title('ChR2');

pval=signrank(mean_distance(chr,1),mean_distance(chr,2));
text(1.2,375,['p = ' num2str(round(pval,2,'significant'))]);

included=zeros(24,1);
included(chr)=1;
included(yfp)=1;
arclog=zeros(24,1);
arclog(chr)=1;

subplot(1,4,4);
hold on;
difference=(mean_distance(:,2)-mean_distance(:,1))./mean_distance(:,1);


boxplot(difference(included==1),arclog(included==1),'colorgroup',[0 1],'colors',[colors{2};colors{2}],'boxstyle','filled');
scatter(rand([length(yfp) 1])/4+1.2,difference(yfp),55,colors{2});
scatter(rand([length(chr) 1])/4+2.2,difference(chr),55,colors{2},'filled');
axis([0.7 2.74 -0.5 0.5]);
plot([0 3],[0 0],'color','k');
pval=ranksum(difference(yfp),difference(chr));
text(1.2,0.4,['p = ' num2str(round(pval,2,'significant'))]);
ylabel('Fold-change over no laser');
xticks([1.22 2.22]);
xticklabels({'YFP','ChR2'});


dataforanova=[mean_distance(chr,1) mean_distance(chr,2);mean_distance(yfp,1) mean_distance(yfp,2)];