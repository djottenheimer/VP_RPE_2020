function [All,First,Last,Custom]=ChooseNs(R,DOI,ROI)

%a code to get "early" neurons and "late" neurons and "all" where no channel is
%repeated across multiple days for any rat
alls=DOI; %which days to look at

Rat={};
Day=[];
Channel=[];
%breaking it down by session
for i=1:length(R.Ninfo)
    A=char(R.Ninfo(i,1));
    Rat(i,1)=cellstr(A(1:3));
    RatNum(i,1)=str2num(A(2:3));
    Day(i,1)=str2num(A(8:9));
    B=char(R.Ninfo(i,2));
    Channel(i,1)=str2num(B(4:5));
end
earlies=min(Day);

C=unique(Rat(:,1));
D=unique(Day(:,1));
E=unique(Channel(:,1));

%count number of neurons on each channel
for i=1:length(C)
    Subj=strcmp(C(i),Rat);
    for j=D(1):D(end)
%        if R.CueRatio(i,j)>0.6 && R.CueRatio(i,j)<0.8
            for k=1:length(E)
                neurons{i,1}(k,j)=sum(Subj&Day==j&Channel==k);
            end
%         else
%             for k=1:length(E)
%                 neurons{i,1}(k,j)=0;
%             end
%         end
    end
end

%pick the session with most neurons for each channel 
for i=1:length(C)
    lates(i,1)=max(Day(strcmp(C(i),Rat)));
    for k=1:length(E)
        [~,earlysess{i,1}(k,1)]=max(neurons{i,1}(k,earlies));
        [~,latesess{i,1}(k,1)]=max(neurons{i,1}(k,lates(i,1)));
        [~,allsess{i,1}(k,1)]=max(neurons{i,1}(k,alls));
    end
end

%reset the logical
Ear=strcmp('nothing',R.Ninfo(:,1)); %this just makes all zeros
Lat=strcmp('nothing',R.Ninfo(:,1));
All=strcmp('nothing',R.Ninfo(:,1));

%create the selection
for i=1:length(C)
    Subj=strcmp(C(i),Rat);
    for k=1:length(E)
        Ear=Ear|(Subj&Channel==k&Day==(earlies(1)-1+earlysess{i,1}(k,1)));
        Lat=Lat|(Subj&Channel==k&Day==(lates(i,1)-1+latesess{i,1}(k,1)));
        All=All|(Subj&Channel==k&Day==(alls(1)-1+allsess{i,1}(k,1)));
    end
end

%which rats
rats=RatNum<0;
for i=1:length(ROI)
    rats=rats | RatNum==ROI(i);
end

All=All & rats;
First=Ear & rats;
Last=Lat & rats;
