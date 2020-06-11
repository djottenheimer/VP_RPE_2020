function  [PSRout Nout]=MakePSR04(Nrast,Erast,Win,PRMS)
%Win is the time window around the event
%PRMS=2 give a PSR cell {neurons,trials} and PRMS=1 collapses all the trials >
%cell{neurons,1}
%N is the number of trials



Otype=PRMS{1};
FLAG=0;Bmax=0;
if length(PRMS)>1
    FLAG=1;
end

if Otype==1
    PSRout=cell(size(Nrast,1),1);
    for i=1:length(Nrast)
        
        for k=1:size(Erast,1)
            tmp=Nrast{i}-Erast(k);
            if FLAG==0
                PSRout{i}=cat(1,PSRout{i},tmp(tmp>Win(1)-Bmax & tmp<Win(2)+Bmax));
            else
                PSRout{i}=cat(1,PSRout{i},tmp(find(tmp>Win(1)-Bmax & tmp<Win(2)+Bmax,1,PRMS{2})));
            end
        end
        
    end
elseif Otype==2
    
    PSRout=cell(size(Nrast,1),size(Erast,1)); %PSRout is a cell(neurons,trials) 
    for i=1:length(Nrast) %loops through the number of neurons
        
        for k=1:size(Erast,1) %loops through the trials
            tmp=Nrast{i}-Erast(k);
            if ~isempty(find(tmp>Win(1) & tmp<Win(2), 1))
                if FLAG==0
                    PSRout{i,k}=tmp(tmp>Win(1)-Bmax & tmp<Win(2)+Bmax);
                else
                    PSRout{i,k}=tmp(find(tmp>Win(1)-Bmax & tmp<Win(2)+Bmax,1,PRMS{2}));
                end
            else
                PSRout{i,k}=NaN;
            end
        end
        
    end
end
Nout=size(Erast,1);