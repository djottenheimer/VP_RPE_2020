function Eind=myfind(EVENTS,Ename)
Eind=[];
if isfield(EVENTS{1},'name')
    for k=1:length(EVENTS)
        
        if strcmpi(EVENTS{k}.name,Ename)
            Eind=[Eind,k];
        end
        
    end
elseif iscell(EVENTS)
    for k=1:length(EVENTS)
        
        if strcmpi(EVENTS{k},Ename)
            Eind=[Eind,k];
        end
        
    end
    
end
end