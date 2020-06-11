function [PTHout,BWout,MSRout]=MakePTH06(PSRin,Numin,PRMS)
%!!!!!!!!CAUTION !!!!!
%this version was modified by Fred and Ali to use DP and NOT Optimal BW
%
%INPUTS:
%PSRin can be collapsed  or not. if not, the output is uncollapsed as well
%Numin= number of trials
%PRMS{1}==1 for making baseline-only PSTH and ==2 for normal, main PSTH %(i.e. -15 to 15s)
%PRMS{2}==0 for fixed binwidth and ==1 for DP binwidth
%PRMS{3}=binwidth used
%OUTPUTS:
%PTHout
%BWout is a vector stores the DP BW
%MSRout is MSR=f(BW)




global Tm Tbase

if PRMS{1}==1
    T0=Tbase;
elseif PRMS{1}==2
    T0=Tm;
end
BWx=[0.01:0.005:1]; % Range of binsizes considered


PTH=[];
if size(PSRin,2)==1 %if TRUE, all trials are collapsed
    if PRMS{2}==0 %fixed/assigned binwidth for all
        if numel(PRMS{3})==1,
            BW=repmat(PRMS{3},length(PSRin),1);
        else
            BW=PRMS{3};
        end
        for k=1:length(PSRin)
            Tn=[-BW(k)/2:-BW(k):T0(1)-BW(k)]; Tp=[BW(k)/2:BW(k):T0(end)+BW(k)];
            T=[Tn(end:-1:1),Tp];%making PSTH centered on zero
            %T=[T0(1)-BW(k):BW(k):T0(end)+BW(k)];
            TMP=hist(PSRin{k},T)/(BW(k)*Numin(k)); %  make PSTH and normalize by number of refs used (Hz)
            %deals with low firing in which there is too few DP bins to do
            %interppolation
            if length(T)>4
                PTH(k,:)=interp1(T(2:end-1),TMP(2:end-1),T0,'nearest','extrap');BWout(k,1)=BW(k);
            else
                PTH(k,:)=repmat(mean(TMP),1,length(T0));BWout(k,1)=T0(end)-T0(1);
            end
            
        end
        PTHout=PTH;MSRout=[];
    else %DP binwidth for all
        for k=1:length(PSRin)
            [~,~,MSR]=optimalBW05(PSRin(k),T0); %  make PSTH and normalize by number of refs used (Hz)
            [~,I]=findDP(BWx,MSR,{2,.9});%{DP=1,DEG/TSH=2,TSHOLD}
            Tn=[-BWx(I)/2:-BWx(I):T0(1)-BWx(I)]; Tp=[BWx(I)/2:BWx(I):T0(end)+BWx(I)];T=[Tn(end:-1:1),Tp];%making PSTH centered on zero
            
            TMP=hist(PSRin{k},T)/(BWx(I)*Numin(k)); %  make PSTH and normalize by number of refs used (Hz)
            %deals with low firing in which there is too few DP bins to do
            %interppolation
            if length(T)>4
                PTH(k,:)=interp1(T(2:end-1),TMP(2:end-1),T0,'nearest','extrap'); BWout(k,1)=BWx(I);
            else
                PTH(k,:)=repmat(mean(TMP),1,length(T0));BWout(k,1)=T0(end)-T0(1);
            end
            
           MSRout(k,:)=MSR;
        end
        PTHout=PTH;
    end
else % trial-by-trial PTH
    if PRMS{2}==0 %fixed/assigned binwidth for all
        if numel(PRMS{3})==1,
            BW=repmat(PRMS{3},length(PSRin),1);
        else
            BW=PRMS{3};
        end
        for k=1:size(PSRin,1) %neuron
            Tn=[-BW(k)/2:-BW(k):T0(1)-BW(k)]; Tp=[BW(k)/2:BW(k):T0(end)+BW(k)];T=[Tn(end:-1:1),Tp];%making PSTH centered on zero
            %  T=[T0(1)-BW(k):BW(k):T0(end)+BW(k)];
           
            for n=1:size(PSRin,2)%trial
                TMP=hist(PSRin{k,n},T)/BW(k); %  make PSTH and normalize by number of refs used (Hz)
                %deals with low firing in which there is too few DP bins to do
                %interppolation
                if length(T)>4
                    PTH(k,:,n)=interp1(T(2:end-1),TMP(2:end-1),T0,'nearest','extrap');BWout(k,1)=BW(k);
                else
                    PTH(k,:,n)=repmat(mean(TMP),1,length(T0));BWout(k,1)=T0(end)-T0(1);
                end
                
            end
        end
        PTHout=PTH;MSRout=[];
        
    else %DP binwidth for all
        for k=1:size(PSRin,1)%neuron
            [~,~,MSR]=optimalBW05(PSRin(k,:),T0);
            [~,I]=findDP(BWx,MSR,{2,.9});%{DP=1,DEG/TSH=2,TSHOLD}
            Tn=[-BWx(I)/2:-BWx(I):T0(1)-BWx(I)]; Tp=[BWx(I)/2:BWx(I):T0(end)+BWx(I)];T=[Tn(end:-1:1),Tp];%making PSTH centered on zero
            for n=1:size(PSRin,2)%trial
                TMP=hist(PSRin{k,n},T)/BWx(I); %  make PSTH and normalize by number of refs used (Hz)
                %deals with low firing in which there is too few DP bins to do
                %interppolation
                if length(T)>4
                    PTH(k,:,n)=interp1(T(2:end-1),TMP(2:end-1),T0,'nearest','extrap');BWout(k,1)=BWx(I);
                else
                    PTH(k,:,n)=repmat(mean(TMP),1,length(T0));BWout(k,1)=T0(end)-T0(1);
                end
                
            end
           MSRout(k,:)=MSR;
        end
        PTHout=PTH;
    end
    
end


