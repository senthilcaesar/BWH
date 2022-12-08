[BreathDataTable2,~,BreathFLDataTable2,~,BreathDataTableFulls]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
      
%%
dt = 1/settings.Fs;
ploton=1;
EvtEnd = Evts.RespT.EventEnd;
EvtStart = Evts.RespT.EventStart;
EvStIdx = nan(length(EvtStart),1);
EvEndIdx = nan(length(EvtStart),1);
Window = nan(length(EvtStart),1);

for ii=1:length(EvtStart)
    if EvtStart(ii) < BreathDataTable2.Time_start(1) ||... %if outside whole table then skip
            EvtStart(ii) > BreathDataTable2.Time_start(end)
        excluderow = ii;
        continue
    end
    
    % find event start index and time
    TimeDiff = BreathDataTable2.Time_start - EvtStart(ii);
    TimeDiff(TimeDiff<0 | BreathDataTable2.FDuplicated2~=0) = nan;
    [~,EvStIdx(ii,1)] = min(TimeDiff);
    
    Windowtemp = BreathDataTable2.Time0(EvStIdx(ii,1));
    Window(ii)=Windowtemp;
    
    % find event end index and time - FIX THIS NEEDS TO BE IN SAME
    % WINDOW AS ABOVE
    TimeDiff = BreathDataTable2.Time_end - EvtEnd(ii);
    TimeDiff(TimeDiff>0 | BreathDataTable2.Time0 ~= Windowtemp) = nan;
    [~,EvEndIdx(ii,1)] = max(TimeDiff);
    if isnan(max(TimeDiff))
        EvEndIdx(ii,1)=NaN;
    end
        
end

EvEndIdx(isnan(EvEndIdx))=[];



%% Skip the unstretched stack method
Nbreaths=0;
for M=1:(Nbreaths*2+1)
    tempT = BreathDataTable2(EvEndIdx-Nbreaths+M-1,:);
    Tbreath = nanmedian(tempT.Ttot)
    Tbreathi = round(Tbreath/dt)
    
    invertflow=0; %fix
    
    SigIn = ((-1)^settings.Pnasaldownisinsp)*((-1)^invertflow)*DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %%% added invert flow
    %SigIn = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1));
    BreathBox_u{M}=nan(height(tempT),Tbreathi);
    
    figure(98); clf(98)
    for i=1:height(tempT)
        
        Time_starti=round((tempT.Time_start(i)-DataEventHypnog_Mat(1,1))/dt)+1;
        I = Time_starti:(Time_starti+round(tempT.Ttot(i)/dt)-1);
        %I = Time_starti:(Time_starti+round(tempT.Ti(i)/dt)-1);
        
        if 1
            exponent = settings.scalingexponent;
            leak = tempT.leak(i);
            leak2 = tempT.leak2(i);
            IEratio = tempT.IEratio(i);
            Vflow_out2 = SigIn(I)-leak;
            Vflow_out2(Vflow_out2>0)=(Vflow_out2(Vflow_out2>0).^(exponent))/(IEratio^0.5);
            Vflow_out2(Vflow_out2<0)=(-(-Vflow_out2(Vflow_out2<0)).^(exponent))*(IEratio^0.5);
            Vflow_out2 = Vflow_out2-leak2;
            Vflow_out2 = Vflow_out2./tempT.Veup(i);
            sigout = Vflow_out2;
        else
            sigout=SigIn(I);
        end
        
        sigout((Tbreathi+1):end)=[];
        ploton=1;
        if ploton
            figure(98)
            subplot(3,1,1);
            plot(sigout)
            hold on
        end
        BreathBox_u{M}(i,1:length(sigout))=(sigout(:))';
        
    end
    
    for jj=1:5
    if ploton
        figure(98)
        subplot(3,1,2);
        plot(nanmean(BreathBox{M}))
    end
    
    
        temp = BreathBox_u{M};
        temp2 = temp - nanmedian(temp);
        %temp2 = temp - prctile(abs(temp),25);
        %temp3 = nanmedian(abs(temp2),2);
        temp3 = nanmean(abs(temp2),2);
        temp3 = max(abs(temp2)')';
        BreathBox{M}(temp3>prctile(temp3,90),:)=NaN;
    
    
    if ploton
       figure(98)
        % plot(nanmean(BreathBox{M}));
        subplot(3,1,3);
        plot(BreathBox_u{M}');
        if 0
            pause(0.5)
        end
    end
    
    end
end

%
figure(80)
sig=[];
for M=1:(Nbreaths*2+1)
    sig = [sig nanmean(BreathBox_u{M})];
    plot(sig);
end

%% stretch and stack method

Nbreaths=5;
for M=1:(Nbreaths*2+1)
    tempT = BreathDataTable2(EvEndIdx-Nbreaths+M-1,:);
    Tbreath = nanmedian(tempT.Ttot)
    Tbreathi = round(Tbreath/dt)
    
    TbreathInsp = nanmedian(tempT.Ti)
    TbreathInspi = round(TbreathInsp/dt)
    
    TbreathExp = nanmedian(tempT.Te)
    TbreathExpi = round(TbreathExp/dt)
    
    invertflow=0; %fix
%     TbreathInspi + TbreathExpi
    
    SigIn = ((-1)^settings.Pnasaldownisinsp)*((-1)^invertflow)*DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %%% added invert flow
    %SigIn = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1));
    BreathBox{M}=nan(height(tempT),TbreathInspi+TbreathExpi);
    
    figure(98); clf(98)
    for i=1:height(tempT)
        
        Time_starti=round((tempT.Time_start(i)-DataEventHypnog_Mat(1,1))/dt)+1;
        I = Time_starti:(Time_starti+round(tempT.Ttot(i)/dt)-1);
        %I = Time_starti:(Time_starti+round(tempT.Ti(i)/dt)-1);
        
        if 1
            exponent = settings.scalingexponent;
            leak = tempT.leak(i);
            leak2 = tempT.leak2(i);
            IEratio = tempT.IEratio(i);
            Vflow_out2 = SigIn(I)-leak;
            Vflow_out2(Vflow_out2>0)=(Vflow_out2(Vflow_out2>0).^(exponent))/(IEratio^0.5);
            Vflow_out2(Vflow_out2<0)=(-(-Vflow_out2(Vflow_out2<0)).^(exponent))*(IEratio^0.5);
            Vflow_out2 = Vflow_out2-leak2;
            Vflow_out2 = Vflow_out2./tempT.Veup(i);
            sigout = Vflow_out2;
        else
            sigout=SigIn(I);
        end
        

        Xin2 = 0:1:([tempT.Ti(i)/dt]-1);
         Xin = Xin2 / (([tempT.Ti(i)/dt-1])/(TbreathInspi-1));
          Xin3 = 0:1:(TbreathInspi-1);
%         
        sigout2 = interp1(Xin,sigout(1:length(Xin)),Xin3);
        
        Xin2 = 0:1:(tempT.Te(i)/dt-1);
         Xin = Xin2 / (([tempT.Te(i)/dt-1])/(TbreathExpi-1));
          Xin3 = 0:1:(TbreathExpi-1);
%         
        sigout3 = interp1(Xin,sigout(end-length(Xin)+1:end),Xin3);
        
        
        
        sigoutS = [sigout2 sigout3];
        ploton=1;
        if ploton
            figure(98)
            subplot(3,1,1);
            %plot(sigout)
            
            plot([sigout2 sigout3])
            hold on
        end
        BreathBox{M}(i,1:length(sigoutS))=(sigoutS(:))';
        
    end
    
    for jj=1:3
        
    if ploton
        figure(98)
        subplot(3,1,2);
        plot(nanmean(BreathBox{M}))
        hold on
    end
    
    
        temp = BreathBox{M};
        temp2 = temp - nanmean(temp);
        %temp2 = temp - prctile(abs(temp),25);
        %temp3 = nanmedian(abs(temp2),2);
        temp3 = nanmean(abs(temp2),2);
        temp3 = nanmean(temp2.^2,2);
        temp3 = max(abs(temp2)')';
        BreathBox{M}(temp3>prctile(temp3,50),:)=NaN;
    
    
    if ploton
       figure(98)
        % plot(nanmean(BreathBox{M}));
        subplot(3,1,3);
        plot(BreathBox{M}');
        if 0
            pause(0.5)
        end
    end
    
    end
end

%
figure(80)
sig=[];
for M=1:(Nbreaths*2+1)
%     sig_breath = nanmean(BreathBox{M}); % temp hack, can be removed
%     sig_breath(isnan(sig_breath))=[];
    sig = [sig nanmean(BreathBox{M})];
    plot(sig);
end

