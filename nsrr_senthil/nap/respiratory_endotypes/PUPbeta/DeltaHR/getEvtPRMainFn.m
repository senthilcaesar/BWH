function Evts=getEvtPRMainFn(Evts,SigT)

minpossibleHR=20;
minvalidsignaldur=300;

if ~isfield(Evts,'RespT')
    Evts.RespT = EventRespSigtoRespT(SigT.EventsResp,SigT.Time);
end
if ~isfield(Evts,'ArT')
    Evts.ArT = EventRespSigtoRespT(SigT.EventsAr,SigT.Time);
end
dt=SigT.Time(2)-SigT.Time(1);

figure(23); clf(23); set(gcf,'color',[1 1 1]);
h=subplot(2,1,1);

try
    HRtemp = SigT.HR;
    dt = SigT.Time(2)-SigT.Time(1);
    ww = 6;
%     SigT.HRfilt = medfilt1(HRtemp,round(ww/dt));
    SigT.HRfilt = movmedian(HRtemp,round(ww/dt),'omitnan');
    
    %start
    HRsignal = 'HRfilt';
    LocalSignalsT = SigT(:,{'Time',HRsignal});
    
    temp = LocalSignalsT{:,2};
    temp(temp<minpossibleHR)=NaN;
    LocalSignalsT{:,2} = temp;
 
    ArMode=0;
    
    validDur = sum(~isnan(LocalSignalsT{:,2}))*dt;
    if validDur<minvalidsignaldur
        h=[];
    end
    
    HRresponse = HRresponseFromSignalT(LocalSignalsT,Evts.RespT,HRsignal,ArMode,h);
    nanmean(HRresponse.DeltaFromMin)
    HRresponse.Properties.VariableNames = strcat('HR',HRresponse.Properties.VariableNames);

    
    I = 1==sum(Evts.RespT.Properties.VariableNames==string(HRresponse.Properties.VariableNames'));
    Evts.RespT(:,I)=[];
    Evts.RespT = [Evts.RespT HRresponse];
catch
    disp('No EKG heart rate data available')
end

try
    HRsignal = 'Pulse';
    hold on
    LocalSignalsT = SigT(:,{'Time',HRsignal});
    
    temp = LocalSignalsT{:,2};
    temp(temp<minpossibleHR)=NaN;
    LocalSignalsT{:,2} = temp;

    ArMode=0;
    
    validDur = sum(~isnan(LocalSignalsT{:,2}))*dt;
    if validDur<minvalidsignaldur
        h=[];
    end
    
    PRresponse = HRresponseFromSignalT(LocalSignalsT,Evts.RespT,HRsignal,ArMode,h);
    nanmean(PRresponse.DeltaFromMin)
    PRresponse.Properties.VariableNames = strcat('PR',PRresponse.Properties.VariableNames);
    

    I = 1==sum(Evts.RespT.Properties.VariableNames==string(PRresponse.Properties.VariableNames'));
    Evts.RespT(:,I)=[];
    Evts.RespT = [Evts.RespT PRresponse];
catch
    disp('No pulse rate data available')
end

HRsignal = 'HRfilt';
h=subplot(2,1,2);
try
    LocalSignalsT = SigT(:,{'Time',HRsignal});
    temp = LocalSignalsT{:,2};
    temp(temp<minpossibleHR)=NaN;
    LocalSignalsT{:,2} = temp;
    
    validDur = sum(~isnan(LocalSignalsT{:,2}))*dt;
    if validDur<minvalidsignaldur
        h=[];
    end
    
    ArMode=1;
    HRresponseAr = HRresponseFromSignalT(LocalSignalsT,Evts.ArT,HRsignal,ArMode,h);
    nanmean(HRresponseAr.DeltaFromMin)
    HRresponseAr.Properties.VariableNames = strcat('HR',HRresponseAr.Properties.VariableNames);
    
    I = 1==sum(Evts.ArT.Properties.VariableNames==string(HRresponseAr.Properties.VariableNames'));
    Evts.ArT(:,I)=[];
    Evts.ArT = [Evts.ArT HRresponseAr];
catch
end

try
    HRsignal = 'Pulse';
    hold on
    LocalSignalsT = SigT(:,{'Time',HRsignal});
    temp = LocalSignalsT{:,2};
    temp(temp<minpossibleHR)=NaN;
    LocalSignalsT{:,2} = temp;
    
    ArMode=1;
    
    validDur = sum(~isnan(LocalSignalsT{:,2}))*dt;
    if validDur<minvalidsignaldur
        h=[];
    end
    
    PRresponseAr = HRresponseFromSignalT(LocalSignalsT,Evts.ArT,HRsignal,ArMode,h);
    nanmean(PRresponseAr.DeltaFromMin)
    PRresponseAr.Properties.VariableNames = strcat('PR',PRresponseAr.Properties.VariableNames);

    
    I = 1==sum(Evts.ArT.Properties.VariableNames==string(PRresponseAr.Properties.VariableNames'));
    Evts.ArT(:,I)=[];
    Evts.ArT = [Evts.ArT PRresponseAr];
catch
end