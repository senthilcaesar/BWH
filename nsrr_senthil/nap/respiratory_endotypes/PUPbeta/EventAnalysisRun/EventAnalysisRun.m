function [Boxes,Ensembles]=EventAnalysisRun(SignalsT,BreathDataTable2,RespT,Method,RespList)
%% Used exclusively as a box and ensemble generator
%no models etc
%[Boxes,Ensembles]=EventAnalysisRun_20190710(SignalsT,BreathDataTable2,Fs(?),Evts)
%%current plan is to handle both data in signal time plus optional for breath time
%%BreathDataTable2 should contain Table subset to be used to get boxes
%%(RespList is then table col headings)
global settings

%%then should be able to run in Analysis before EventRespTable to
%%generate signal (not breath level) boxes for Desats and Arousals minimum)
if ~exist('RespList') || isempty(RespList)
    %Add breath-domain variables to this list that may ever be used
    RespList = {'VI','DeltaEdi','PeakEdi','VdriveEdiNorm','GGpeak','GGtonic', ...
        'DeltaPes','VdrivePesNorm','DeltaThorax','DeltaAb','DeltaPmus','VdrivePmusNorm','FlowDrive',...
        'VdriveEst','Time_start','Time_mid','Time_end','VT','VdrivePdriveKNorm','VdrivePdrive30Norm',...
        'PdriveK','Pdrive30','Pdrive10','PdriveKraw','VdrivePdriveKrawNorm','DeltaPepi','VdrivePepiNorm','AA_ExpFlat90Te_O','Tpeak1_Ti_T',...
        'YFOTmeanI','YFOTmeanIW','YFOTmeanE','YFOTmeanEW'};
end

dT=1; %0.25
rangeX=100;
if ~exist('Method')
    Method=1;
end

%% Misc
%     EvtData.TST=length(find(SignalsT.Epochs>=0 & SignalsT.Epochs<4))/(Fs*60);
%     EvtData.TST_REM=length(find(SignalsT.Epochs==3))/(Fs*60);
%     EvtData.TST_NREM=length(find(SignalsT.Epochs>=0 & SignalsT.Epochs<3))/(Fs*60);
%     EvtData.TIB=length(SignalsT.Epochs)/(Fs*60);

%% Prepare Drive data (move outside)
%     try
%         if 1
%             [G,~] = VdriveFromEdi(BreathDataTable2,'VE','DeltaEdi');
%             BreathDataTable2.VdriveEdi = BreathDataTable2.DeltaEdi.*G; %L/min
%             BreathDataTable2.VdriveEdiNorm = BreathDataTable2.VdriveEdi./BreathDataTable2.Veup; %Fraction of eupnea
%         end
%     catch
%     end
%     try
%         if 1
%             [G,~] = VdriveFromEdi(BreathDataTable2,'VE','DeltaPes');
%             BreathDataTable2.VdrivePes = BreathDataTable2.DeltaPes.*G; %L/min
%             BreathDataTable2.VdrivePesNorm = BreathDataTable2.VdrivePes./BreathDataTable2.Veup; %Fraction of eupnea
%         end
%     catch
%     end

%% New way to event start and stop times
%     evtStart_i = Evts.RespT.starttimesi;
EvtStart = RespT.EventStart;
EvtEnd = RespT.EventEnd;
%     evtEnd_i = round(Evts.RespT.EventEnd*Fs);

if  isfield(settings,'EventAnalysisForArousal') && settings.EventAnalysisForArousal==1
    disp('SS: sorry I broke your code. Try passing a modified version of Evts.Table1 in place of the RespT table directly')
%     if  isfield(settings,'EventAnalysisForArousal') & settings.EventAnalysisForArousal==1
%         EvtEnd = Evts.Table1.EventStart(strcmp(Evts.Table1.EventCodesList,'AR'));
%         EvtStart = EvtEnd-30;
%     end
end

%% Resp Event Analysis using SignalsT
TimeEnsemble=-rangeX:dT:rangeX;
%     Ydata=BreathDataTable2.VI;
%     Pesdata=BreathDataTable2.DeltaPes;

%     TotalAHI=60*length(EvtEnd)/EvtData.TST;
%
%     EvtData.TotalAHI=TotalAHI;
%     EvtData.Start=EvtStart-(EvtEnd-rangeX);
%     EvtData.End=EvtEnd-(EvtEnd-rangeX);
%     EvtData.SyncTime=EvtEnd;

Boxes.Time=repmat(TimeEnsemble,length(EvtEnd),1)+repmat(EvtEnd,1,length(TimeEnsemble));

if ~isempty(SignalsT)
%edit / add additional signals to Table
try
%     EventsResp=0*SignalsT.Time; % MOVED TO EVENT ANALYSIS
%     for m=1:length(RespT.EventStart) %Resp Events
%         lefti=round((RespT.EventStart(m)-SignalsT.Time(1))*settings.Fs)+1;
%         righti=lefti+round((RespT.EventDuration(m))*settings.Fs);
%         if lefti<1, lefti=1; end
%         if righti>length(SignalsT.Time), righti=length(SignalsT.Time); end
%         if RespT.EventCodes(m)>1
%             EventsResp(lefti:righti)=RespT.EventCodes(m);
%         end
%     end
%     SignalsT.EventsResp = 1*(EventsResp>0);
end

List = SignalsT.Properties.VariableNames;
Boxes.Time=repmat(TimeEnsemble,length(EvtEnd),1)+repmat(EvtEnd,1,length(TimeEnsemble));
for i=1:length(List)
    try
        temp=interp1(SignalsT.Time,SignalsT{:,List{i}},Boxes.Time,'nearest');
        Boxes.(List{i}) = temp;
    catch
    end
end
end

%% Resp. Event Analysis using BreathDataTable
try
    %         BrTime=BreathDataTable2.Time_start;
    %         EType=BreathDataTable2.Etype;
    %         EType(EType>0)=1;
    %         EType(end)=0;
    %         IdxUp=find([0; diff(EType)>0]);
    %         IdxDn=find([diff(EType)<0; 0]);
    %
    %         EvStIdx2=IdxUp(BreathDataTable2.FDuplicated2(IdxUp)==0);
    %         EvIdx=IdxDn(BreathDataTable2.FDuplicated2(IdxUp)==0);
    %         Window2=BreathDataTable2.Time0(EvIdx);
    %         evtEnd=BreathDataTable2.Time_end(EvIdx);
    %         evtSt=BreathDataTable2.Time_start(EvStIdx);
    %
    %         EvtSig=zeros(size(SignalsT.Time));
    %         for ii=1:length(EvIdx)
    %            EvtSig(SignalsT.Time>=evtSt(ii) & SignalsT.Time<=evtEnd(ii))=1;
    %         end
    
    
    % Make event start and stop indices, lookup in BreathDataTable2
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
    end
    
    % remove nans
    EvStIdx(isnan(EvStIdx)) = [];
    EvEndIdx(isnan(EvEndIdx)) = [];
    
    % Get associated time variables (e.g. event start time)
    EvtStart2 = BreathDataTable2.Time_start(EvStIdx);
    EvtEnd2 = BreathDataTable2.Time_end(EvEndIdx);
    %Window = BreathDataTable2.Time0(EvStIdx);
    
    
    %RespList = BreathDataTable2.Properties.VariableNames;
    for i=1:length(RespList)
        try
            tempa = NaN*ones(length(EvtStart),length(TimeEnsemble)); %initialize with NaN:
            Boxes.(RespList{i}) = tempa;
            for ii=1:length(EvtStart)
                try
                    if ~isnan(Window(ii))
                    I=(BreathDataTable2.Time0==Window(ii));
                    temp1 = BreathDataTable2{I,RespList{i}};
                    else 
                        continue
                    end
               catch
                   continue
                end
                switch Method
                    case 1
                        Time_temp=BreathDataTable2.Time_end(I)-EvtEnd(ii); % what do we want to use for EvtEnd??
                        temp2 = interp1(Time_temp,temp1,TimeEnsemble,'linear');
                    case 2
                        Time_temp=BreathDataTable2.Time_start(I)-EvtEnd(ii); % what do we want to use for EvtEnd??
                        temp2 = interp1(Time_temp,temp1,TimeEnsemble,'linear');
                end               
                Boxes.(RespList{i})(ii,:) = temp2;
            end
        catch
            disp([RespList{i}, ' is not available from BreathDataTable'])
        end
    end
    
    for i=1:length(RespList)
        try
            %initialize with NaN:
            tempa = NaN*ones(length(EvtStart),length(TimeEnsemble));
            Boxes.([RespList{i} 'Stairs']) = tempa;
            for ii=1:length(EvtStart)
                try
                    if ~isnan(Window(ii))
                    I=(BreathDataTable2.Time0==Window(ii));
                    temp1 = BreathDataTable2{I,RespList{i}};
                    else 
                        continue
                    end
               catch
                   continue
                end
                Time_temp=BreathDataTable2.Time_start(I)-EvtEnd(ii); % what do we want to use for EvtEnd??
                temp2 = interp1(Time_temp,temp1,TimeEnsemble,'previous');
                Boxes.([RespList{i} 'Stairs'])(ii,:) = temp2;
            end
        catch
            disp([[RespList{i} 'Stairs'], ' is not available from BreathDataTable'])
        end
    end
    
catch
    disp('BreathDataTable does not exist')
end

try %move?
Boxes.Ti=Boxes.Time_mid-Boxes.Time_start;
Boxes.Te=Boxes.Time_end-Boxes.Time_mid;
Boxes.Ttot=Boxes.Time_end-Boxes.Time_start;
end

clear temp temp1 temp2 Time_temp temprs tempa

Ensembles=[];
BoxesList=fieldnames(Boxes);
for i=1:length(BoxesList)
    temp=nanmean(Boxes.(BoxesList{i}),1); %SignalsT.Time,EvtSig,Boxes.Time,'nearest');
    Ensembles.(BoxesList{i}) = temp(:);
end
Ensembles.Time = TimeEnsemble(:);
Ensembles = struct2table(Ensembles);

%%
%     [Ftrs]=EventFeaturesRun(Boxes,Ensembles,Evts.RespT,TotalAHI,Models);

end

