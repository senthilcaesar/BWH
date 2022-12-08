function [Evts]=EventRespTable(SigT,Evts,ChannelsList,EventsChannelName,ArChannelName)

global settings n
%%currently called from Analysis.m
%Will also be called from EventAnalysisAA 
%will input boxes. separate custom subroutines called from here based on existance
%of necsesary signal in Boxes
%will output events table same as now
%Inputs: Boxes,Evts

%Also keep structure variable, single valued, to store patient level data
%e.g. circ delay, time to start arousal, time to peak HR, event duration;
%name these starting with the signal e.g. SpO2Delay, SpO2Baseline, HRBaseline, 

if ~exist('EventsChannelName')
    EventsChannelName = 'EventsResp';
end

% if ~exist('UseEventsRespAuto')
%     UseEventsRespAuto=0; %not used yet
% end
% 
% if UseEventsRespAuto==1
%     EventsChannelName = 'EventsRespAuto';
% end



if exist('ArChannelName')==1 && ~isempty(ArChannelName)
     %do nothing
else
       
    if isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==1
        ArChannelName = 'EventsArWS';
    elseif isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==2
        ArChannelName = 'EventsArWSB';
    else
        ArChannelName = 'EventsAr';
    end
end

disp(EventsChannelName);
disp(ArChannelName);

plotfig=0;
if isfield(settings,'plotHBfigs') 
    plotfig=settings.plotHBfigs;
end

pausemanually=0; %debug mode
if pausemanually
    verbose=1;
    plotfig=1;
else
    verbose=0;
end

% addpath(genpath('E:\PUPbeta_git\PUPbeta20190629'));
Time = SigT.Time;
dt = Time(2) - Time(1);
Fs = 1/dt;

%% Get the events start and end times either from Evts table or events signal
% From table version
if 1
    if ~isfield(Evts,'RespT')
        T = Evts.Table1;
        I = T.EventCodes>1;
        T=T(I,:);
    else
        T = Evts.RespT;
    end
    
    T.EventEnd = T.EventStart + T.EventDuration;
    T = sortrows(T,{'EventEnd','EventStart'}); %added 1/18/2020 3:06 pm; data before this may have events out of order (based on recording system export).

% else
%     % From signal version
%     if Flag==1 % for using with automated scoring
%         RespCh = find(strcmp(ChannelsList,'EventsRespAuto')==1);
%     else % for using with clinical scoring from PSG
%         RespCh = find(strcmp(ChannelsList,'EventsResp')==1);
%     end
%     
%     EventStarti=find([diff(DataEventHypnog_Mat(:,RespCh))>0]);
%     EventEndi=find([diff(DataEventHypnog_Mat(:,RespCh))<0]);
%     EventDuration = (EventEndi - EventStarti)/Fs;
%     EventEndt = Time(EventEndi);
%     %     IdxUp=find([0; diff(EType)>0]);
end

Nevents = size(T.EventEnd,1);
if Nevents==0
    T.SpO2DeltaE = zeros(0);
    T.ArE = zeros(0);
    T.SpO2NadirE = zeros(0);
    T.SpO2NadirEt = zeros(0);
    T.SpO2BaselineE = zeros(0);
    T.SpO2BaselineEt = zeros(0);
    T.HBarea = zeros(0);
    T.HBbaseline = zeros(0);
    T.ArEt = zeros(0);
    T.ArDurationE = zeros(0);
    T.starttimesi = zeros(0);
    T.Epochs = zeros(0);
    T.EpochsMidEvent = zeros(0);
    T.Position = zeros(0);
    T.PositionRaw = zeros(0);
    Evts.RespT = T;
    return
end
%%


%% Get the ensemble average
rangeX=120;
dT=0.25; %0.25
TimeEnsemble=-rangeX:dT:rangeX;
TimeBox=repmat(TimeEnsemble,length(T.EventEnd),1)+repmat(T.EventEnd,1,length(TimeEnsemble));

SigList = {EventsChannelName,'SpO2',ArChannelName};
SigDir = [1 -1 1];
%% Remove SpO2 Artifacts (5 s)
if 1
    SigCh = find(strcmp(ChannelsList,'SpO2')==1);
    if ~isempty(SigCh)
    SpO2=SpO2ArtifactReject(SigT.SpO2,1/Fs);
    %replace SpO2 in DataEventHypnog_Mat
    SigT.SpO2=SpO2ArtifactReject(SigT.SpO2,1/Fs);
    end
end

%% Make Signal Box and Plot

Nplots = length(SigList);
if plotfig||1
        figure(76); clf(76);
end
for i=1:length(SigList)
    if ~isempty(find(strcmp(ChannelsList,SigList{i})==1))
        SigCh(i) = find(strcmp(ChannelsList,SigList{i})==1);
        SigTemp=SigT{:,SigCh(i)};
        Sig=interp1(Time,SigTemp,TimeBox,'nearest');
    else
        Sig=nan(size(TimeBox));
    end
        
        if strcmp(EventsChannelName,SigList{i})
            Sig = 1*(Sig>0);
        end
        if strcmp('SpO2',SigList{i})
            Sig(Sig<40)=NaN;
        end

    eval([SigList{i} 'Box=Sig;']);
    
    % SigEnsmb=nanmean(Sig,1); % ensemble averaged signal across time;
    try
        if plotfig||1
            subplot(Nplots,1,i);
            plot(TimeEnsemble,Sig,'color',0.8*[1 1 1])
            hold on
            box off
            plot(TimeEnsemble,nanmean(Sig),'k','linewidth',2)
        end
    catch me
        disp('');
    end
end

%% Test
if 1
    windowandfilter=1;
    
    if windowandfilter
    [B,A] = butter(1,((1/120))/((1/dT)/2),'high');
    win = hann(size(SpO2Box,2));
    else
    B=1; A=1;
    win = ones(size(SpO2Box,2),1);    
    end
    figure(22); clf(22);
    temp2 = SpO2Box-nanmean(SpO2Box,2);
    I = ones(size(SpO2Box,1),1)==1;
    temp = (nanmean(SpO2Box(I,:)));
    temp = temp-nanmean(temp);
    plot(TimeEnsemble,temp,'linewidth',2);
    hold on
    temp=nanfilter(B,A,temp(:),1)'; %replaces temp=filtfilt(B,A,temp);
    plot(TimeEnsemble,temp,'linewidth',2);

    if 1
        for i=1:2
            temp = (nanmean(SpO2Box(I,:),1));
            temp = temp - nanmean(temp); 
                temp=nanfilter(B,A,temp(:),1)'; %replaces temp=filtfilt(B,A,temp);
            temp = filter121(temp(:).*win,round(3/dT));
            plot(TimeEnsemble,temp,'linewidth',2);
            diffSpO2 = temp2 - temp(:)';
            Err=nanmean(abs(diffSpO2),2);
            thresX = prctile(Err,90);
            Evts.RespSearchInfo.ErrCentile=thresX;
            if thresX<1,thresX=1; end
            thres2 = max([2 prctile(Err,50)]);
            if thresX>thres2,thresX=thres2; end
            I = Err<thresX;
            Evts.RespSearchInfo.ErrCentile2=thresX;
            Evts.RespSearchInfo.FErr=1-mean(I);
            hold on
        end
        if 0
            figure(2); clf(2);
            hist(Err)
        end
    end
    I1 = I;
    if 1
        for i=1:3
            temp = (nanmean(SpO2Box(I,:),1));
            temp = temp - nanmean(temp); 
            temp=nanfilter(B,A,temp(:),1)'; %replaces temp=filtfilt(B,A,temp);
            temp = filter121(temp(:).*win,round(3/dT));
            clear CorrX
            for i=1:size(SpO2Box,1)
                CorrX(i,1)=corr(temp2(i,:)',temp(:),'rows','complete');
            end
            if 0
                figure(2); clf(2);
                hist(CorrX)
            end
            thresX = 0.5; 
            CorrX2=CorrX; CorrX2(I1==0)=-1; thresX2 = prctile(CorrX2,50);
            thresX = min([thresX thresX2]);
            thresX = max([thresX 0]);
            Evts.RespSearchInfo.FLowCorr=mean(CorrX<=thresX);
            I = CorrX>thresX & I;
            Evts.RespSearchInfo.FLowQ=1-mean(I);
        end
        
    end
    temp = (nanmean(SpO2Box(I,:),1));
            temp = temp - nanmean(temp); 
            temp=nanfilter(B,A,temp(:),1)'; %replaces temp=filtfilt(B,A,temp);
            plot(TimeEnsemble,temp,'linewidth',2);
            temp = filter121(temp(:).*win,round(3/dT));
    plot(TimeEnsemble,temp,'linewidth',2);
    Evts.RespSearchInfo.NboxesUsed = sum(I);
    Evts.RespSearchInfo.NboxesAvailable = length(I);
    
end
%%
if 1
    SpO2Box = round(SpO2Box);
end

%% try this
UsePeriodToDefineSearchRange=1;
if UsePeriodToDefineSearchRange
    [~,~,~,~,T1] = WaveformBreakdown(temp(:),TimeEnsemble(:),[],0,1);
    Evts.RespSearchInfo.Period = T1;
    uppersearchrange = max([20 min([T1*1.33 45])]);
    Evts.RespSearchInfo.uppersearchrange = uppersearchrange;
else
    uppersearchrange=45
end

%% Defaults
minNevents = 8;
delaySpO2 = 12;
MeanEventDuration = nanmean(T.EventDuration);
searchBaseline = -15; %N seconds before event start
searchNadir = 25; %N seconds after delay
searchAr = 15;
delayAr = 3;

%%
if sum(I)>minNevents%length(T.EventDuration)>minNevents
    %% Find Delay
    initialsearchranget=[0 uppersearchrange]; %for SpO2
    initialsearchrangei=round((initialsearchranget--rangeX)/dT+1);
    clear delay
    for i=1:length(SigList)
        SigCh = find(strcmp(ChannelsList,SigList{i})==1);
        SigDir(i);
        
        %TimeEnsemble(initialsearchrangei)
        tempSig = eval([SigList{i} 'Box']);
        if i==2
            tempSig = nanmean(tempSig(I,:));
            tempSig = filter121(tempSig,round(3/dT));
        else
           tempSig = nanmean(tempSig);
        end
        if 0
            tempSig = detrend(tempSig);
        end
        [maxv,maxi]=max(SigDir(i)*tempSig(initialsearchrangei(1):initialsearchrangei(2)));
        maxv=SigDir(i)*maxv;
        maxi=maxi + initialsearchrangei(1)-1;
        
        delay(i) = TimeEnsemble(maxi);
    end
    
    %% Find Left Search
    initialsearchranget=[delay(2)-MeanEventDuration*2 delay(2)]; %for SpO2
    initialsearchrangei=round((initialsearchranget--rangeX)/dT+1);
    clear baselinet
    for i=1:length(SigList)
        SigCh = find(strcmp(ChannelsList,SigList{i})==1);
        tempSig = eval([SigList{i} 'Box']);
        if i==2
            tempSig = nanmean(tempSig(I,:));
            tempSig = filter121(tempSig,round(3/dT));
        else
           tempSig = nanmean(tempSig);
        end
        [maxv,maxi]=min(SigDir(i)*tempSig(initialsearchrangei(1):initialsearchrangei(2)));
        maxv=SigDir(i)*maxv;
        maxi=maxi + initialsearchrangei(1)-1;
        baselinet(i) = TimeEnsemble(maxi);
    end
    
    %% Find Recovery
    initialsearchranget=[delay(2) delay(2)+30]; %for SpO2
    initialsearchrangei=round((initialsearchranget--rangeX)/dT+1);
    clear peakt
    for i=1:length(SigList)
        SigCh = find(strcmp(ChannelsList,SigList{i})==1);
        tempSig = eval([SigList{i} 'Box']);
        if i==2
            tempSig = nanmean(tempSig(I,:));
            tempSig = filter121(tempSig,round(3/dT));
        else
           tempSig = nanmean(tempSig);
        end
        [maxv,maxi]=min(SigDir(i)*tempSig(initialsearchrangei(1):initialsearchrangei(2)));
        maxv=SigDir(i)*maxv;
        maxi=maxi + initialsearchrangei(1)-1;
        peakt(i) = TimeEnsemble(maxi);
    end
    
    
    %% Find individual desats
    delaySpO2 = delay(2);
    
    DesatPeakTroughPeakTimes = [baselinet(2) delaySpO2 peakt(2)];
    SpO2atDesatPeakTroughPeakTimes = interp1(TimeEnsemble,nanmean(SpO2Box),DesatPeakTroughPeakTimes);
    
    SuspiciousSpO2Window = mean(abs(diff(SpO2atDesatPeakTroughPeakTimes)))<1.5;

    if 0&&SuspiciousSpO2Window %new
        disp('warning: SuspiciousSpO2Window');
        searchBaseline1 = -10 -(delaySpO2 - baselinet(2)) + MeanEventDuration; %N seconds before event start
        searchBaseline2 = -10; %N seconds before event start
        searchBaseline = min([searchBaseline1 searchBaseline2]);
    else
        searchBaseline = -10 -(delaySpO2 - baselinet(2)) + MeanEventDuration; %N seconds before event start    
    end
    
    %(delaySpO2 - baselinet(2)) is desat event duration
    searchNadir = peakt(2) - delaySpO2;
    delayAr = delay(3);
    
end
%% Find individual desats loop
if exist('SuspiciousSpO2Window')
    Evts.RespSearchInfo.SuspiciousSpO2Window = SuspiciousSpO2Window;
end
if exist('SpO2atDesatPeakTroughPeakTimes')
Evts.RespSearchInfo.SpO2atDesatPeakTroughPeakTimes = SpO2atDesatPeakTroughPeakTimes;
end
if exist('DesatPeakTroughPeakTimes')
Evts.RespSearchInfo.DesatPeakTroughPeakTimes = DesatPeakTroughPeakTimes;
end


delaySpO2i = round(delaySpO2/dT);
searchNadiri = round(searchNadir/dT);
toli = (2/dT);
SpO2BaselineE = nan(Nevents,1);
SpO2BaselineEt = nan(Nevents,1);
SpO2NadirE = nan(Nevents,1);
SpO2NadirEt = nan(Nevents,1);
ArE = nan(Nevents,1);
ArEt = nan(Nevents,1);
ArDurationE = nan(Nevents,1);
HBbaseline = nan(Nevents,1);
HBarea = nan(Nevents,1);

EventsRespBox = eval([EventsChannelName 'Box']);
EventsArBox = eval([ArChannelName 'Box']);
% pausemanually=1
for i=1:Nevents
    
    tempSig = SpO2Box(i,:);
    tempEvent = 1*(EventsRespBox(i,:)>0);
    tempEventsAr = 1*(EventsArBox(i,:)>0);
    tempTime = TimeBox(i,:);
    duration = T.EventDuration(i);
    if plotfig
    clear ax
    figure(13); clf(13);
    ax(1)=subplot(3,1,1);
    plot(TimeEnsemble,[tempEvent;tempEventsAr*0.9+0.05]);
    ylim([0 1.1]);
   
    ax(3)=subplot(3,1,3);
    end
    if 1
        I=find(Time>tempTime(1)&Time<tempTime(end));
        ChFlow=2;
        if plotfig
        plot(Time(I)-mean(Time(I)),SigT{I,ChFlow});
        end
    end
    if plotfig
    ax(2)=subplot(3,1,2);
    plot(TimeEnsemble,tempSig);
    
    linkaxes(ax,'x')
    end
    
    %get nadir SpO2
    %tempRange = 0 + delaySpO2 + [-duration*0.5 searchNadir];
    tempRange = 0 + delaySpO2 + [-duration*0.5 searchNadir];
    tempRangei = round((tempRange-TimeEnsemble(1))/dT + 1);
    tempRangei(tempRangei<1)=1;
    searchrangefull(2)=tempRangei(2);
    
    
    %code to avoid bumping into next event's baseline
    EventEndi = round((0-TimeEnsemble(1))/dT + 1);
    Inexteventstart = find(tempEvent(EventEndi+1:end)==1,1,'first');
    toli2 = toli;%round((toli + duration/2)/dT); %note using current duration as surrogate for next
    if ~isempty(Inexteventstart)
        Inexteventstart = EventEndi + Inexteventstart + delaySpO2i;
        if tempRangei(2)>(Inexteventstart + toli2)
            if verbose
                disp('run code to avoid next event');
            end
            (Inexteventstart + toli2) - tempRangei(2);
            tempRangei(2)=Inexteventstart + toli2;
        end
    end
    
    if plotfig
    hold('on');
    
    plot(TimeEnsemble(tempRangei(1):tempRangei(2)),tempSig(tempRangei(1):tempRangei(2)),'g:');
    end
    I_ = tempRangei(1):tempRangei(2);
    
    %TimeEnsemble(I_)
    [minval,~]=min(tempSig(I_));
    
    minanddiff = [[(tempSig(I_)==minval)]' [diff(tempSig(I_))' ;NaN]];
    I2_=find(minanddiff(:,1)==1 & minanddiff(:,2)>0);
    if length(I2_)==1
        if verbose
            disp('using left corner as nadir time');
        end
        
        maxi = I_(I2_(1));
    elseif length(I2_)==2
        if verbose
            disp('using left corner as nadir time');
        end
        [~,minx]=min(abs(I_(I2_) - (EventEndi+delaySpO2i)));
        maxi = I_(I2_(minx));
    else
        [~,maxi]=min(fliplr(tempSig(tempRangei(1):tempRangei(2)))); %fliplr is to get max to choose rightmost max on plateau
        maxi = tempRangei(2) - maxi + 1;
    end
    
    if plotfig
    hold on;
    plot(TimeEnsemble(maxi),tempSig(maxi),'go');
    end
    
    SpO2NadirE(i,1)=tempSig(maxi);
    SpO2NadirEt(i,1)=TimeEnsemble(maxi);
    %
    HBsearchi(2) = searchrangefull(2);
    %attempt to extend searchrangefull(2) to the right if no reduction in SpO2
    
    tempextrasearchrange = searchrangefull(2):searchrangefull(2)+round(10/dT);
    tempextrasearchrange(tempextrasearchrange>length(tempSig))=[];
    temp = tempSig(tempextrasearchrange);
    
    Ix = find(diff(temp)==-1,1);
    if isempty(Ix)
        Ix=round(10/dT);
    end
    if ~(sum(abs(diff(temp)))==0) %but don't extend if SpO2 is flat
        searchrangefull(2)=searchrangefull(2)+Ix;
    end
    
    %standard version
    [maxval,maxHBrecoveryi]=max(fliplr(tempSig(maxi:searchrangefull(2))));
        maxHBrecoveryi = searchrangefull(2) - maxHBrecoveryi + 1; %finds right most max
    if maxval >= SpO2NadirE(i,1) + 2
        if sum(abs(diff(tempSig(maxHBrecoveryi:searchrangefull(2)))))==0
            maxHBrecoveryi=searchrangefull(2);
        end
    HBsearchi(2) = maxHBrecoveryi;
    end
    
    if plotfig
    hold on;
    plot(TimeEnsemble(HBsearchi(2)),tempSig(HBsearchi(2)),'k.');
    end
    
    
%     [~,maxHBrecoveryi]=max(tempSig(maxi:maxi+searchNadiri));
%     maxHBrecoveryi = maxi + maxHBrecoveryi - 1;
%     if sum(abs(diff(tempSig(maxHBrecoveryi:maxi+searchNadiri))))==0
%         maxHBrecoveryi=maxi+searchNadiri;    
%     end
%     HBsearchi(2) = maxHBrecoveryi;
    
    
    
    %get max pre-event SpO2
    EventStarti = round((T.EventStart(i)-tempTime(1))/dT + 1);
    if EventStarti<1 %if event is actually longer than the window length
        EventStarti=1;
    end
    tempEventStartt = TimeEnsemble(EventStarti);
    
    tempRange = tempEventStartt + delaySpO2 + [searchBaseline duration*0.75];
    tempRangei = round((tempRange-TimeEnsemble(1))/dT + 1);
    if tempRangei(1)<1
        tempRangei(1)=1;
    end
    searchrangefull(1)=tempRangei(1);
    
    %code to avoid baseline after nadir:
    if tempRangei(2)>maxi
        
        if verbose
            disp('run code to baseline after nadir');
        end
        tempRangei(2)=maxi;
    end
    
    %code to avoid bumping current baseline way back into previous event's baseline
    Ilasteventend = find(tempEvent(1:EventStarti-1)==1,1,'last');
    if ~isempty(Ilasteventend)
        Ilasteventend = Ilasteventend + delaySpO2i + 1;
        tempRangei(1) - (Ilasteventend - toli);
        if tempRangei(1)<(Ilasteventend - toli)
            if verbose
                disp('run code to avoid previous baseline');
            end
            
            tempRangei(1)=Ilasteventend - toli;
        end
    end
    
    if plotfig
    hold('on')
    plot(TimeEnsemble(tempRangei(1):tempRangei(2)),tempSig(tempRangei(1):tempRangei(2)),'r');
    end
    %code to avoid setting baseline before previous event's nadir SpO2
    if i>1 && ~isnan(SpO2NadirEt(i-1))
        tprev = SpO2NadirEt(i-1) - (T.EventEnd(i) - T.EventEnd(i-1));
        tprevi = round((tprev-TimeEnsemble(1))/dT + 1);
        if tempRangei(1)<tprevi
            if verbose
                disp('run code to avoid previous baseline before last nadir');
            end
            tempRangei(1)=tprevi;
            tprev_val = tempSig(tprevi); %%%%% caused error, exceeds bounds
        end
    end
    if tempRangei(1)>tempRangei(2) %seems unlikely but has occurred
        tempRangei(1)=tempRangei(2)-1;
    end
    
    [~,maxi]=max(fliplr(tempSig(tempRangei(1):tempRangei(2)))); %fliplr is to get max to choose rightmost max on plateau
    
    maxi = tempRangei(2) - maxi + 1;
    
    %workaround to make baseline time exactly equal to last nadir time if appropriate
    if exist('tprev_val') && tempSig(tprevi) == tempSig(maxi)
        maxi = tprevi;
    end
     if plotfig
    hold on;
    plot(TimeEnsemble(maxi),tempSig(maxi),'ro');
     end
    SpO2BaselineE(i,1)=tempSig(maxi);
    SpO2BaselineEt(i,1)=TimeEnsemble(maxi);
    HBsearchi(1)=maxi;
    clear tprev_val
    
    I3_ = searchrangefull(1):searchrangefull(2);
    if any(isnan(tempSig(I3_)))
        SpO2NadirE(i,1)=NaN;
        SpO2NadirEt(i,1)=NaN;
        SpO2BaselineE(i,1)=NaN;
        SpO2BaselineEt(i,1)=NaN;
        if verbose
            disp('nans in search range')
        end
    end
    if min(diff(tempSig(I3_)))<-3
        SpO2NadirE(i,1)=NaN;
        SpO2NadirEt(i,1)=NaN;
        SpO2BaselineE(i,1)=NaN;
        SpO2BaselineEt(i,1)=NaN;
        if verbose
            disp('spo2 jumps in search range')
        end
    end
    
    %% Hypoxic burden
    HBbaseline(i) = max(tempSig(HBsearchi(1):HBsearchi(2)));
    HBarea(i) = sum(HBbaseline(i) - tempSig(HBsearchi(1):HBsearchi(2)))*dT; % %sec
    if plotfig
    hold on;    
    h=fill([TimeEnsemble(HBsearchi(1):HBsearchi(2)) fliplr(TimeEnsemble(HBsearchi(1):HBsearchi(2)))],[tempSig(HBsearchi(1):HBsearchi(2)) tempSig(HBsearchi(1):HBsearchi(2))*0+HBbaseline(i)],'r','facealpha',0.3,'Edgecolor','None');
%    pause(10);
    end
    
    %recheck for nadir
    if 0
         [~,maxi]=min(fliplr(tempSig(HBsearchi(1):HBsearchi(2)))); %fliplr is to get max to choose rightmost max on plateau
         maxi = HBsearchi(2) - maxi + 1;
         SpO2NadirE(i,1)=tempSig(maxi);
         SpO2NadirEt(i,1)=TimeEnsemble(maxi);
         
         if plotfig
             hold on;
             plot(TimeEnsemble(maxi),tempSig(maxi),'g.');
         end
         
    end
     
    %% get new onset arousal
    %     tempRange = 0 + delay(3) + [-duration/2 searchAr]; % used originally
    tempRange(1) = 0 + [-duration/2]; % to increase search time range within the evnt for Arousals
    tempRange(2) = 0 + delayAr + searchAr;
    tempRangei = round((tempRange-TimeEnsemble(1))/dT + 1);
    
    tempRangei(tempRangei<1)=1;
    %code to avoid bumping into next event's baseline
    Inexteventstart = find(tempEvent(EventEndi+1:end)==1,1,'first');
    toli3 = 0;%round((toli + duration/2)/dT); %note using current duration as surrogate for next
    if ~isempty(Inexteventstart)
        Inexteventstart = EventEndi + Inexteventstart;
        if tempRangei(2)>(Inexteventstart + toli3)
            if verbose
                disp('run code to avoid next event');
            end
            
            (Inexteventstart + toli3) - tempRangei(2);
            tempRangei(2)=Inexteventstart + toli3;
        end
    end
    ArOnsetSig = diff(tempEventsAr);
    Iar = find(ArOnsetSig(tempRangei(1):tempRangei(2))==1,1); %%%%%%error here
    if ~isempty(Iar)
        maxi = tempRangei(1) + Iar - 1;
        ArEt(i,1)=TimeEnsemble(maxi);
        Iaroff = find(ArOnsetSig(maxi:end)==-1,1);
        if isempty(Iaroff)
            Iaroff = length(ArOnsetSig);
        end
        ArDurationE(i,1) = (Iaroff-1)*dT;
        if ArDurationE(i,1)>30
            ArDurationE(i,1)=30;
        end
         if plotfig
        hold on;
        axes(ax(1));hold on;
        plot(TimeEnsemble(maxi:maxi+ArDurationE(i,1)/dT-1),[tempEventsAr((maxi:maxi+ArDurationE(i,1)/dT-1))*0.9+0.05],'m');
         end
    end
    
    ArE(i,1)=length(Iar);
    
    arinfo = [ArE(i) ArEt(i) ArDurationE(i)];
    
    if pausemanually
        pause()
    end
end


SpO2DeltaE = SpO2BaselineE - SpO2NadirE;

SpO2DelayMean = nanmean(SpO2NadirEt);
SpO2DelayMean2 = nanmean(SpO2NadirEt(SpO2DeltaE>=2));
SpO2DelayMean3 = nanmean(SpO2NadirEt(SpO2DeltaE>=3));
SpO2DelayMean4 = nanmean(SpO2NadirEt(SpO2DeltaE>=4));

SpO2BaselineMean = nanmean(SpO2BaselineE);
SpO2BaselineMedian = nanmedian(SpO2BaselineE);

%test if negative desats
NnegativeDesats = sum(SpO2DeltaE<0);
SpO2DeltaE(SpO2DeltaE<0)=0;


T.SpO2DeltaE = SpO2DeltaE;
T.ArE = ArE;

T.SpO2NadirE = SpO2NadirE;
T.SpO2NadirEt = SpO2NadirEt;
T.SpO2BaselineE = SpO2BaselineE;
T.SpO2BaselineEt = SpO2BaselineEt;
T.HBarea = HBarea;
T.HBbaseline = HBbaseline;

T.ArEt = ArEt;
T.ArDurationE = ArDurationE;

T.starttimesi=round((T.EventStart-SigT{1,1})./dt+1);
I = T.starttimesi > size(SigT,1);
if sum(I)>0 && verbose
    disp('warning: events scored after end of recording');
end
T(I,:)=[];
I = T.starttimesi < 1;
if sum(I)>0 && verbose
    disp('warning: events scored before start of recording');
end
T(I,:)=[];
            
EpochsChan = find(strcmp(ChannelsList,'Epochs')==1);
EpochsXHz=SigT.Epochs;

T.Epochs=EpochsXHz(T.starttimesi); % get state from start of event
tempi = T.starttimesi+round((T.EventDuration/2)*Fs);
    tempi(tempi>length(EpochsXHz))=length(EpochsXHz); %to stop errors if the middle of the last event is off the PSG recording
T.EpochsMidEvent=EpochsXHz(tempi); % get state from middle of event

% adding the following if loop to prevent overwriting if
% settings.supinepositioncode exists already.

if ~isfield (settings,'supinepositioncode')
    settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
    settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
    settings.supinepositioncode = settings.positioncodes(1);
end

[Position,PositionRaw] = getPos(SigT,ChannelsList,settings);

T.Position=Position(T.starttimesi);
T.PositionRaw=PositionRaw(T.starttimesi);

% if UseEventsRespAuto==1
%     Evts.RespAutoT = T;
% else
    Evts.RespT = T;
% end











