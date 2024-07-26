% This version uses LG calculation software neatened for paper publication....
function [WinInfo,BreathDataTable] = ...
    LGfromFlowBetaPart2(SigT,BreathDataTable,Vflow_out)
%% global variables
global n winNum settings ChannelsList

%% local settings
% clocktimestart = clock;
settings.plotfiguresqrtscaling=0;
settings.plotfigureLGparameters=0;
WindowDuration = settings.windowlength*60;

if settings.verbose
    clocktimestart = clock;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if settings.plotfigure==1
    LGfig = 2;
end

%%
if ~istable(BreathDataTable) % 
    disp('WARNING: data appears to be incorrectly formatted');
    % if we are here, then BreathDataTable is NOT a table, and things go
    % pear-shaped if trying dot indexing...
end

    rowofnans = nan(size(BreathDataTable,1),1);
    BreathDataTable.E1 = rowofnans; %Unable to perform assignment because dot indexing is not supported for variables of this type.
    BreathDataTable.AR3 = rowofnans;
    BreathDataTable.Ecentralapnea = rowofnans;
    BreathDataTable.Vdr_est = rowofnans;
    BreathDataTable.VAr_est = rowofnans;
    BreathDataTable.Error = rowofnans;


%%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Control parameters for operation of algorithm
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ArousalShift=settings.ARmodel;

if ArousalShift==0 %value of 1 should lead to normal use of VRA; was inverted from ?? until 10/7/19 -SS
    % was correct in 20190205 code version, was incorrect in 20190404 code version and 20190629 until 10/7/19 as above. (DLM)
    VraOn=0;
else
    VraOn=1;
end

% mark_onset_of_arousal_only=0;
delete_arousal_in_events=1;
seteventsbasedonvlessthaneupnea2=7;
eventslessthanX=0.7; %default<0.7
findcentralhypopneasandapneas=settings.findcentralhypopneasandapneas;
if settings.havescoredcentralhypops
    removecentralhypops=0;
else
    removecentralhypops=1;
end

numArousalBreaths=2;
dt=1/settings.Fs;%DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1);
% filter_LFcutoff_butter1 = 0.03;
% filter_HFcutoff_butter1 = 10;
% filter_order1 = 2;
% [B_butter1,A_butter1] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Extraction and Preprocessing of flow data for calculation of LG:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%Flow=((-1)^settings.Pnasaldownisinsp)*((-1)^invertflow)*DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %%% added invert flow

PosChan = find(strcmp(ChannelsList,'Position')==1);
if ~isempty(PosChan)
    PositionSig = 0.1*round(10*SigT.Position);
    position_mode=mode(PositionSig); %%% added round position
    N_position_changes=sum(abs(diff(SigT.Position)));
    Percent_position=mean(position_mode==SigT.Position); %does this work(?)
else
    PositionSig = NaN*SigT.Time;
    position_mode = NaN;
    N_position_changes = NaN;
    Percent_position = NaN;
end

Time=SigT.Time;

if isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==1
    Arousal=SigT.EventsArWS; %EEG selection uses scoring
    ArousalOrig=SigT.EventsAr; %For plots
elseif isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==2
    Arousal=SigT.EventsArWSB; %EEG selection does not use scoring
    ArousalOrig=SigT.EventsAr; %For plots
else
    Arousal=SigT.EventsAr; %Use original arousal scoring
end

hypnog=SigT.Epochs;

spo2Chan=find(strcmp(ChannelsList,'SpO2')==1);
if ~isempty(spo2Chan)
    try
        spo2=SigT.SpO2;
        filter_HFcutoff_butter0 = 2; filter_order0 = 2;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
        spo2 = filtfilt(B_butter0,A_butter0,spo2);
        spo2 = round(spo2);
        spo2(spo2>100)=100;
        spo2(spo2<45)=NaN;
    catch
        spo2 = NaN*SigT.Time;
    end
else
    spo2 = NaN*SigT.Time;
end

%% Pull frequently used data from Table
BB_i_start = BreathDataTable.BB_i_start;
BB_i_mid = BreathDataTable.BB_i_mid;
BB_i_end = BreathDataTable.BB_i_end;
VI = BreathDataTable.VI;
AR = BreathDataTable.AR;
BB_t = Time(BB_i_start);
Ttot_B = BreathDataTable.Time_end - BreathDataTable.Time_start;
meanVIbeforenormalizing = BreathDataTable.Veup;
DeltaPes = BreathDataTable.DeltaPes;
DeltaPmus = BreathDataTable.DeltaPmus;
VIpes = BreathDataTable.VIpes;
DeltaEdi = BreathDataTable.DeltaEdi;
PeakEdi = BreathDataTable.PeakEdi;
VIedi = BreathDataTable.VIedi;
GGpeak = BreathDataTable.GGpeak;
GGtonic = BreathDataTable.GGtonic;
FlowPes_VI = BreathDataTable.FlowPes_VI;
FlowEdi_VI = BreathDataTable.FlowEdi_VI;


%% Load in the scored events.
%obstructive if obstrucive apnea / hypopnea / or mixed apnea

EventsChan = find(strcmp(ChannelsList,'EventsResp')==1);

if isfield(settings,'AutoScoreRespEvents') && settings.AutoScoreRespEvents==1 && settings.UseAutoScoredRespEventsForLG==1
    EventsChan = find(strcmp(ChannelsList,'EventsRespAuto')==1);
    
    Etype = 0*BB_i_start;
    for k=1:length(BB_i_start)
        % 3=CentralA 2=ObstructiveA 4=OHypopnea 5=MixedA 6=CentralHypopnea
        temp = SigT{BB_i_start(k):BB_i_end(k),EventsChan};
        F = sum(temp~=0)/length(temp);
        if settings.eventsarebreathsfullywithinmargins==1
            if F==1
                Etype(k)=mode(temp(temp>0));
            end
        else
            if F>0
                Etype(k)=mode(temp(temp>0));
            end
        end
    end
    BreathDataTable.Etype=Etype;
    
end

Events=sum(SigT{:,EventsChan}'==[2 4 5]')';
EventsScoredCentralApnea=1*[SigT{:,EventsChan}==3];
EventsScoredCentralHyp=1*[SigT{:,EventsChan}==6];

%first half of mixed event is central apnea
if settings.handlemixedeventsseparately
    try
        EventsMixed = 1*[SigT{:,EventsChan}==5];
        if sum(EventsMixed)>0
            Imixed = find([NaN;diff(EventsMixed)]==1); %index of arousal onset
            Imixed2 = find(diff(EventsMixed)==-1); %index of pre arousal onset
            if Imixed2(1)<Imixed(1) %tidy up in case they don't match
                Imixed = [1;Imixed];
            end
            if Imixed2(end)<Imixed(end)
                Imixed2 = [Imixed2;length(EventsMixed)];
            end
            lengthmixedE=Imixed2-Imixed+1;
            for i=1:length(Imixed)
                ImixedMID = round(lengthmixedE(i)/2)+Imixed(i);
                EventsScoredCentralApnea(Imixed(i):ImixedMID)=1;
                Events(Imixed(i):ImixedMID)=0;
            end
        end
    catch me
    end
end

%% Convert events to breath domain
E=zeros(1,length(BB_i_start)); %Events
E_ScoredCentralApnea=zeros(1,length(BB_i_start)); %EventsScoredCentralApnea
E_ScoredCentralHyp=zeros(1,length(BB_i_start)); %EventsScoredCentralHyp
Te_i1 = round(median(BB_i_end-BB_i_start));
for k=1:length(BB_i_start)
        if k==1
            I=max([BB_i_start(k)-Te_i1 1]):BB_i_mid(k);
        else
            I=BB_i_mid(k-1):BB_i_mid(k);
        end
        I2=BB_i_start(k):BB_i_end(k);
        E(k)=1-(max(Events(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        if settings.eventsarebreathsfullywithinmargins
            E(k)=1-(min(Events(BB_i_start(k):BB_i_end(k)))); %if breath is entirely inside an event
            E_ScoredCentralApnea(k)=(min(EventsScoredCentralApnea(BB_i_start(k):BB_i_end(k)))); %if event is within breath
            E_ScoredCentralHyp(k)=(min(EventsScoredCentralHyp(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        else
            E(k)=1-(max(Events(BB_i_start(k):BB_i_end(k)))); %if breath is even slightly inside an event
            E_ScoredCentralApnea(k)=(max(EventsScoredCentralApnea(BB_i_start(k):BB_i_end(k)))); %if event is within breath
            E_ScoredCentralHyp(k)=(max(EventsScoredCentralHyp(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        end
end

%% RIP analysis
% included as (1) placeholder for future development, (2) interesting, (3) cool plot
% assign variables here, they are used below, even if not here for TAA
try
    Ab = SigT.Abdomen;
    Th = SigT.Thorax;
catch
    Ab=[]; Th=[];
end

%if there are no Th or Ab data the channel will be set to zeros
if isempty(Th)
    Th = 0*Time;
end
if isempty(Ab)
    Ab = 0*Time;
end

%% Pes analysis (if it exists)
if sum(strcmp(ChannelsList,'Pes')==1)
    Pes = SigT.Pes; % not indexed, 'find' stays
    %Filtered Pes (gentle low pass)
    [B_butter0,A_butter0] = butter(2,5/(1/dt/2),'low');
    Pes = filtfilt(B_butter0,A_butter0,Pes);
end

%% Edi analysis (if it exists)
if sum(strcmp(ChannelsList,'Edi')==1)
    Edi = SigT.Edi;
end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calculate LG from ventilation data:

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
TimeB = BB_t-BB_t(1); % The time of each breath, from the first breath....[breath time series]

% *********************************************************************
% if there is a missing section of samples either at the end or the
% beginning of the period, break.
% *********************************************************************
if TimeB(end)<0.8*WindowDuration
    WinInfo(1:22)=NaN;
    disp('exiting 4: ');
    return
end

% *********************************************************************

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Update events based on eupnea, and detection of "detected" central
% apneas... i.e. any breath greater than eupnoea has events cleared,
% any breath at less than 70% of eupnoea is set as an event. Also, any
% breath detected as a central event (depending on criteria) has event
% cleared for purpose of model fitting.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    plot_on2=settings.plotfigure;
else
    plot_on2=0; %Override - don't plot plots from ObstructBreathScore_V14_0
end

[E1,Ecentralapnea,Ecentralhypop] = ObstructBreathScore_V14_1(1,E,VI,BB_i_start,removecentralhypops,eventslessthanX,plot_on2);
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if settings.usescoredcentralapneas %Use scored central events (if they exist)
    E1(E_ScoredCentralHyp|E_ScoredCentralApnea)=1; %1=unobstructed
    Ecentralapnea(E_ScoredCentralApnea==1)=1; %1=central apnea
    Ecentralhypop(E_ScoredCentralHyp==1)=1; %1=central hypopnea
end

if isfield(settings,'lowVEareObstructiveEvents') && settings.lowVEareObstructiveEvents ~= 0
    %try additional events manipulation based on flow trace:
    %reset:
    E1(VI<=settings.lowVEareObstructiveEvents)=0;
    %E1(VI>settings.lowVEareObstructiveEvents)=1;
    E1(Ecentralapnea==1)=1;
    E1(Ecentralhypop==1)=1;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Alternative events, purely for the purpose of arousal marking
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
E_recover=zeros(1,length(VI));
E_Terminate=zeros(1,length(VI));
E_recover(VI>1)=1;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    if settings.plotfigure
        figure(100),
        ax100(1)=subplot(311);plot(SigT{BB_i_start,1},VI,SigT{BB_i_start,1},[zeros(1,length(BB_i_start));ones(1,length(BB_i_start))],':k');
        ax100(2)=subplot(312);plot(SigT{BB_i_start,1},VEThorax,SigT{BB_i_start,1},VEAbdo,SigT{BB_i_start,1},[zeros(1,length(BB_i_start));ones(1,length(BB_i_start))],':k');
        ax100(3)=subplot(313);plot(SigT{BB_i_start,1},Ecentralapnea,'.',SigT{BB_i_start,1},Ecentralhypop,'.');
        linkaxes(ax100,'x');
    end
end


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This section identifies the duration of each arousal and links it to
% the breath domain arousal markings
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR_Rise=AR(2:end)-AR(1:end-1);
AR_Rise(AR_Rise<1)=0;
ArousalDuration=zeros(1,length(AR));
for jk=1:1:length(AR_Rise)
    if AR_Rise(jk)==1
        time_start=BB_i_start(jk+1);
        time_count=time_start;
        while Arousal(time_count)==0&&time_count<length(Time)
            time_count=time_count+1;
        end
        while Arousal(time_count)==1&&time_count<length(Time)
            time_count=time_count+1;
        end
        ArousalDuration(jk+1)=Time(time_count)-Time(time_start);
        
        %AR(k)=max(Arousal(BB_i_start(k):BB_i_mid(k)));
    end
end

ArousalCount=1;
for jk=2:1:length(Arousal)
    if Arousal(jk)==1&&Arousal(jk-1)==0 % Arousal Start
        ArousalDat(ArousalCount,1)=Time(jk);
    elseif Arousal(jk)==0&&Arousal(jk-1)==1 %Arousal End
        ArousalDat(ArousalCount,2)=Time(jk);
        ArousalCount=ArousalCount+1;
    end
end
if ArousalCount==1
    ArousalDat(ArousalCount,2)=NaN;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%This section places notations on event markings to allow
%implementation of the desired VRA modelling and data fitting.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%%
AR1=AR;

if delete_arousal_in_events
 %   [ AR1, numArousalBreaths, E_Terminate] = ArousalShift_V14_0b(ArousalShift,AR1,E1,E_Terminate,BB_i_start,numArousalBreaths);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Run the autoregressive loop gain calculation for a range of delays,
% and pick the one that results in the best (lowest) cost function.
%**********************************************************************
%**********************************************************************
forcedelays=[1:settings.maxdelaybreaths NaN]; % delay from 1 to 5 breaths.
Nforcedelays=length(forcedelays)-1; % number of iterations of delay before searching for the lowest resultant cost function.


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Identify the best starting breath from the window, based on: starting
% after the maximum possible delay, and being the start of the longest
% section of data in the first 2 minutes that is costed in the AR
% function (i.e: no event)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Ttot_ser=(diff(TimeB(:))); % T-tot series
Ttot_med=median(Ttot_ser);
Ttot_ser_prev=[Ttot_med;Ttot_ser(1:end)];

%% Arousal decay model
if delete_arousal_in_events
    AR1b=AR1;
    AR1b(E1==0|E_ScoredCentralApnea==1|E_ScoredCentralHyp==1)=0;
    Ndeleted=sum(AR1)-sum(AR1b);
    AR1=AR1b;
    clear AR1b Ndeleted;
end

arousaldecaymethod=0;
if arousaldecaymethod
    AR3=AR1';
    try
        if sum(AR3)>0
            
            Iar = find([NaN;diff(AR3)]==1); %index of arousal onset
            Iar2 = find(diff(AR3)==-1); %index of pre arousal onset
            if Iar2(1)<Iar(1) %tidy up in case they don't match
                Iar = [1;Iar];
            end
            if Iar2(end)<Iar(end)
                Iar2 = [Iar2;length(AR3)];
            end
            lengtharousalB=Iar2-Iar+1;
            for i=1:length(Iar)
                if 1 %create signal for length of arousal
                    data=2.^(-1*(0:lengtharousalB(i)-1)); %[1 0.5 0.25 etc]
                else
                    data=(lengtharousalB(i):-1:1)/lengtharousalB(i); %linear reduction with no effect at offset
                end
                AR3(Iar(i):(Iar(i)+lengtharousalB(i)-1))=data; %save signal as appropriate in AR3
            end
        end
        AR3=AR3'*1;
    catch me
    end
else
    AR3=AR1;
end

%still seems to be arousals in events so writing code to make it actually happen here
%% Calculate Veupnea
Veupnea=1;

%% MAIN PLOT
while 1
    Ventilation_EventScoring_data=[VI(:)-Veupnea E1(:) AR3(1:length(E1)) TimeB(:) Ecentralapnea(:) [Ttot_ser(:); 0] Ttot_ser_prev(:)];
    
    %
    
    
    %find arousal onset
    ARonset=AR1;
    ARonset(1:end)=0;
    for i=2:1:length(AR1)
        if AR1(i)==1&&AR1(i-1)==0
            ARonset(i)=1;
        end
    end
    
    %counting events, note this only counts obstructive events.
    N_events=0;
    for k=2:length(E)
        if E(k)==0&&E(k-1)==1
            N_events=N_events+1;
        end
    end
    
    
    modelpolyfitorder=3;
    temp_i_i=10;
    if sum(AR3((temp_i_i+1):end)==0)==0
        VraOn=0;
    end
    [Error,Vdr_est,VAr_est,LoopGainParameters, Optimal_vent_control_parameters,Optimal_MSE,FitQual,i_1,i_end,...
        lowerSEM,upperSEM,CI_parameters] = ...
        FindBestModelParametersCI(Ventilation_EventScoring_data,VraOn,Veupnea,modelpolyfitorder);
    %[Error,Vdr_est,VAr_est,LoopGainParameters,Optimal_vent_control_parameters,Optimal_MSE,FitQual,i_1,i_end] = ...
    %   FindBestModelParameters(Ventilation_EventScoring_data,VraOn,Veupnea,modelpolyfitorder);
    
    if VraOn==0
        LoopGainParameters(8)=NaN; %VRA
    end
    
    %Vdr_est=Vdr_est_X
    
    if isnan(LoopGainParameters(1))==1
        WinInfo(1:22)=NaN;
        disp('exiting 6');
        return
    end
    
    Arthressample=Vdr_est+Veupnea;
    Arthressample([1:i_1 i_end:end])=[];
    ARonset([1:i_1 i_end:end])=[];
    Arthressample(ARonset==0)=[];
    Arthres=mean(Arthressample);
    N_arousals_used=sum(ARonset);
    
    % calculate this outside of plot if statement
    meanSEM = nanmean(upperSEM(i_1:end)-lowerSEM(i_1:end));
    
    %% 
    RIPinfo='';
    
    %% Plot
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %If plotfigure=1, create plot
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if settings.plotfigure==1
        figure(LGfig); clf(LGfig); LGfigID = gcf;
        LGfigID.Color = [1 1 1];
        LGfigID.Units = 'Inches';
        %LGfigID.Position = [8 3 11 8]; %removed, irritating :)
        
        if 1
            plotflow=Vflow_out; plotflow_t=Time;
        else
            plotflow=Vflow; plotflow_t=time;
        end
        if 1
            RIP_Thorax_plot=Th;
            RIP_Abdo_plot=Ab;
            plotfilteredRIP=0;
            if plotfilteredRIP
                RIP_Thorax_plot = filtfilt(B_butter1,A_butter1,RIP_Thorax_plot);
                RIP_Abdo_plot = filtfilt(B_butter1,A_butter1,RIP_Abdo_plot);
            end
        end
        set(gcf,'Color',[1 1 1]);
        ax1(1)=subplot(3,1,1);
        
        plot(Time,hypnog,'k');
        
        %plot(BB_t(AR1>0),AR1(AR1>0),'r.',BB_t(E1<1),E1(E1<1),'k.',BB_t(E<1),E(E<1)+1,'c.',BB_t(Ecentralhypop>0),Ecentralhypop(Ecentralhypop>0)-2,'b.',...
        %    BB_t(Ecentralapnea>0),Ecentralapnea(Ecentralapnea>0)-3,'g.');
        set(gca,'XColor',[1 1 1],'YColor',[0 0 0],'ylim',[-0.1 5.1],...
            'FontSize',7,...
            'Position',[0.131 0.782 0.775 0.071],...
            'YTickLabel',{'3','2','1','R','W'},...
            'YTick',[0:4],...
            'TickDir','out',...
            'TickLength',[0 0],...
            'box','off',...
            'FontName','Arial');
        %line([Time(1) Time(1)],[-0.1 5.1]);
        xlim([min(Time) max(Time)]);
        title(['Patient ' num2str(n) ', Window ' num2str(winNum) ' ' RIPinfo],'fontname','Arial Narrow','fontsize',12);
        
        ax1(2)=subplot(3,1,2);
        
        box('off');
        xlim([min(Time) max(Time)]);
        hold('on')
        
        alpha_=98.5; % Mesa version set this to 99.5
        nflow=1.5*(prctile(plotflow,alpha_)-prctile(plotflow,100-alpha_));    mflow=mean(plotflow);
        nTH=2*(prctile(RIP_Thorax_plot,alpha_)-prctile(RIP_Thorax_plot,100-alpha_));    mTH=mean(RIP_Thorax_plot);
        nAB=2*(prctile(RIP_Abdo_plot,alpha_)-prctile(RIP_Abdo_plot,100-alpha_));     mAB=mean(RIP_Abdo_plot);
        minspo2=min(spo2);
        if minspo2>90
            minspo2=90;
        end
        nspo2=1.5*(100-minspo2);
        mspo2=mean([100,minspo2]);
        
        %plot EEG and arousals scored over the top in green
        if max(Arousal)>0
            tempds=4; tempsize=0.2; tempoffset=1.2-0.0;
            tempx=downsample(Time,tempds);
            tempy=tempsize*downsample(Arousal,tempds)+tempoffset;
            tempxx=[tempx;flipud(tempx)];
            tempyy=[tempy;tempoffset+zeros(length(tempy),1)];
            fill(tempxx,tempyy,[0 1 0],'EdgeColor','none','FaceAlpha',0.75);
            %fill(tempxx,tempyy,[0 1 0],'EdgeColor','none');
        end
        if exist('ArousalOrig') && max(ArousalOrig)>0 && 1 %compare scoring, plot original scoring higherup
            tempds=4; tempsize=0.2; tempoffset=1.2 - 0.0 + 0.1;
            tempx=downsample(Time,tempds);
            tempy=tempsize*downsample(ArousalOrig,tempds)+tempoffset;
            tempxx=[tempx;flipud(tempx)];
            tempyy=[tempy;tempoffset+zeros(length(tempy),1)];
            fill(tempxx,tempyy,[0.05 0.2 0.9],'EdgeColor','none','FaceAlpha',0.75);
            %fill(tempxx,tempyy,[0 1 0],'EdgeColor','none');
        end
        if sum(strcmp(ChannelsList,'EEG')==1)&&1
            EEGplot = SigT.EEG;
            plot(Time,1.1+0.1*EEGplot/(prctile(EEGplot,99.5)-prctile(EEGplot,0.05)),'color',[0.1 0.5 0.3],'linewidth',0.25);
        end
        
        if sum(strcmp(ChannelsList,'WakeSleep')==1)
            WakeSleep = SigT.WakeSleep;
            plot(Time,1.2+0.2*WakeSleep,'k');
        end
        if sum(strcmp(ChannelsList,'WPr')==1)
            WPr = SigT.WPr;
            ArPr = SigT.ArPr;
            WAPr = max([WPr ArPr]')';
            plot(Time,1.2+0.2*WAPr,'r');
        end
        if sum(strcmp(ChannelsList,'WPrB')==1)
            WPr = SigT.WPrB;
            ArPr = SigT.ArPrB;
            WAPr = max([WPr ArPr]')';
            plot(Time,1.2+0.2*WAPr,'r');
        end
        
        plotbreathgrid=1;
        flowplotpos=0.5;
        if plotbreathgrid
            gridheight=0.5;
            xgridlist = Time(BB_i_start)';
            plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^2);
            xgridlist = Time(BB_i_mid)';
            plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^4);
            %zero flow baseline:
            plot(plotflow_t,flowplotpos+0*plotflow_t,'color',0.8*[1 0.5 0.5]);
        end
        
        %back to plotting flow, RIP, SpO2 (over EEG)
        plot(plotflow_t,0.5+(plotflow-mflow)/nflow,'k',...
            Time,-0.25+(RIP_Thorax_plot-mTH)/nTH,'k',...)
            Time,-1+(RIP_Abdo_plot-mAB)/nAB,'k',...
            Time,-1.75+(spo2-mspo2)/nspo2,'r');
        plotpos = -1.75;
        
        set(gca,'YTick',[-1.75 -1 -0.25 0.5 1.2],...
            'XColor',[1 1 1],...
            'TickDir','out',...
            'TickLength',[0 0],...
            'Position',[0.131 0.455 0.775 0.321],'FontName','Arial','FontSize',7);
        box('off');
        ylim([-2.35 1.6]);
        set(gca,'YDir','normal','FontSize',7);
        set(gca,'YTickLabel',{'SpO2', 'Abdomen', 'Thorax','Flow','EEG'});
        
        optionalchannelslist = {'Edi','Pes','GGpmax'};
        if 1 %clip Pes upper
            try
                Pes(Pes>Pesbaselineest+20)=NaN; %overwrite just for plot, careful not to analyze more below
            catch me
            end
        end
        for aa=1:length(optionalchannelslist)
            try
            channeltoplot = optionalchannelslist{aa};
            if sum(strcmp(ChannelsList,channeltoplot)==1)
                hold('on')
                tempupper=prctile(eval(channeltoplot),alpha_);
                templower=prctile(eval(channeltoplot),100-alpha_);
                tempdelta=tempupper-templower;
                tempupper=tempupper-tempdelta*0.1;
                templower=templower+tempdelta*0.1;
                tempdelta=tempupper-templower;
                tempmid = templower+tempdelta*0.5;
                %nX=1.5*(prctile(eval(channeltoplot),alpha_)-prctile(eval(channeltoplot),100-alpha_)); mX=nanmean(eval(channeltoplot));
                plotpos=plotpos-0.75;
                plot(plotflow_t,plotpos+(eval(channeltoplot)-tempmid)/tempdelta/2.2,'k')
                hold('off')
                temp = get(gca,'YTick');
                set(gca,'YTick',[temp(1)-0.75 temp]);
                temp = get(gca,'Ylim');
                set(gca,'Ylim',[temp(1)-0.75 temp(2)]);
                temp = get(gca,'YTickLabel');
                set(gca,'YTickLabel',[channeltoplot;temp]);
            end
            catch
            end
        end
        
        if 1
            %plot final event breaths over the top of flow in blue
            hold('on')
            dt_new=0.25;
            startT=ceil(BB_t(1));
            endT=ceil(BB_t(length(BB_t)));
            time_dt=startT:dt_new:endT;
            E1_rs=0*time_dt;
            for i=1:length(time_dt)
                E1_rs(i) = 1-E1(find(Time(BB_i_start)<=time_dt(i),1,'last'));
            end
            C1_rs=0*time_dt;
            for i=1:length(time_dt)
                C1_rs(i) = Ecentralapnea(find(Time(BB_i_start)<=time_dt(i),1,'last'));
            end
            E1_rs=E1_rs&(~C1_rs);
            if sum(E1_rs)>0
                tempsize=0.25; tempoffset=0.375;
                tempy=tempsize*(E1_rs)+tempoffset;
                tempxx=[time_dt fliplr(time_dt)];
                tempyy=[tempy tempoffset+zeros(1,length(tempy))];
                fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
            end
            hold('on')
            if sum(C1_rs)>0
                tempsize=0.25; tempoffset=0.375;
                tempy=tempsize*(C1_rs)+tempoffset;
                tempxx=[time_dt fliplr(time_dt)];
                tempyy=[tempy tempoffset+zeros(1,length(tempy))];
                fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
            end
            %figure(980); plot(time_dt,E1_rs,'r',BB_t,E1,'r.')
        end
        
        ax1(3)=subplot(313);
        hold('on');
        %plot event breaths over the top of flow in blue
        tempsize=0.5; tempoffset=0.75;
        tempy=tempsize*(E1_rs)+tempoffset;
        tempxx=[time_dt fliplr(time_dt)];
        tempyy=[tempy tempoffset+zeros(1,length(tempy))];
        fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
        tempsize=0.5; tempoffset=0.75; %updated 6/2014
        tempy=tempsize*(C1_rs)+tempoffset;
        tempxx=[time_dt fliplr(time_dt)];
        tempyy=[tempy tempoffset+zeros(1,length(tempy))];
        fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
        
        if 0
            try
                realVdrive = FlowEdi_VI/(meanVIbeforenormalizing);
                stairs(BB_t(1:i_end),realVdrive(1:i_end),'color',[1 0.6 0.6]);
                realVdrive2 = FlowPes_VI/(meanVIbeforenormalizing);
                stairs(BB_t(1:i_end),realVdrive2(1:i_end),'color',[0.6 1 1]);
            catch me
            end
        end
        
        stairs(BB_t(1:i_end),Ventilation_EventScoring_data(1:i_end,1)+Veupnea,'color',[0.4 0.4 0.7]);
        plot(BB_t(1:i_end),[zeros(1,length(BB_t(1:i_end))); ones(1,length(BB_t(1:i_end)))],'k:');
        stairs(BB_t(i_1:i_end)+0.3,Vdr_est(i_1:end)'+Veupnea+VAr_est(i_1:i_end)','color',[0.1 0.8 0.1]);
        plot(BB_t(i_1:i_end),Vdr_est(i_1:end)'+Veupnea,'color',[0 0 0]);
        
        % meanSEM = nanmean(upperSEM(i_1:end)-lowerSEM(i_1:end));
        
        %         if exist('lowerSEM')&&meanSEM<1&&0
        %         plot(BB_t(i_1:i_end),lowerSEM(i_1:end)'+Veupnea,'color',[0.6 0.6 0.6]);
        %         plot(BB_t(i_1:i_end),upperSEM(i_1:end)'+Veupnea,'color',[0.6 0.6 0.6]);
        %         end
        %plot(BB_t(i_1:i_end),W(i_1:end)','r');
        %plot(BB_t(i_1:i_end),VAr_est(i_1:i_end),'b');
        %xlabel(num2str(LGplusinfo(4),2))
        %LGplusinfo=[Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA_est];
        if ~isnan(Arthres)
            arthrestext = [' ArThr=' num2str(Arthres,2)];
        else
            arthrestext=[];
        end
        %xlabel([ 'LG1=' num2str(LoopGainParameters(6),2) ' Tn=' num2str(LoopGainParameters(5),2) ' LGn=' num2str(LoopGainParameters(4),2) ' delay=' num2str(LoopGainParameters(3),2) arthrestext ' RsqTotal=' num2str(FitQual(2),2)])
        xlabel([ 'Pos=', num2str(position_mode,1), '(', num2str(Percent_position,0), ')',...
            ' meanSEM=', num2str(meanSEM,2),...
            ' LG1=', num2str(LoopGainParameters(6),2),'(', num2str(CI_parameters(1),2), ')',...
            ' Tn=', num2str(LoopGainParameters(5),2), ' LGn=', num2str(LoopGainParameters(4),2), ...
            ' delay=', num2str(LoopGainParameters(3),2), '(', num2str(CI_parameters(5),2), ')',...
            arthrestext,' RsqTotal=', num2str(FitQual(2),2)]);
        
        hold('off');
        set(gca,'FontSize',7,...
            'Position',[0.131 0.221 0.775 0.216],...
            'TickLength',[0.0 0.0],...
            'TickDir','out',...
            'xtick',[Time(1) max(Time)],'xticklabel',{datestr(Time(1)/86400,'HH:MM:SS'),datestr(max(Time)/86400,'HH:MM:SS')},...
            'FontName','Arial');
        box('off');
        
        xlim([min(Time) max(Time)]);
        
        linkaxes(ax1,'x');
        
        if settings.manualscoringtouchups
            disp('click left of data if no touch ups needed, or left then right to add events in range, or right then left to remove events in range');
            [t(1),~,~]=ginput(1);
            if t(1)<min(Time)
                manualscoringtouchupsincomplete=0;
            else
                [t(2),~,button]=ginput(1);
                manualscoringtouchupsincomplete=1;
                if button==1
                    if t(1)<t(2)
                        E1(BB_t>t(1)&BB_t<t(2))=0; %add events...
                    elseif t(1)>t(2)
                        E1(BB_t>t(2)&BB_t<t(1))=1; %remove events...
                    end
                elseif button>1 %right click
                    if t(1)<t(2)
                        Ecentralapnea(BB_t>t(1)&BB_t<t(2))=1; %add events...
                    elseif t(1)>t(2)
                        Ecentralapnea(BB_t>t(2)&BB_t<t(1))=0; %remove events...
                    end
                end
            end
        else
            manualscoringtouchupsincomplete=0;
        end
    else
        manualscoringtouchupsincomplete=0;
    end
    
    %%
    if manualscoringtouchupsincomplete==0
        break
    end
    
end %end while manualscoringtouchups loop

%% save loop gain plot
if settings.saveplots
    % old
    if 0
        savefig(gcf,[settings.OutputDataDirectory '\' settings.savename '_n=' num2str(n) '_w=' num2str(winNum)],'compact');
    else 
        % save in folder for each patient
        saveloc=[settings.OutputDataDirectory settings.savename filesep settings.filename];
        if ~(exist(saveloc, 'dir') == 7)
            mkdir(saveloc);
        end
        % save LGfigID in the file format specified.
        switch settings.savefigas
            case 'saveasTIFF'
                print(LGfigID, [saveloc, '\LG_window=',num2str(winNum)], '-dtiff', '-r300');               
            case 'saveasPNG'
                saveas(LGfigID, [saveloc, '\LG_window=',num2str(winNum)], 'png');                
            case 'saveasFIG' 
                savefig(LGfigID, [saveloc, '\LG_window=',num2str(winNum)],'compact');                
        end
    end
end

%%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%LGplusinfo=[Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA_est];
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% LG output vector which includes some other key results parameters.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
LoopGainParameters=[BB_t(i_1) LoopGainParameters(1:2) 0 LoopGainParameters(4:7) LoopGainParameters(3) LoopGainParameters(8) 0 Arthres Optimal_MSE mean(Ttot_B) std(Ttot_B) median(Ttot_B) iqr(Ttot_B)];
LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position CI_parameters' meanSEM];

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%'LGplusinfo=[breath Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA1_est VRA2 Arthres1 Optimal_MSE mean(Ttot_B) std(Ttot_B) median(Ttot_B) iqr(Ttot_B)];'
% old
%'LGplusinfo=[breath Gtot_est tau1_est delay_est LG_n_est T_n_est LG_1_est LG_2_est VRA_est Arthres Optimal_MSE mean(Ttot_B) std(Ttot_B) median(Ttot_B) iqr(Ttot_B)]'
% Current.
%LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position CI_parameters' meanSEM];


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calculate the Events Information Parameters:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
E=Ventilation_EventScoring_data(:,2);
AR=Ventilation_EventScoring_data(:,3);
eventHist=zeros(length(E),1);
eventCurrent=0;
eventCount=0;
for b=1:1:length(E)
    if eventCurrent==0
        if E(b)==0
            eventCurrent=1;
            eventCount=eventCount+1;
        elseif E(b)==1
            
        end
    elseif eventCurrent==1
        if E(b)==0
            eventCount=eventCount+1;
        elseif E(b)==1
            eventHist(eventCount)=eventHist(eventCount)+1;
            eventCurrent=0;
            eventCount=0;
        end
    end
    if (b==length(E))&&(eventCurrent==1)
        eventHist(eventCount)=eventHist(eventCount)+1;
    end
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Breath by Breath Events:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR_Rise=AR(i_1+1:end)-AR(i_1:end-1);
AR_Rise(AR_Rise<1)=0;
Eup=1-E;
EventRise=Eup(i_1+1:end)-Eup(i_1:end-1);
EventRise(EventRise<1)=0;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Time Based Events:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR_Time=Arousal(ceil(TimeB(i_1)*dt)+1:ceil(TimeB(i_1)*dt+600*dt))-Arousal(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt-1));
AR_Time(AR_Time<1)=0;
EventTime=Events(ceil(TimeB(i_1)*dt)+1:ceil(TimeB(i_1)*dt+600*dt))-Events(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt-1));
EventTime(EventTime<1)=0;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%EventsInfo can be removed:
EventsInfo=[length(Events)/100 mean(Arousal(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt))) sum(AR_Time)  mean(Events(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt)))  sum(EventTime) length(VI) mean(AR(i_1+1:end)) sum(AR_Rise)  mean(Eup)  sum(EventRise) ];


%
EventHistOut{1}=eventHist;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
%     size(BB_i_start(i_1:i_end))
%     size(BB_i_mid(i_1:i_end))
%     size(BB_i_end(i_1:i_end))
%     size(E1(i_1:i_end))
%     size(E_recover(i_1:i_end))
%     size(W(i_1:i_end))
%     size(Vdr_est(i_1:i_end))
%
%DataOut=[Time(BB_i_start(i_1:i_end)) Time(BB_i_mid(i_1:i_end)) Time(BB_i_end(i_1:i_end)) VI(i_1:i_end).' Vdr_est(i_1:i_end).' E1(i_1:i_end).' E_recover(i_1:i_end).' E_Terminate(i_1:i_end).' W(i_1:i_end).' ];
pos_B = PositionSig(BB_i_start);
pos_B = pos_B(:);
hypnog_B = SigT{BB_i_start(:),find(strcmp(ChannelsList,'Epochs')==1)};
hypnog_B = hypnog_B(:);

% DataOut=[Time(BB_i_start(i_1:i_end)) Time(BB_i_mid(i_1:i_end)) Time(BB_i_end(i_1:i_end)) ...
%         VI(i_1:i_end).' Vdr_est(i_1:i_end) E1(i_1:i_end).' E_recover(i_1:i_end).' E_Terminate(i_1:i_end).' Error(i_1:i_end) AR(i_1:i_end) ...
%         meanVIbeforenormalizing+0*AR(i_1:i_end) VAr_est(i_1:i_end) pos_B hypnog_B ...
%         meanVIbeforenormalizing*VI(i_1:i_end).' DeltaPes(i_1:i_end) DeltaPmus(i_1:i_end) VIpes(i_1:i_end) ...
%         DeltaEdi(i_1:i_end) VIedi(i_1:i_end) GGpeak(i_1:i_end) GGtonic(i_1:i_end) FlowPes_VI(i_1:i_end) FlowEdi_VI(i_1:i_end)];
DataOut=[Time(BB_i_start) Time(BB_i_mid) Time(BB_i_end) ...
    VI(:) Vdr_est(:) E1(:) E_recover(:) E_Terminate(:) Error(:) AR(:) ...
    meanVIbeforenormalizing(:) VAr_est(:) pos_B(:) hypnog_B(:) ...
    meanVIbeforenormalizing.*VI(:) DeltaPes(:) DeltaPmus(:) VIpes(:) ...
    DeltaEdi(:) VIedi(:) GGpeak(:) GGtonic(:) FlowPes_VI(:) FlowEdi_VI(:) PeakEdi(:)];

if settings.verbose
    delta_t = etime(clock, clocktimestart); % delta in seconds
    D_win = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
    str = ['Analysis of window complete: ', char(D_win), ' (hh:mm:ss)']; disp(str);
end

    BreathDataTable.E1 = E1(:);
    BreathDataTable.AR3 = AR3(:);
    BreathDataTable.Ecentralapnea = Ecentralapnea(:);
    BreathDataTable.Vdr_est = Vdr_est(:);
    BreathDataTable.VAr_est = VAr_est(:);
    BreathDataTable.Error = Error(:);

    %% upsample VE
    VI1 = interp1(BB_t-BB_t(1),VI,0:(BB_t(end)-BB_t(1)),'previous')';
    VIsd = nanstd(detrend(VI1));
    
    %%
    MoreInfo = [1-FitQual(2),VIsd];
    %plot(filter121(detrend(VI1),10))
    
    %% WindowLevelData is here
WinInfo = [LoopGainParameters,LG_QualityInfo(1:3),MoreInfo];
%%
% 1 Time(BB_i_start)
% 2 Time(BB_i_mid)
% 3 Time(BB_i_end)
% 4 VI'
% 5 Vdr_est
% 6 E1'
% 7 E_recover'
% 8 E_Terminate'
% 9 Error
% 10 AR
% 11 meanVIbeforenormalizing+0*AR
% 12 VAr_est
% 13 pos_B
% 14 hypnog_B
% 15 meanVIbeforenormalizing*VI'
% 16 DeltaPes
% 17 DeltaPmus
% 18 VIpes
% 19 DeltaEdi
% 20 VIedi
% 21 GGpeak
% 22 GGtonic
% 23 FlowPes_VI
% 24 FlowEdi_VI];
end


%% ObstructBreathScore_V14_1
function [E1,Ecentralapnea,Ecentralhypop] = ObstructBreathScore_V14_1(Veupnea,E,VI,BB_i_start,removecentralhypops,eventslessthanX,plot_on)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Update events based on eupnea, and detection of "detected" central
% apneas... i.e. any breath greater than eupnoea has events cleared,
% any breath at less than 70% of eupnoea is set as an event. Also, any
% breath detected as a central event (depending on criteria) has event
% cleared for purpose of model fitting.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<

EventBreaths=VI;
EventBreaths(E==1)=[];
ClearBreaths=VI;
ClearBreaths(E==0)=[];
VI = VI(:)';

    % Only edit events at the edge of existing manually defined events, and
    % only using the >1 to clear, and <0.7 to mark. for between 0.7 and 1,
    % events left as-is
    E1=E;
    E_Trans=E;
    E_Trans(1:end)=0;
    EventMod=1;
    while EventMod==1
        % Identify the breaths before and after events, check
        % threshold, and mark as events
        E_temp=E1;
        E_Trans(1:end)=0;
        for k=1:length(BB_i_start)-1
            if (E1(k)==1)&&(E1(k+1)==0) % k is breath before event starts; k+1 is the first breath of event
                E_Trans(k)=1;
                E_Trans(k+1)=-1;
            elseif (E1(k)==0)&&(E1(k+1)==1) % k is last breath of event; k+1 is first breath after event
                E_Trans(k)=-1;
                E_Trans(k+1)=1;
            end
        end
        E1((VI<Veupnea*eventslessthanX)&(E_Trans>0))=0; % for breaths identified as the breath before or after scored event (even if it is a single non event breath) and ventilation below threshold, score as an event breath
        
        % Identify first and last breaths of event, check
        % threshold, and clear events as appropriate.
        E_Trans(1:end)=0;
        for k=1:length(BB_i_start)-1
            if (E1(k)==1)&&(E1(k+1)==0) % k is breath before event starts; k+1 is the first breath of event
                E_Trans(k)=E_Trans(k);
                E_Trans(k+1)=E_Trans(k)-1;
            elseif (E1(k)==0)&&(E1(k+1)==1) % k is last breath of event; k+1 is first breath after event
                E_Trans(k)=E_Trans(k)-1;
                E_Trans(k+1)=E_Trans(k)+1;
            end
        end
        E1((VI>Veupnea)&(E_Trans==-1))=1; % For first, or last breath of event, unless they are the same breath (when E_trans=-2, not -1), and ventilation greater than threshold, clear event.
        
        if sum(abs(E_temp-E1))==0 %i.e. no changes to events
            EventMod=0;
        else
            EventMod=1;
        end
    end

EventBreaths1=VI;
EventBreaths1(E1==1)=[];
ClearBreaths1=VI;
ClearBreaths1(E1==0)=[];
N_EB = hist(EventBreaths,0.05:0.10:2);
N_CB = hist(ClearBreaths,0.05:0.10:2);
N_EB1 = hist(EventBreaths1,0.05:0.10:2);
N_CB1 = hist(ClearBreaths1,0.05:0.10:2);
if plot_on==1
    figure
    subplot(1,2,1)
    hold on
    plot(0.05:0.10:2,N_EB,'r')
    plot(0.05:0.10:2,N_CB,'b')
    subplot(1,2,2)
    hold on
    plot(0.05:0.10:2,N_EB1,'r')
    plot(0.05:0.10:2,N_CB1,'b')
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Identify any central events based on the derived VE data, and
% the abdo and thorax RIP data:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Ecentralhypop=0*E;
Ecentralapnea=0*E;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Identify any single breath central's, and any single non-centrals
% (within a central cluster) and remove.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for k=2:(length(VI)-1)
    %if in a sequence of obstructed breaths... (three breaths of
    %event... k is middle breath)
    if ((E(k)==0)&&(E(k-1)==0)&&(E(k+1)==0))
        if Ecentralhypop(k)==1&&Ecentralhypop(k-1)==0&&Ecentralhypop(k+1)==0
            Ecentralhypop(k)=0;
        elseif Ecentralhypop(k)==0&&Ecentralhypop(k-1)==1&&Ecentralhypop(k+1)==1
            Ecentralhypop(k)=1;
        end
        if Ecentralapnea(k)==1&&Ecentralapnea(k-1)==0&&Ecentralapnea(k+1)==0
            Ecentralapnea(k)=0;
        elseif Ecentralapnea(k)==0&&Ecentralapnea(k-1)==1&&Ecentralapnea(k+1)==1
            Ecentralapnea(k)=1;
        end
        
    end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Clear "Events" which have been marked as central events
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<

for k=1:length(BB_i_start)
    if (Ecentralhypop(k)||Ecentralapnea(k))&&removecentralhypops
        E1(k)=1;
        %         elseif (Ecentralapnea(k))&&removecentralaps
        %             E1(k)=1;
    end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

