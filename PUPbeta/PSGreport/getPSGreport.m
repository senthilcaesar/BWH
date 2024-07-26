function A_Summary_T=getPSGreport(SigT,Evts,ChannelsList,settings,fname,CPAPoff)
%% Spike Style
%% Initial Set up

Epochs=SigT.Epochs;
EventsResp=SigT.EventsResp;
EventsAr=SigT.EventsAr;
Time=SigT.Time;
try
    SnoreDB=SigT.SnoreDB;
end
try
    NoxAudio=SigT.NoxAudio;
end

if ~exist('CPAPoff')
    try
        CPAP=SigT.Pmask;
        [CPAPoff,CPAP] = getCPAP(SigT,settings,ChannelsList,0);
    catch
        CPAPoff=ones(size(SigT,1),1);
    end
end

%%
PSGReportfilename = fname;
settings.onlykeepdataafterCPAPswitchedoff=0;
settings.tableformat=0;
settings.writexls=1;
savexlsdatato=[settings.workdir 'Source'];
settings.printpdfviaexcel=0;
settings.savesummaryfigure=1;
settings.Nsubplots=6;

%%

dt=Time(2)-Time(1);
Fs=1/dt;

% added code for position codes depending on protocol
PosChan = find(strcmp(ChannelsList,'Position')==1);
if ~isempty(PosChan)
    PositionRaw=SigT.Position;
    Position = PositionTranslator(settings.positioncodes,settings.positioncodesout,PositionRaw);
else  % DLM added this case to handle NaN ColumnHeads(10), 
      % It artificially sets the entire night to supine
    if length(settings.supinepositioncode) > 1
        supinecode = settings.supinepositioncode(2);
    else
        supinecode = settings.supinepositioncode;
    end
    Position=ones(size(SigT,1),1)*supinecode;
    PositionRaw = Position; %to avoid code to break if no position channel available
end
issupine = PositionSelector(Position,'Supine');

%% CPAP OFF
if settings.onlykeepdataafterCPAPswitchedoff %consider all time leading up to CPAP on as CPAP on
    CPAPoff_original = CPAPoff;

    while 1
        % Generate indices for 1)CPAP switch on TO switch off 
        % and 2) CPAP switch off to switch on
        
        % 1)
        CPAPoff(end) = 0; CPAPoff(1) = 0;
        CPAPoff_diff = diff(CPAPoff);
        CPAPswitchOff_1 = find(CPAPoff_diff == 1);
        CPAPswitchOn_1 = find(CPAPoff_diff == -1);
        durationCPAPon = CPAPswitchOn_1 - CPAPswitchOff_1; %duration bw off and on

        % 2)
        CPAPoff(1) = 1; CPAPoff(end) = 1;
        CPAPon = ~CPAPoff;   
        CPAPon_diff = diff(CPAPon);
        CPAPswitchOn_2 = find(CPAPon_diff == 1);
        CPAPswitchOff_2 = find(CPAPon_diff == -1);
        durationCPAPoff =  CPAPswitchOff_2 - CPAPswitchOn_2; %duration bw off and on

        starti = [CPAPswitchOff_1; CPAPswitchOn_2];
        stopi = [CPAPswitchOn_1; CPAPswitchOff_2];
        durationAll = [durationCPAPon; durationCPAPoff];
        
        [minDuration, min_i] = min(durationAll);
        
        if minDuration > 200*Fs
            break
        end
        
        CPAPoff(starti(min_i):stopi(min_i)) = ...
                ~CPAPoff(round((starti(min_i)+stopi(min_i))/2));
            
    end
    
    % Find end of long period of CPAPon (i.e. not off)    
    CPAPoff_diff = diff(CPAPoff);
    CPAPswitchOff = find(CPAPoff_diff == 1);
    CPAPswitchOn = find(CPAPoff_diff == -1);
    
    durationCPAPon = CPAPswitchOff - CPAPswitchOn; %duration bw off and on
    for ii = 1:length(durationCPAPon)
        if durationCPAPon(ii) > 30*60*Fs && CPAPswitchOn(ii) < length(CPAPoff)/2
            CPAPoff(1:CPAPswitchOff(ii)) = false;
        end
    end
end

%% Import lights off lights on

sleepcodes = [0:3];
startsleepi = find(sum(Epochs==sleepcodes,2)>0,1);
endsleepi = find(sum(Epochs==sleepcodes,2)>0,1,'last');

% Assign lights off/on times if available 
textfilename=[fname(1:end-4) '_lights.txt'];
if exist(textfilename,'file')==2
    'found lightson/lightsoff file'
    lights = textread(textfilename,'%d%*[^\n]','delimiter',' ');
else
    xminprepost=300;
    lefti = startsleepi-round((xminprepost/dt));
    righti = endsleepi+round((xminprepost/dt));
    if lefti<1, lefti=1; end
    if righti>length(Time), righti=length(Time); end
    clear lights
    lights(1) = Time(lefti);
    lights(2) = Time(righti);
end

if isfield(settings,'lightsoffonarestartendfile') & settings.lightsoffonarestartendfile==1
    lights(1) = Time(1);
    lights(2) = Time(end);
end

clear textfilename lefti righti

lightsoffXHz = 0*Time; lightsoffXHz(Time>=lights(1)&Time<=lights(2))=1;

%% SLEEP METRICS Total

NREMEpochs=Epochs>=0&Epochs<3;
REMEpochs=Epochs==3;
AllSleepEpochs=Epochs>=0&Epochs<=3;
WakeEpochs=Epochs==4;

Duration_total=sum((lightsoffXHz==1))*dt/60;
Duration_sleep_min=sum((lightsoffXHz==1&AllSleepEpochs))*dt/60;
Duration_NREM_min=sum((lightsoffXHz==1&NREMEpochs))*dt/60;
Duration_REM_min=sum((lightsoffXHz==1&REMEpochs))*dt/60;
Duration_wake_min=sum((lightsoffXHz==1&WakeEpochs))*dt/60;
SleepEfficiency=Duration_sleep_min/Duration_total;

Duration_total_CPAPoff=sum((lightsoffXHz==1)&CPAPoff)*dt/60;
Duration_sleep_min_CPAPoff=sum((lightsoffXHz==1&AllSleepEpochs&CPAPoff))*dt/60;
Duration_NREM_min_CPAPoff=sum((lightsoffXHz==1&NREMEpochs&CPAPoff))*dt/60;
Duration_REM_min_CPAPoff=sum((lightsoffXHz==1&REMEpochs&CPAPoff))*dt/60;
Duration_wake_min_CPAPoff=sum((lightsoffXHz==1&WakeEpochs&CPAPoff))*dt/60;
SleepEfficiency_CPAPoff=Duration_sleep_min_CPAPoff/Duration_total_CPAPoff;


%% SpO2 data --start from here

SaO2XHz = SigT.SpO2;

% Find Nadir SpO2
if 1 %Nadir SpO2 in sleep only
    lastEpochs=[NaN*zeros(round(30/dt),1);Epochs((1:(length(Epochs)-round(30/dt))))];
    [nadir_desat_sleep,tempi] = min(SaO2XHz((Epochs<=3|lastEpochs<=3)&~isnan(SaO2XHz)));
    temptime=Time((Epochs<=3|lastEpochs<=3)&~isnan(SaO2XHz));
else
    [nadir_desat_sleep,tempi] = min(SaO2XHz(~isnan(SaO2XHz))); %forget sleep as a criteria...
    temptime=Time(~isnan(SaO2XHz));
end

nadir_desat_sleep_time=temptime(tempi);
clear tempi temptime
% Tidy up desat

%Plot
if 1
    if settings.plotfigs
        figure(1111),
        ax1(5)=subplot(settings.Nsubplots,1,1);
        plot(Time,SaO2XHz,nadir_desat_sleep_time,nadir_desat_sleep,'ro');  %,Time,SaO2XHz3,
        set(gca,'FontName','Arial Narrow','FontSize',10,'box','off','XTick',[],'Xcolor',[1 1 1]);
        ylabel('SpO2 %');
        linkaxes(ax1,'x');
    end
end

SpO2 = SaO2XHz;
clear SaO2XHz;

%% CalcODI

%By Ali Azarbarzin and Scott Sands

Include(:,1) = lightsoffXHz==1&~isnan(SpO2)&CPAPoff;
Include(:,2) = lightsoffXHz==1&~isnan(SpO2)&NREMEpochs&CPAPoff; %NREM
Include(:,3) = lightsoffXHz==1&~isnan(SpO2)&NREMEpochs&(issupine)&CPAPoff; %NREMSupine
Include(:,4) = lightsoffXHz==1&~isnan(SpO2)&REMEpochs&CPAPoff; %REM
[ODI3,ODI4] = CalcODI(SpO2,dt,Include);

%% Desat info

%nadir_desat_sleep
%nadir_desat_sleep_time
SpO2mean_wake=mean(SpO2(WakeEpochs&~isnan(SpO2))); %lightsoffXHz==1
nadir_desat_sleep_prct_0_1_5_25_50_75_95_99 = prctile(SpO2(lightsoffXHz==1&AllSleepEpochs&CPAPoff),[0 1 5 25 50 75 95 99 100]);
SpO2mean_sleep=mean(SpO2(lightsoffXHz==1&AllSleepEpochs&~isnan(SpO2)&CPAPoff));
SpO2mean_sup=mean(SpO2(lightsoffXHz==1&AllSleepEpochs&~isnan(SpO2)&(issupine)&CPAPoff));

SpO2below90p_sleep_prct = 100*sum((lightsoffXHz==1&AllSleepEpochs&SpO2<90))/sum((lightsoffXHz==1&AllSleepEpochs&~isnan(SpO2))&CPAPoff);
SpO2below90p_sleep_prct_NREMsup = 100*sum((lightsoffXHz==1&NREMEpochs&(issupine)&SpO2<90&CPAPoff))/sum((lightsoffXHz==1&NREMEpochs&(issupine)&~isnan(SpO2)&CPAPoff));
SpO2below95p_sleep_prct = 100*sum((lightsoffXHz==1&AllSleepEpochs&SpO2<95)&CPAPoff)/sum((lightsoffXHz==1&AllSleepEpochs&~isnan(SpO2)&CPAPoff));

SpO2below70to100p_sleep_prct = 100*sum((lightsoffXHz==1&AllSleepEpochs&SpO2<[70:100]&CPAPoff))/sum((lightsoffXHz==1&AllSleepEpochs&~isnan(SpO2)&CPAPoff));

nadir_desat_sup_prct_0_1_5_25_50_75_95_99 = prctile(SpO2(lightsoffXHz==1&AllSleepEpochs&(issupine)&CPAPoff),[0 1 5 25 50 75 95 99 100]);
nadir_desat_NREMsup_prct_0_1_5_25_50_75_95_99 = prctile(SpO2(lightsoffXHz==1&NREMEpochs&(issupine)&CPAPoff),[0 1 5 25 50 75 95 99 100]);
nadir_desat_REMsup_prct_0_1_5_25_50_75_95_99 = prctile(SpO2(lightsoffXHz==1&REMEpochs&(issupine)&CPAPoff),[0 1 5 25 50 75 95 99 100]);
nadir_desat_NREMlat_prct_0_1_5_25_50_75_95_99 = prctile(SpO2(lightsoffXHz==1&NREMEpochs&~(issupine)&CPAPoff),[0 1 5 25 50 75 95 99 100]);
SpO2mean_NREMsup=mean(SpO2(lightsoffXHz==1&NREMEpochs&(issupine)&~isnan(SpO2)&CPAPoff));
SpO2mean_REMsup=mean(SpO2(lightsoffXHz==1&REMEpochs&(issupine)&~isnan(SpO2)&CPAPoff));
SpO2mean_NREMlat=mean(SpO2(lightsoffXHz==1&NREMEpochs&~(issupine)&~isnan(SpO2)&CPAPoff));

clear A_SpO2data
A_SpO2data = [...
    SpO2mean_wake ...
    nadir_desat_sleep ...
    SpO2below90p_sleep_prct ...
    SpO2mean_sleep ...
    SpO2mean_sup ...
    SpO2mean_NREMsup ...
    SpO2mean_REMsup ...
    SpO2mean_NREMlat ...
    SpO2below90p_sleep_prct_NREMsup ...
    nadir_desat_NREMsup_prct_0_1_5_25_50_75_95_99 ...
    ]';


%% Total AHI 
%uppercode=12; %double check

uppercode=7; %First code to exclude in regular AHI... currently code 7 is "other" i.e. dip in VE without desat/arousal.

starteventindices=round((Evts.Table1.EventStart-Time(1))/dt);

%Events during sleep: Total AHI
N_events_allresp=0;
N_arousals_all=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    if temp_currentepoch>=0&&temp_currentepoch<=3 %All sleep
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp=N_events_allresp+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all=N_arousals_all+1;
        end
    end
end


%Events during NREM
N_events_allresp_NREM=0;
N_arousals_all_NREM=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    if temp_currentepoch>=0&&temp_currentepoch<3 %NREM
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_NREM=N_events_allresp_NREM+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_NREM=N_arousals_all_NREM+1;
        end
    end
end

%Events during REM
N_events_allresp_REM=0;
N_arousals_all_REM=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    if temp_currentepoch==3 %REM
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_REM=N_events_allresp_REM+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_REM=N_arousals_all_REM+1;
        end
    end
end

%AHIs
TotalAHI = N_events_allresp/(Duration_sleep_min_CPAPoff/60);
NREMAHI = N_events_allresp_NREM/(Duration_NREM_min_CPAPoff/60);
REMAHI = N_events_allresp_REM/(Duration_REM_min_CPAPoff/60);
%ArI
TotalArI = N_arousals_all/(Duration_sleep_min_CPAPoff/60);
NREMArI = N_arousals_all_NREM/(Duration_NREM_min_CPAPoff/60);
REMArI = N_arousals_all_REM/(Duration_REM_min_CPAPoff/60);


%% Calculate supine AHI
%Events during sleep

N_events_allresp_sup=0;
N_arousals_all_sup=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if (temp_currentepoch>=0&&temp_currentepoch<=3)&&(pos_currentepoch==1) %All sleep
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_sup=N_events_allresp_sup+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_sup=N_arousals_all_sup+1;
        end
    end
end

%Events during NREM
N_events_allresp_NREM_sup=0;
N_arousals_all_NREM_sup=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if (temp_currentepoch>=0&&temp_currentepoch<3)&&(pos_currentepoch==1) %NREM
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_NREM_sup=N_events_allresp_NREM_sup+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_NREM_sup=N_arousals_all_NREM_sup+1;
        end
    end
end

%Events during REM
N_events_allresp_REM_sup=0;
N_arousals_all_REM_sup=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if temp_currentepoch==3&&(pos_currentepoch==1) %REM
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_REM_sup=N_events_allresp_REM_sup+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_REM_sup=N_arousals_all_REM_sup+1;
        end
    end
end


%duration of sleep
Duration_total_sup=sum(lightsoffXHz==1&(issupine))*dt/60;
Duration_sleep_min_sup=sum((lightsoffXHz==1&(AllSleepEpochs)&(issupine)))*dt/60;
Duration_NREM_min_sup=sum((lightsoffXHz==1&(NREMEpochs)&(issupine)))*dt/60;
Duration_REM_min_sup=sum((lightsoffXHz==1&REMEpochs&(issupine)))*dt/60;
Duration_wake_min_sup=sum((lightsoffXHz==1&WakeEpochs&(issupine)))*dt/60;
SleepEfficiency_sup=Duration_sleep_min_sup/Duration_total_sup;

Duration_total_sup_CPAPoff=sum(lightsoffXHz==1&(issupine)&CPAPoff)*dt/60;
Duration_sleep_min_sup_CPAPoff=sum((lightsoffXHz==1&(AllSleepEpochs)&(issupine)&CPAPoff))*dt/60;
Duration_NREM_min_sup_CPAPoff=sum((lightsoffXHz==1&(NREMEpochs)&(issupine)&CPAPoff))*dt/60;
Duration_REM_min_sup_CPAPoff=sum((lightsoffXHz==1&REMEpochs&(issupine)&CPAPoff))*dt/60;
Duration_wake_min_sup_CPAPoff=sum((lightsoffXHz==1&WakeEpochs&(issupine)&CPAPoff))*dt/60;
SleepEfficiency_sup_CPAPoff=Duration_sleep_min_sup_CPAPoff/Duration_total_sup_CPAPoff;

%AHI supine
TotalAHI_sup = N_events_allresp_sup/(Duration_sleep_min_sup_CPAPoff/60);
NREMAHI_sup = N_events_allresp_NREM_sup/(Duration_NREM_min_sup_CPAPoff/60);
REMAHI_sup = N_events_allresp_REM_sup/(Duration_REM_min_sup_CPAPoff/60);
%ArI supine
TotalArI_sup = N_arousals_all_sup/(Duration_sleep_min_sup_CPAPoff/60);
NREMArI_sup = N_arousals_all_NREM_sup/(Duration_NREM_min_sup_CPAPoff/60);
REMArI_sup = N_arousals_all_REM_sup/(Duration_REM_min_sup_CPAPoff/60);


%% Calculate non-supine AHI
%Events during sleep
N_events_allresp_lat=0;
N_arousals_all_lat=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if (temp_currentepoch>=0&&temp_currentepoch<=3)&&(pos_currentepoch~=1) %All sleep
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_lat=N_events_allresp_lat+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_lat=N_arousals_all_lat+1;
        end
    end
end

%Events during NREM
N_events_allresp_NREM_lat=0;
N_arousals_all_NREM_lat=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if (temp_currentepoch>=0&&temp_currentepoch<3)&&(pos_currentepoch~=1) %NREM
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_NREM_lat=N_events_allresp_NREM_lat+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_NREM_lat=N_arousals_all_NREM_lat+1;
        end
    end
end

%Events during REM
N_events_allresp_REM_lat=0;
N_arousals_all_REM_lat=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if temp_currentepoch==3&&(pos_currentepoch~=1) %REM
        if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
            N_events_allresp_REM_lat=N_events_allresp_REM_lat+1;
        end
        if Evts.Table1.EventCodes(x,1)==1
            N_arousals_all_REM_lat=N_arousals_all_REM_lat+1;
        end
    end
end

%duration of sleep
Duration_total_lat=sum(lightsoffXHz==1&(~issupine))*dt/60;
Duration_sleep_min_lat=sum((lightsoffXHz==1&(AllSleepEpochs)&(~issupine)))*dt/60;
Duration_NREM_min_lat=sum((lightsoffXHz==1&(NREMEpochs)&(~issupine)))*dt/60;
Duration_REM_min_lat=sum((lightsoffXHz==1&REMEpochs&(~issupine)))*dt/60;
Duration_wake_min_lat=sum((lightsoffXHz==1&WakeEpochs&(~issupine)))*dt/60;
SleepEfficiency_lat=Duration_sleep_min_lat/Duration_total_lat;

Duration_total_lat_CPAPoff=sum(lightsoffXHz==1&(~issupine)&CPAPoff)*dt/60;
Duration_sleep_min_lat_CPAPoff=sum((lightsoffXHz==1&(AllSleepEpochs)&(~issupine)&CPAPoff))*dt/60;
Duration_NREM_min_lat_CPAPoff=sum((lightsoffXHz==1&(NREMEpochs)&(~issupine)&CPAPoff))*dt/60;
Duration_REM_min_lat_CPAPoff=sum((lightsoffXHz==1&REMEpochs&(~issupine)&CPAPoff))*dt/60;
Duration_wake_min_lat_CPAPoff=sum((lightsoffXHz==1&WakeEpochs&(~issupine)&CPAPoff))*dt/60;
SleepEfficiency_lat_CPAPoff=Duration_sleep_min_lat_CPAPoff/Duration_total_lat_CPAPoff;

%AHI lateral
TotalAHI_lat = N_events_allresp_lat/(Duration_sleep_min_lat_CPAPoff/60);
NREMAHI_lat = N_events_allresp_NREM_lat/(Duration_NREM_min_lat_CPAPoff/60);
REMAHI_lat = N_events_allresp_REM_lat/(Duration_REM_min_lat_CPAPoff/60);
%ArI lateral
TotalArI_lat = N_arousals_all_lat/(Duration_sleep_min_lat_CPAPoff/60);
NREMArI_lat = N_arousals_all_NREM_lat/(Duration_NREM_min_lat_CPAPoff/60);
REMArI_lat = N_arousals_all_REM_lat/(Duration_REM_min_lat_CPAPoff/60);

% Event types in Supine NREM
clear N_Etypes_events_allresp_NREM_sup;
N_Etypes_events_allresp_NREM_sup=zeros(1,10);
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    if (temp_currentepoch>=0&&temp_currentepoch<3)&&(pos_currentepoch==1) %NREM Supine
        if Evts.Table1.EventCodes(x,1)>0
            N_Etypes_events_allresp_NREM_sup(Evts.Table1.EventCodes(x,1))=N_Etypes_events_allresp_NREM_sup(Evts.Table1.EventCodes(x,1))+1;
        end
    end
end
Etypes_events_allresp_NREM_sup_per_hour=N_Etypes_events_allresp_NREM_sup/(Duration_NREM_min_sup_CPAPoff/60);



% Event types in TOTAL
clear N_Etypes_events_allresp_TOTAL;
N_Etypes_events_allresp_TOTAL=zeros(1,10);
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    if (temp_currentepoch>=0&&temp_currentepoch<=3) %SLEEP
        if Evts.Table1.EventCodes(x,1)>0
            N_Etypes_events_allresp_TOTAL(Evts.Table1.EventCodes(x,1))=N_Etypes_events_allresp_TOTAL(Evts.Table1.EventCodes(x,1))+1;
        end
    end
end
Etypes_events_allresp_TOTAL_per_hour=N_Etypes_events_allresp_TOTAL/(Duration_sleep_min_CPAPoff/60);

%% First half vs second half of night
if 0
    time50percentTST=Time(round(length(Time)/2)); %in seconds according to Spike time
else %split by NREM supine midpoint -- or split by NREM sleep regardless of position.
    % figure(3);plot(Time,cumsum(Epochs>=0&Epochs<3)/max(cumsum(Epochs>=0&Epochs<3)));
    if 0 %split by NREMsupine midtime
        temp=cumsum(CPAPoff&lightsoffXHz==1&Epochs>=0&Epochs<3&(issupine))/max(cumsum(CPAPoff&Epochs>=0&Epochs<3&(issupine)));
    else %split by NREM midtime, changed 2/18/2015
        temp=cumsum(CPAPoff&lightsoffXHz==1&Epochs>=0&Epochs<3)/max(cumsum(CPAPoff&Epochs>=0&Epochs<3));
    end
    time50percentTST=Time(find(temp>0.5,1));
end
%Events during sleep: Total AHI
N_events_allresp_1sthalfnight=0;
N_arousals_all_1sthalfnight=0;
N_events_allresp_2ndhalfnight=0;
N_arousals_all_2ndhalfnight=0;
N_events_allresp_NREM_sup_1sthalfnight=0;
N_arousals_all_NREM_sup_1sthalfnight=0;
N_events_allresp_NREM_sup_2ndhalfnight=0;
N_arousals_all_NREM_sup_2ndhalfnight=0;
for x=1:length(Evts.Table1.EventStart)
    if Evts.Table1.EventStart(x)<lights(1)||Evts.Table1.EventStart(x)>lights(2)||CPAPoff(starteventindices(x))==0
        continue
    end
    time_index = round((Evts.Table1.EventStart(x)-Time(1))/dt+1);
    temp_currentepoch=Epochs(time_index);
    pos_currentepoch=(Position(time_index));
    %currenttime=Time(time_index);
    if Time(time_index)<time50percentTST %first half of night
        if temp_currentepoch>=0&&temp_currentepoch<=3 %All sleep
            if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
                N_events_allresp_1sthalfnight=N_events_allresp_1sthalfnight+1;
            end
            if Evts.Table1.EventCodes(x,1)==1
                N_arousals_all_1sthalfnight=N_arousals_all_1sthalfnight+1;
            end
        end
        if (temp_currentepoch>=0&&temp_currentepoch<3)&&(pos_currentepoch==1) %NREM SUP
            if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
                N_events_allresp_NREM_sup_1sthalfnight=N_events_allresp_NREM_sup_1sthalfnight+1;
            end
            if Evts.Table1.EventCodes(x,1)==1
                N_arousals_all_NREM_sup_1sthalfnight=N_arousals_all_NREM_sup_1sthalfnight+1;
            end
        end
    else %2nd half of night
        if temp_currentepoch>=0&&temp_currentepoch<=3 %All sleep
            if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
                N_events_allresp_2ndhalfnight=N_events_allresp_2ndhalfnight+1;
            end
            if Evts.Table1.EventCodes(x,1)==1
                N_arousals_all_2ndhalfnight=N_arousals_all_2ndhalfnight+1;
            end
        end
        if (temp_currentepoch>=0&&temp_currentepoch<3)&&(pos_currentepoch==1) %NREM SUP
            if Evts.Table1.EventCodes(x,1)>=2&&Evts.Table1.EventCodes(x,1)<uppercode
                N_events_allresp_NREM_sup_2ndhalfnight=N_events_allresp_NREM_sup_2ndhalfnight+1;
            end
            if Evts.Table1.EventCodes(x,1)==1
                N_arousals_all_NREM_sup_2ndhalfnight=N_arousals_all_NREM_sup_2ndhalfnight+1;
            end
        end
    end
end

%all CPAPoff
Duration_min_1sthalfnight=sum((CPAPoff&lightsoffXHz==1&(Epochs>=0)&(Epochs<=4)&(Time<=time50percentTST)))*dt/60;
Duration_min_2ndhalfnight=sum((CPAPoff&lightsoffXHz==1&(Epochs>=0)&(Epochs<=4)&(Time>time50percentTST)))*dt/60;
Duration_sleep_min_1sthalfnight=sum((CPAPoff&lightsoffXHz==1&(AllSleepEpochs)&(Time<=time50percentTST)))*dt/60;
Duration_sleep_min_2ndhalfnight=sum((CPAPoff&lightsoffXHz==1&(AllSleepEpochs)&(Time>time50percentTST)))*dt/60;
Duration_sleep_min_1sthalfnight_NREMsup=sum((CPAPoff&lightsoffXHz==1&(NREMEpochs)&(issupine)&(Time<=time50percentTST)))*dt/60;
Duration_sleep_min_2ndhalfnight_NREMsup=sum((CPAPoff&lightsoffXHz==1&(NREMEpochs)&(issupine)&(Time>time50percentTST)))*dt/60;

AHI_1stHN_total=N_events_allresp_1sthalfnight/Duration_sleep_min_1sthalfnight*60;
AHI_2ndHN_total=N_events_allresp_2ndhalfnight/Duration_sleep_min_2ndhalfnight*60;
ArI_1stHN_total=N_arousals_all_1sthalfnight/Duration_sleep_min_1sthalfnight*60;
ArI_2ndHN_total=N_arousals_all_2ndhalfnight/Duration_sleep_min_2ndhalfnight*60;
AHI_1stHN_NREMsup=N_events_allresp_NREM_sup_1sthalfnight/Duration_sleep_min_1sthalfnight_NREMsup*60;
AHI_2ndHN_NREMsup=N_events_allresp_NREM_sup_2ndhalfnight/Duration_sleep_min_2ndhalfnight_NREMsup*60;
ArI_1stHN_NREMsup=N_arousals_all_NREM_sup_1sthalfnight/Duration_sleep_min_1sthalfnight_NREMsup*60;
ArI_2ndHN_NREMsup=N_arousals_all_NREM_sup_2ndhalfnight/Duration_sleep_min_2ndhalfnight_NREMsup*60;

%% Percent Stable Breathing
%Outputs are SBpercent SBpercentNREM, SBpercentNREMSupine
minlength=(60/dt)*[1,2,3,5,10]; %Count 3 minute periods of stable/event-free sleep...

%0=no event, 1=has event, -1=excluded.
NonEventSleep=-ones(1,length(Time));
NonEventSleep((AllSleepEpochs)&CPAPoff&lightsoffXHz==1)=0;
NonEventSleep(EventsResp>0|EventsAr>0&CPAPoff&lightsoffXHz==1)=1; %overwrite

%Actual script
den=sum(NonEventSleep~=-1);
diffleft=[NaN diff(NonEventSleep)];
diffright=[diff(NonEventSleep) NaN];
start_I1 = NonEventSleep==0&diffleft~=0; %no event, and change just happened
end_I1 = NonEventSleep==0&diffright~=0; %no event, and change about to happen
start_I2 = find(start_I1>0); %list of start indices
end_I2 = find(end_I1>0); %list of end indices
lengthsofeventfreeperiods=end_I2-start_I2+1;
num=zeros(1,length(minlength));
for i=1:length(minlength)
    num(i)=sum(lengthsofeventfreeperiods(lengthsofeventfreeperiods>minlength(i)));
end
%Result array
SBpercent=100*num/den;

%0=no event, 1=has event, -1=excluded.
NonEventNREMSleep=-ones(1,length(Time));
NonEventNREMSleep(NREMEpochs)=0;
NonEventNREMSleep(EventsResp>0|EventsAr>0)=1; %overwrite

%Actual script
den=sum(NonEventNREMSleep~=-1);
diffleft=[NaN diff(NonEventNREMSleep)];
diffright=[diff(NonEventNREMSleep) NaN];
start_I1 = NonEventNREMSleep==0&diffleft~=0; %no event, and change just happened
end_I1 = NonEventNREMSleep==0&diffright~=0; %no event, and change about to happen
start_I2 = find(start_I1>0); %list of start indices
end_I2 = find(end_I1>0); %list of end indices
lengthsofeventfreeperiods=end_I2-start_I2+1;
clear num;
for i=1:length(minlength)
    num(i)=sum(lengthsofeventfreeperiods(lengthsofeventfreeperiods>minlength(i)));
end
%Result array
SBpercentNREM=100*num/den;



%0=no event, 1=has event, -1=excluded.
NonEventNREMSupineSleep=-ones(1,length(Time));
NonEventNREMSupineSleep((NREMEpochs)&(issupine))=0;
NonEventNREMSupineSleep(EventsResp>0|EventsAr>0)=1; %overwrite

%Actual script
den=sum(NonEventNREMSupineSleep~=-1);
diffleft=[NaN diff(NonEventNREMSupineSleep)];
diffright=[diff(NonEventNREMSupineSleep) NaN];
start_I1 = NonEventNREMSupineSleep==0&diffleft~=0; %no event, and change just happened
end_I1 = NonEventNREMSupineSleep==0&diffright~=0; %no event, and change about to happen
start_I2 = find(start_I1>0); %list of start indices
end_I2 = find(end_I1>0); %list of end indices
lengthsofeventfreeperiods=end_I2-start_I2+1;
clear num;
for i=1:length(minlength)
    num(i)=sum(lengthsofeventfreeperiods(lengthsofeventfreeperiods>minlength(i)));
end
%Result array
SBpercentNREMSupine=100*num/den;

%% Stages and Times (CPAP on included)
S1duration=dt*sum(lightsoffXHz==1&Epochs==2)/60;
S2duration=dt*sum(lightsoffXHz==1&Epochs==1)/60;
S3duration=dt*sum(lightsoffXHz==1&(Epochs==0))/60;

%% Event information
Eventinfolist=zeros(length(Evts.Table1.EventStart),7);
Eventinfolist(:,1)=Evts.Table1.EventCodes;
Eventinfolist(:,2)=Evts.Table1.EventStart;
Eventinfolist(:,3)=Evts.Table1.EventStart+Evts.Table1.EventDuration; % end of events
Eventinfolist(:,4)=Evts.Table1.EventDuration;
Eventinfolist(:,5)=interp1(Time,Epochs,Eventinfolist(:,2),'nearest','extrap');
Eventinfolist(:,6)=interp1(Time,Epochs,Eventinfolist(:,3),'nearest','extrap');
Eventinfolist(:,7)=interp1(Time,Position,Eventinfolist(:,2),'nearest','extrap');
%AR|ApO|ApC|HypO|M|HypC|HypO2%|HypC2%|HypOx|HypCx|
tempI=(Eventinfolist(:,1)==2)|(Eventinfolist(:,1)==4)&(Eventinfolist(:,5)>=0&Eventinfolist(:,5)<3)&(Eventinfolist(:,7)==1);
tempI2=(Eventinfolist(:,1)==3)|(Eventinfolist(:,1)==6)&(Eventinfolist(:,5)>=0&Eventinfolist(:,5)<3)&(Eventinfolist(:,7)==1);
tempIall=(Eventinfolist(:,1)>1)&(Eventinfolist(:,5)>=0&Eventinfolist(:,5)<3)&(Eventinfolist(:,7)==1);
OAHdurations=Eventinfolist(tempI,4);
CAHdurations=Eventinfolist(tempI2,4);

histbinwidth=5; upper=60; bins=0:histbinwidth:upper; Xticks=[0:15:upper];
H2=hist(OAHdurations,bins);
H3=hist(CAHdurations,bins);
subplot(settings.Nsubplots,4,(settings.Nsubplots-1)*4+3);
bar(bins(1:end-1),H2(1:end-1)/sum(H2+H3)+H3(1:end-1)/sum(H2+H3),'BarWidth',1,'EdgeColor','none');
hold('on');
bar(bins(1:end-1),H3(1:end-1)/sum(H2+H3),'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]);
set(gcf,'color',[1 1 1]);
set(gca,'XTick',Xticks,'FontName','Arial Narrow','FontSize',10,'box','off','XLim',[-histbinwidth/2 upper+histbinwidth/2]);
ylabel('Proportion of Events'); xlabel('Event Duration (s)');

[temp,temp_i]=max(H2(1:end-1));
modeduration_obs=bins(temp_i);
[temp,temp_i]=max(H3(1:end-1));
modeduration_cen=bins(temp_i);


eventdurationmedian_obs = median(Eventinfolist(tempI,4)); % NREM Supine
eventdurationmedian_cen = median(Eventinfolist(tempI2,4));
eventdurationmedian_obscen = median(Eventinfolist(tempIall,4));

%% Interevent times: ALL
histbinwidth=5; upper=90; bins=0:histbinwidth:upper;
Xticks=[0:30:upper];
%0=no event, 1=in resp event, -1=excluded.
Temp=-ones(1,length(Time));
Temp((NREMEpochs)&(Position==1)&CPAPoff&lightsoffXHz==1)=0; % NREM Supine
Temp(EventsResp>0&CPAPoff&lightsoffXHz==1)=1; %overwrite
%Temp=[-1 -1 -1 -1 -1 -1 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0  1 1 1 1 1 0 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 -1 ];
%Actual script
diffleft=[NaN diff(Temp)];
start_I1 = (Temp==1&diffleft~=0); %no event, and change just happened
start_I2 = find(start_I1>0); %list of start indices
durationsbetweenevents=dt*diff(start_I2);
if length(start_I2)>1
    clear temp3
    for i=1:(length(start_I2)-1)
        temp3(i)=1-length(find(Temp(start_I2(i):start_I2(i+1))==-1,1));
    end
    durationsbetweenevents=durationsbetweenevents(temp3>0);
end
H=hist(durationsbetweenevents,bins);

[temp,modecyclingperiod_i]=max(H(1:end-1));
modecyclingperiod=bins(modecyclingperiod_i);
mediancyclingperiod=median(durationsbetweenevents(durationsbetweenevents<upper));
meancyclingperiod=mean(durationsbetweenevents(durationsbetweenevents<upper));


% Interevent times: Obstructive

%0=no event, 1=in resp event, -1=excluded.
Temp=-ones(1,length(Time));
Temp((NREMEpochs)&(Position==1))=0;
Temp(EventsResp==2|EventsResp==4)=1; %overwrite / obstructive only...

%Actual script
diffleft=[NaN diff(Temp)];
start_I1 = (Temp==1&diffleft~=0); %no event, and change just happened
start_I2 = find(start_I1>0); %list of start indices
durationsbetweenevents=dt*diff(start_I2);
if length(start_I2)>1
    clear temp3
    for i=1:(length(start_I2)-1)
        temp3(i)=1-length(find(Temp(start_I2(i):start_I2(i+1))==-1,1));
    end
    durationsbetweenevents=durationsbetweenevents(temp3>0);
end

H1=hist(durationsbetweenevents,bins);

[temp,modecyclingperiod_i]=max(H1(1:end-1));
modecyclingperiod_obs=bins(modecyclingperiod_i);
mediancyclingperiod_obs=median(durationsbetweenevents(durationsbetweenevents<upper));
meancyclingperiod_obs=mean(durationsbetweenevents(durationsbetweenevents<upper));

% Interevent times: Central

%0=no event, 1=in resp event, -1=excluded.
Temp=-ones(1,length(Time));
Temp((NREMEpochs)&(Position==1))=0;
Temp(EventsResp==3|EventsResp==6)=1; %overwrite

%Actual script
diffleft=[NaN diff(Temp)];
start_I1 = (Temp==1&diffleft~=0); %no event, and change just happened
start_I2 = find(start_I1>0); %list of start indices
durationsbetweenevents=dt*diff(start_I2);
if length(start_I2)>1
    clear temp3
    for i=1:(length(start_I2)-1)
        temp3(i)=1-length(find(Temp(start_I2(i):start_I2(i+1))==-1,1)); %#ok<*SAGROW>
    end
    durationsbetweenevents=durationsbetweenevents(temp3>0);
end

H=hist(durationsbetweenevents,bins); %central events
%subplot(Nsubplots,3,(Nsubplots-1)*3+3);


subplot(settings.Nsubplots,4,(settings.Nsubplots-1)*4+4);
bar(bins(1:end-1),(H1(1:end-1)/sum(H1)),'BarWidth',1,'EdgeColor','none');
set(gcf,'color',[1 1 1]);
set(gca,'XTick',Xticks,'FontName','Arial Narrow','FontSize',10,'box','off','XLim',[-histbinwidth/2 upper+histbinwidth/2]);
ylabel('Proportion of Events'); xlabel('Inter-event Duration (s)');

hold('on')
bar(bins(1:end-1),(H(1:end-1))/(sum(H1)),'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]);
set(gcf,'color',[1 1 1]);
set(gca,'XTick',Xticks,'FontName','Arial Narrow','FontSize',10,'box','off','XLim',[-histbinwidth/2 upper+histbinwidth/2]);
ylabel('Proportion of Events'); xlabel('Inter-event Duration (s)');

[temp,modecyclingperiod_i]=max(H(1:end-1));
modecyclingperiod_cen=bins(modecyclingperiod_i);
mediancyclingperiod_cen=median(durationsbetweenevents(durationsbetweenevents<upper));
meancyclingperiod_cen=mean(durationsbetweenevents(durationsbetweenevents<upper));

%%
if 1&&exist('SnoreDB','var') 
      %SnoreRMSsq = 10.^((SnoreDB/10))*(0.00002^2);
      %SnoreDB = 10*(log10(SnoreRMSsq/0.00002^2))
       criteria = (sum(Epochs==sleepcodes,2)>0) & SnoreDB>0 & ~isnan(SnoreDB);
       if 0
        criteria = EventsResp~=2&(sum(Epochs==sleepcodes,2)>0);
       end
       %criteria = EventsResp~=2&(sum(Epochs==sleepcodes,2)>0);
       %criteria = EventsResp==2&(sum(Epochs==sleepcodes,2)>0);
       SnoreRMSsq = 10.^((SnoreDB/20))*(0.00002);
       SnoreRMS_sleep = SnoreRMSsq(criteria==1);
       
       SnoreRMSsq = 10.^((SnoreDB/10))*(0.00002^2);
       SnoreRMSsq_sleep = SnoreRMSsq(criteria==1);
       %histogram(SnoreDB(criteria==1),20);
       SnoreDBprctiles = prctile(SnoreDB(criteria==1),[0:100]');
       if 0
            mean_SnoreRMSsq_sleep = nanmean(SnoreRMS_sleep);
       else
            mean_SnoreRMSsq_sleep = nanmean(SnoreRMS_sleep.^2).^0.5;
            mean_SnoreRMSsq_sleep = nanmean(SnoreRMSsq_sleep);
       end
       %mean_SnoreDB_sleep_ = 20*(log10(mean_SnoreRMSsq_sleep/0.00002));
       mean_SnoreDB_sleep = 10*(log10(mean_SnoreRMSsq_sleep/(0.00002^2)));
       SnoreTemp = [mean_SnoreDB_sleep SnoreDBprctiles([1 51 76 91 100 101])']
       disp('1');
else
    mean_SnoreDB_sleep=NaN;
    SnoreDBprctiles=ones(101,1)*NaN;
end
%%
if 1&&exist('NoxAudio','var') 
      %SnoreRMSsq = 10.^((SnoreDB/10))*(0.00002^2);
      %SnoreDB = 10*(log10(SnoreRMSsq/0.00002^2))
       criteria = (sum(Epochs==sleepcodes,2)>0) & NoxAudio>0 & ~isnan(NoxAudio);
       if 0
        criteria = EventsResp~=2&(sum(Epochs==sleepcodes,2)>0);
       end
       %criteria = EventsResp~=2&(sum(Epochs==sleepcodes,2)>0);
       %criteria = EventsResp==2&(sum(Epochs==sleepcodes,2)>0);
       
        SnoreDBNox = 10*log10((NoxAudio).^2./(0.00002^2)); %10*(log10(mean_SnoreRMSsq_sleep/(0.00002^2)));
        SnoreRMSsqNox = 10.^((SnoreDBNox/10))*(0.00002^2);
    
%        SnoreRMSsq = 10.^((SnoreDB/20))*(0.00002);
%        SnoreRMS_sleep = SnoreRMSsq(criteria==1);
       
       %SnoreRMSsq = 10.^((SnoreDB/10))*(0.00002^2);
       SnoreRMSsq_sleepNox = SnoreRMSsqNox(criteria==1);
       %histogram(SnoreDB(criteria==1),20);
       SnoreDBprctilesNox = prctile(SnoreDBNox(criteria==1),[0:100]');
%        if 0
%             mean_SnoreRMSsq_sleep = nanmean(SnoreRMS_sleep);
%        else
%             mean_SnoreRMSsq_sleep = nanmean(SnoreRMS_sleep.^2).^0.5;
            mean_SnoreRMSsq_sleepNox = nanmean(SnoreRMSsq_sleepNox);
%        end
       %mean_SnoreDB_sleep_ = 20*(log10(mean_SnoreRMSsq_sleep/0.00002));
       mean_SnoreDB_sleepNox = 10*(log10(mean_SnoreRMSsq_sleepNox/(0.00002^2)));
       SnoreTempNox = [mean_SnoreDB_sleepNox SnoreDBprctilesNox([1 51 76 91 100 101])']
       disp('1');
else
    mean_SnoreDB_sleepNox=NaN;
    SnoreDBprctilesNox=ones(101,1)*NaN;
end

%% Data list

A_Summary1 = [...
    Duration_total, ...
    Duration_sleep_min, ...
    Duration_NREM_min, ...
    Duration_REM_min, ...
    Duration_wake_min, ...
    SleepEfficiency, ...
    TotalAHI, ...
    NREMAHI, ...
    REMAHI, ...
    TotalArI, ...
    NREMArI, ...
    REMArI, ...
    Duration_total_sup, ...
    Duration_sleep_min_sup, ...
    Duration_NREM_min_sup, ...
    Duration_REM_min_sup, ...
    Duration_wake_min_sup, ...
    SleepEfficiency_sup, ...
    TotalAHI_sup, ...
    NREMAHI_sup, ...
    REMAHI_sup, ...
    TotalArI_sup, ...
    NREMArI_sup, ...
    REMArI_sup, ...
    Duration_total_lat, ...
    Duration_sleep_min_lat, ...
    Duration_NREM_min_lat, ...
    Duration_REM_min_lat, ...
    Duration_wake_min_lat, ...
    SleepEfficiency_lat, ...
    TotalAHI_lat, ...
    NREMAHI_lat, ...
    REMAHI_lat, ...
    TotalArI_lat, ...
    NREMArI_lat, ...
    REMArI_lat, ...
    Etypes_events_allresp_NREM_sup_per_hour, ...
    AHI_1stHN_total, ...
    AHI_2ndHN_total, ...
    ArI_1stHN_total, ...
    ArI_2ndHN_total, ...
    AHI_1stHN_NREMsup, ...
    AHI_2ndHN_NREMsup, ...
    ArI_1stHN_NREMsup, ...
    ArI_2ndHN_NREMsup, ...
    SBpercent, ...
    SBpercentNREMSupine ...
    ]';

More = [...
    S1duration,...
    S2duration,...
    S3duration,...
    modecyclingperiod,...
    mediancyclingperiod,...
    modecyclingperiod_obs,...
    mediancyclingperiod_obs,...
    modecyclingperiod_cen,...
    mediancyclingperiod_cen,...
    Etypes_events_allresp_TOTAL_per_hour, ...
    eventdurationmedian_obs, ... %NREM SUPINE
    eventdurationmedian_cen, ... %NREM SUPINE
    eventdurationmedian_obscen, ... %NREM SUPINE
    ODI3, ...
    ODI4 ...
    mean_SnoreDB_sleep, ...
    SnoreDBprctiles(2), ...
    SnoreDBprctiles(50+1), ...
    SnoreDBprctiles(75+1), ...
    SnoreDBprctiles(90+1), ...
    SnoreDBprctiles(95+1), ...
    SnoreDBprctiles(99+1), ...
    SnoreDBprctiles(100+1), ...
    mean_SnoreDB_sleepNox, ...
    SnoreDBprctilesNox(2), ...
    SnoreDBprctilesNox(50+1), ...
    SnoreDBprctilesNox(75+1), ...
    SnoreDBprctilesNox(90+1), ...
    SnoreDBprctilesNox(95+1), ...
    SnoreDBprctilesNox(99+1), ...
    SnoreDBprctilesNox(100+1), ...
    ]';

A_Summary = [A_Summary1;A_SpO2data;More];

A_Summary_T = A_Summary';

% Same as above, but table format
if settings.tableformat
A_Summary_T = table(...
    Duration_total_CPAPoff, ...
    Duration_sleep_min_CPAPoff, ...
    Duration_NREM_min_CPAPoff, ...
    Duration_REM_min_CPAPoff, ...
    Duration_wake_min_CPAPoff, ...
    SleepEfficiency_CPAPoff, ...
    TotalAHI, ...
    NREMAHI, ...
    REMAHI, ...
    TotalArI, ...
    NREMArI, ...
    REMArI, ...
    Duration_total_sup_CPAPoff, ...
    Duration_sleep_min_sup_CPAPoff, ...
    Duration_NREM_min_sup_CPAPoff, ...
    Duration_REM_min_sup_CPAPoff, ...
    Duration_wake_min_sup_CPAPoff, ...
    SleepEfficiency_sup_CPAPoff, ...
    TotalAHI_sup, ...
    NREMAHI_sup, ...
    REMAHI_sup, ...
    TotalArI_sup, ...
    NREMArI_sup, ...
    REMArI_sup, ...
    SpO2mean_wake, ...
    nadir_desat_sleep, ...
    SpO2below90p_sleep_prct, ...
    SpO2mean_sleep, ...
    SpO2mean_sup, ...
    SpO2mean_NREMsup, ...
    SpO2mean_REMsup, ...
    SpO2mean_NREMlat, ...
    SpO2below90p_sleep_prct_NREMsup, ...
    S1duration,...
    S2duration,...
    S3duration...
    );
end

%% Finish Plot
figure(1111)
%set(gcf,'color',[1 1 1])
subplot(settings.Nsubplots,4,(settings.Nsubplots-1)*4+1)
X=[Duration_wake_min S1duration S2duration S3duration Duration_REM_min];
X=X/sum(X)*100;
bar(X,'EdgeColor','none');
ylabel('%Night')
set(gca,'XTickLabel',{'W', '1', '2', '3','R'})
set(gca,'FontName','Arial Narrow','FontSize',8,'box','off');
% colormap jet
subplot(settings.Nsubplots,4,(settings.Nsubplots-1)*4+2)
X=Etypes_events_allresp_NREM_sup_per_hour(2:6);
bar([X(1) X(3) X(4) X(2) X(5)],'EdgeColor','none');
ylabel('/hour')
set(gca,'XTickLabel',{'OAI','OHI','MAI','CAI','CHI'})
% colormap jet
set(gca,'FontName','Arial Narrow','FontSize',8,'box','off');

%% Write to xls

if settings.writexls
    try
    statustemp = copyfile([settings.codedir '\PSGReport\ReportTemplate3.xls'],[savexlsdatato '\' PSGReportfilename '.xls'])
    cell0='AJ89'; %for filename
    cell1='AJ90'; %for analysis
    cell2='AM90'; %for events lists
    xlswrite([savexlsdatato '\' PSGReportfilename '.xls'],{PSGReportfilename},'PSGReport',cell0);
    xlswrite([savexlsdatato '\' PSGReportfilename '.xls'],A_Summary,'PSGReport',cell1);
    xlswrite([savexlsdatato '\' PSGReportfilename '.xls'],Eventinfolist,'PSGReport',cell2);
    catch me
        disp('Failed write to Excel:')
        disp(me.message)
    end
end

   %% Print via Excel
    
   if settings.printpdfviaexcel
       Excel = actxserver('excel.application');
       Excel.visible = 1;
       Workbooks = Excel.Workbooks;
       Excel.Visible=1;
       Workbook=Workbooks.Open([savexlsdatato '\' PSGReportfilename '.xls']);
       Excel.ActiveWorkbook.PrintOut(1,1,1,'False','Adobe PDF');
       Excel.Quit;
   end
%%
if settings.savesummaryfigure==1
    suffix='_hypnogram';
    set(gcf,'Position',[278    39   915   958]);
    %saveas(1,[savexlsdatato filename suffix],'fig')
    saveas(1111,[savexlsdatato '\' PSGReportfilename suffix],'jpeg')
end

%% Clear large variables
S = whos;
for i=1:length(S)
    try
        if S(i).bytes>200000
            clear(S(i).name);
        end
    catch me
        %not all variables have properties name and bytes
    end
end
clear S

%% Save data to temp
if 0
    save([savexlsdatato '\' PSGReportfilename '_PSGReportInfo.mat']);
end

