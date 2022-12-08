function NoxRIPtoSpike3()

%Using new Nox reader, and broader sync range (15 s) since 7/7/2018

%Before 9/19/17, K was selected based on an average, but this should really have
%been a average of log values. Corrected on 9/19/17. Rerun those prior to
%this date for optimal K calculation.

%User guide:
%0. Before use: Open Nox files in Noxturnal. Go to respiratory tab. Exit.
%1. User selects spike file.
%2. Nox RIP folder(s) inside the spike file directory is automatically found.
%3. New Spike and .mat files with RIP data are output to the same folder.

%Interpretation:
%Equation is Vol = alpha*[Th/(k^0.5) + Ab*(k^0.5)] = Vth + Vab, equal to
%Vol = [alpha/(k^0.5)]*[Th + Ab*k];

%Outputs:
%Overnight (moving-time) median k and alpha values are used for the exported Vab and Vth
%ThNoxRIP/AbNoxRIP are raw inductance in H
%k and alpha are time varying parameters.

%CED library not loading/working in Matlab R2017a [R2016b working fine]
%Nox functions also not loading/working in Matlab R2017a [R2016b working
%fine]
%%
warning('off');
close all;
minfigs=1;

%% Find Spike file, have Nox directory inside
%default = 'J:\PEOPLE\FACULTY\SANDS\O2PSG\_Scored\*.smr';
%default = 'C:\Users\szs88\Dropbox (Partners HealthCare)\MAD-OX\STUDIES\ToScore\InPreparation\*.smr'
default = 'C:\Users\szs88\Dropbox (Partners HealthCare)\MAD-OX\STUDIES\Scored\*.*'

% [FileName,PathName] = uigetfile([default],'Select Spike file with flow data');
% file=[PathName FileName];

a_fnamelastpath = 'lastpathtorec.n2m'; %file stores location of most recently opened directory
if exist(a_fnamelastpath,'file')==2
    fid = fopen(a_fnamelastpath);
    txt = textscan(fid,'%s','delimiter','\n');
    default = txt{1};
end
[FileName,PathName] = uigetfile([default],'Select Spike file with flow data');
file=[PathName FileName];

fid = fopen(a_fnamelastpath,'wt');
fprintf(fid, '%s', [PathName '*.*']);
fclose(fid);

%% Load spike library

addpath(genpath(cd));
cedpath = [cd '\CEDS64ML2017'];

%%
addpath( cedpath ); % so CEDS64LoadLib.m is found
CEDS64LoadLib(cedpath); % load ceds64int.dll


%% channels to load
clear channelnameoptions
flowtype = input('Flow [enter 1] or Pnasal [enter 2]: ');
Pnasallist={'Pnasal','PNasal','Pmask','PMask'}; %'Vflow'
if flowtype==1
    channelnameoptions.Flow={'Flow','Vflow'}; %'Vflow'
elseif flowtype==2
    channelnameoptions.Pnasal=Pnasallist; %'Vflow'
end
if isempty(flowtype)
    disp('Error: no flow type selected');
end
desiredchannel = char(fieldnames(channelnameoptions));

forcesyncPmask = 1; %but use Flow (desiredchannel) for calibration
if forcesyncPmask
channelnameoptions2.Pnasal = Pnasallist;
syncchannel = char(fieldnames(channelnameoptions2));
else
syncchannel = desiredchannel;    
end

%% Open File

fhand1 = CEDS64Open(file);
if fhand1==-21 || fhand1==-1
    disp('Spike file is already open, close it and then try again');
end
if (fhand1 <= 0); unloadlibrary ceds64int; return; end
maxChan=CEDS64MaxChan(fhand1);

%% Get full channel list
clear titles
for ch=1:maxChan
    [~,titles{ch}] = CEDS64ChanTitle( fhand1, ch);
end
SpikeBlank = find(strcmp(titles,'')==1);

%% Search for known channel names
clear sTitle iOk iType
channelnamestemp=fieldnames(channelnameoptions);
Channels=NaN*ones(1,size(channelnamestemp,1));
foundamatch=zeros(1,size(channelnamestemp,1));
for ii=1:length(channelnamestemp)
    temp=eval(['channelnameoptions.' char(channelnamestemp(ii))]);
    for n=1:length(temp)
        for ch=1:maxChan
            [ iOk, sTitle ] = CEDS64ChanTitle( fhand1, ch);
            %[ iType ] = CEDS64ChanType( fhand1, ch ); %0=channel unused, 1=Adc, 2=EventFall,
            %3=EventRise, 4=EventBoth, 5=Marker, 6=WaveMark, 7=RealMark, 8=TextMark,9=RealWave.
            if strcmp(sTitle,temp{n})
                foundamatch(ii)=1;
                disp(['found ' temp{n}])
                if foundamatch(ii)
                    Channels(ii)=ch;
                end
                break
            end
        end
        if foundamatch(ii)
            break
        end
    end
end

%% Get Flow signal data
ch = Channels(find(strcmp(channelnamestemp,desiredchannel)));

%% Load signal
Fs = CEDS64IdealRate(fhand1,ch);
maxTimeTicks = CEDS64ChanMaxTime( fhand1, ch)+2; % +1 so the read gets the last point
NSecs = CEDS64TicksToSecs( fhand1, maxTimeTicks);
%Napprox1=Fs*NSecs;

TimeIntervalInTicks=100000000; %10000000=20 s
TimeIntervalInTicks=min(TimeIntervalInTicks,maxTimeTicks);
[~,testvalues] = CEDS64ReadWaveF(fhand1,ch,TimeIntervalInTicks,1,TimeIntervalInTicks);
Nsamples = length(testvalues);
tickspersample = TimeIntervalInTicks/Nsamples;
Napprox = round(maxTimeTicks/tickspersample);
ChannelsSig=NaN*zeros(Napprox,1);

%% Load signal section by section, avoids overload
tic
Nwindows = ceil(maxTimeTicks/TimeIntervalInTicks);
disp(['importing Flow/Pnasal (~25 s)'])
for i=1:Nwindows
    
    TimeIntervalInTicks_=TimeIntervalInTicks;
    strt=Nsamples*(i-1)+1;
    StrtTick=TimeIntervalInTicks_*(i-1);
    EndTick = StrtTick+TimeIntervalInTicks_-1;
    if EndTick>maxTimeTicks
        TimeIntervalInTicks_=EndTick-StrtTick;
        EndTick=maxTimeTicks;
    end
    values=[];
    [~,values] = CEDS64ReadWaveF(fhand1,ch,TimeIntervalInTicks_,StrtTick,EndTick);
    ChannelsSig(strt:strt+length(values)-1,1)=values;
    disp([num2str(round(StrtTick/maxTimeTicks*1000)/10),'% complete'])
end
disp([num2str(100),'% complete'])
toc

N=length(ChannelsSig);
eval([desiredchannel '=ChannelsSig;']);
clear ChannelsSig

%% Spike time
[~,TimeDateOut] = CEDS64TimeDate(fhand1);
StarttimeSpike=sum(double(TimeDateOut).*[0.01 1 60 3600 0 0 0]);
if StarttimeSpike<43200, StarttimeSpike=StarttimeSpike+86400; end

dt=1/Fs;
Time = StarttimeSpike + (0:dt:(N-1)*dt)';

%% Plot

figure(1);
if flowtype==1
    plot(Time,Flow);
    ylabel('Flow');
    invertflow = 0;
    asymmetricPnasal=0;
    
    Flow(isnan(Flow))=0;
elseif flowtype==2
    plot(Time,Pnasal);
    ylabel('Pnasal');
    invertflow = 1;
    asymmetricPnasal=1;
    %was added to deal with asymmetric flow (Pnasal)
    %Q: is failure to sync due to the DC offset. set to zero to use fast
    %baseline search. answer: no, was due more to asymmetry.
    
    Pnasal(isnan(Pnasal))=0;
end



%% Find periods of noise or absent signal
if 1
    %defaults:
    global settings
    settings.plotfigure=0;
    settings.sqrt_scaling=0;
    settings.scalingexponent=0.5;
    if flowtype==1
        noisewav = FlowSignalToNoise(Time,Flow,1);
    elseif flowtype==2
        noisewav = FlowSignalToNoise(Time,Pnasal,1);
    end
end
%% Option: Estimate baseline flow (All options)
if 1%1flowtype==2&&asymmetricPnasal==0
    tempflow=(-2*invertflow+1)*eval(desiredchannel);
    
    excludenoiseperiods=0;
    N=length(Time);
    
    Twindow=[10 180 1200];
    TeTtotLims=[33 67]; % normal Ti/Ttot% min and max value
    plotfig=1;
    [XbaselineMTArs,IEestimate]=flowbaseline(tempflow,Time,Twindow,TeTtotLims,excludenoiseperiods,noisewav,plotfig);
    
    tol=0.1*1.96*std(tempflow-XbaselineMTArs);
    
    
    if 0
        tic
        [time,Vflow,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot,leak,...
            IEratio,VT,Vpeak,Vpeakmean,Apnea_B,Vflow_out,VTi,VTe,VTi_leakreverse,VTe_leakreverse,leaksig] = VEfromFlowX(time,tempflow-XbaselineMTArs);
        figure(2); clf(2);
        tempvol=cumsum(tempflow)*dt;
        ax(1)=subplot(3,1,1);
        plot(time,tempvol); hold('on');
        plot(time(BB_i_start),tempvol(BB_i_start),'k.');
        plot(time(BB_i_mid),tempvol(BB_i_mid),'r.');
        ax(2)=subplot(3,1,2);
        plot(time,tempflow); hold('on');
        plot(time,XbaselineMTArs); hold('on');
        plot(time,leaksig); hold('on');
        plot(time(BB_i_start),tempflow(BB_i_start),'k.');
        plot(time(BB_i_mid),tempflow(BB_i_mid),'r.');
        ax(3)=subplot(3,1,3);
        stairs(time(BB_i_start),VTi_leakreverse,'k'); hold('on');
        plot(time,0*time,'k:');
        stairs(time(BB_i_start),VTe_leakreverse,'r');
        stairs(time(BB_i_start),VTi,'k:');
        stairs(time(BB_i_start),VTe,'r:');
        linkaxes(ax,'x');
        Ti = (BB_i_mid-BB_i_start)*dt;
        TionTtot=median(Ti./Ttot');
        toc
    end
end
%% Run full breath detection all night to generate a transformed Pnasal signal (Flow2) with insp/exp more symmetrical
if 1
    if abs(log10(IEestimate))<log10(1.5) %3 fold error, 1 gives zero
        asymmetricPnasal=0;
    end
    
    if flowtype==2&&asymmetricPnasal
        clear I
        settings.sqrt_scaling=1;
        settings.scalingexponent=2/3;
        Flow2 = 0*Pnasal;
        Count = 0*Pnasal;
        N=length(Time);
        Foverlap=0.25;
        Twindow=180;
        Nlength = round(Twindow/dt);
        Nstep = round(Nlength*(1-Foverlap));
        M = floor((N-Nlength)/Nstep);
        for m=1:M
            try
                disp([num2str(round(m/M*100*10)/10),'% complete'])
                Li = 1+(m-1)*Nstep;
                Ri = Li+Nlength-1;
                if Ri>N, Ri=N; end
                Time_ = Time(Li:Ri);
                Flow_ = -Pnasal(Li:Ri);
                noisewav_ = noisewav(Li:Ri);
                if sum(noisewav_==3)/length(noisewav_)>0.1
                    continue
                end
                if 0
                    figure(90);
                    plot(Time_,Flow_);
                end
                [~,Flow2_,I.starti,I.midi,I.endi,~,~,~,~,~,~,VT,~,~,~,~,~]=VEfromFlow_sqrt_V16(Time_,Flow_);
                if 0
                    plot(Time_,[Flow_ Flow2_]);
                    hold('on')
                end
                Flow2(Li:Ri)=(Count(Li:Ri).*Flow2(Li:Ri)+Flow2_)./(Count(Li:Ri)+1);
                Count(Li:Ri)=Count(Li:Ri)+1;
                %pause(1)
                if 0
                    figure(90);
                    plot(Time,Flow2);
                end
            catch me
                disp('failed to run VEfromFlow analysis');
                disp(me.message);
            end
        end
    end
end
%% Find Nox Data -- Searching for a Nox signal within the Spike data folder
directory = [file(1:max(strfind(file,'\')))]

yy=num2str(TimeDateOut(7));
mm=num2str(TimeDateOut(6));
if length(mm)==1
    mm=['0' mm];
end
dd=num2str(TimeDateOut(5));
if length(dd)==1
    dd=['0' dd];
end
textstr = [yy mm dd]; %not used yet

x = dir(directory);
x(1:2)=[];
isnoxdir=zeros(length(x),1);
for m=1:length(x)
    if x(m).isdir==0
        continue
    end
    if length(x(m).name)==23
        if x(m).name(1)=='2'&&x(m).name(9)=='T'
            isnoxdir(m)=1;
        end
    end
end
x(isnoxdir==0)=[];

%folders = uipickfiles('FilterSpec',directory);


%% New: Import Nox Signals (only folder k=1 at present)
k=1;
noxmatfilepathfilename = [directory x(k).name '\Matlab-Files\NoxMatv1.mat']
if 1&&~(exist(noxmatfilepathfilename)==2)
    disp('No NoxMat file found, 2trying Cloud conversion tool');
    currentdir = pwd;
    try
        if 1
            cd([pwd '\NoxDataTransferTool\']);
        else
            cd('E:\NoxDataTransferTool\');
        end
        AConvertNoxMat(directory);
    catch me
        disp(me);
        disp('failed conversion tool');
    end
    cd(currentdir);
end

filehandle = matfile(noxmatfilepathfilename);
w = whos('-file',noxmatfilepathfilename);

ChannelsList={'Channel_7','Inductance_Thorax','Inductance_Abdomen','Audio_Volume',...
    'X_Axis','Y_Axis','Z_Axis'...
    };

for i=1:length(ChannelsList)
    foundamatch=0;
    try
        eval([ChannelsList{i} '=filehandle.' ChannelsList{i},';']);
        foundamatch=1;
    catch me
    end
end

Channels_Fs = NaN*ones(length(ChannelsList),1);
Channels_StartTime = NaN*ones(length(ChannelsList),1);
Channels_sigStartTime = [];
ChannelFound = zeros(length(ChannelsList),1);

%Signals = filehandle.Signals;
TimeInfo = filehandle.TimeInfo;

%get sampling rates
for i=1:length(ChannelsList)
    if exist(ChannelsList{i},'var')
        ChannelFound(i)=1;
        try
            %I=find(strcmp(channelnamestemp,ChannelsList{i}));
            %I2=find(strcmp(Signals.labelsNoSpace,ChannelTitlesTemp{I}));
            %Channels_Fs(i)=eval([ChannelsList{i} '.fs']);
            Channels_Fs(i,1)=eval([ChannelsList{i} '.fsTrue']);
            Channels_StartTime(i,1)=eval([ChannelsList{i} '.StartTime']);
            Channels_sigStartTime{i,1}=eval([ChannelsList{i} '.sigStartTime']);
        catch me
            disp(['missing sampling rate']);
        end
    else
        % display to command window, but not to PUPbeta gui
        disp(['strewth, no ', ChannelsList{i}, ' found']);
    end
end

%strip data from struct
for i=1:length(ChannelsList)
    if ChannelFound(i)==1
        temp = eval([ChannelsList{i} '.data']);
        eval([ChannelsList{i} '=temp(:);']);
    end
end

%% NoxPositionSignals
if 1 %save this for later at present
    Pos_AccMagnitude = sqrt(X_Axis.^2+Y_Axis.^2+Z_Axis.^2);
    %supine: Y=-1,X=0,Z=0; code:1/(pi/2)*atan2(0,1)
    %upright: Y=-0,X=0,Z=-1; code:1/(pi/2)*atan2(0,1)
    NoxPosition = 1/(pi/2)*atan2(X_Axis./Pos_AccMagnitude,-Y_Axis./Pos_AccMagnitude); %quarter turns
    %Incline = 1/(pi/2)*atan2(-Z_Axis./Pos_AccMagnitude,sqrt(Y_Axis.^2+X_Axis.^2)./Pos_AccMagnitude); %quarter turns
    NoxUpright = -Z_Axis./Pos_AccMagnitude;
    
    I = find(strcmp(ChannelsList,'X_Axis')==1)
    %Fs_temp = Channels_Fs(I);
    %ChannelFound(I+2)=0; %assumes X,Y,Zaxis are together
    ChannelsList=[ChannelsList,'NoxPosition'];
    ChannelsList=[ChannelsList,'NoxUpright'];
    Channels_Fs=[Channels_Fs;Channels_Fs(I);Channels_Fs(I)];
    ChannelFound = [ChannelFound;1;1];
end

%%
StartTime = TimeInfo.StartTime;
StopTime = TimeInfo.EndTime;

if StartTime<43200
    StartTime=StartTime+86400; 
    StopTime=StopTime+86400; 
end

%% Start loop



%
%
%             AbNoxRIP.values = 0*Time;
%             ThNoxRIP.values = 0*Time;
%             %kNoxRIP.values = 0*TimeNoxHz;
%             kFlow.values = 0*Time;
%             alphaFlow.values = 0*Time;
%             NoxTimeShift.values = 0*Time;


%%
k=1


%% Replace zeros in RIP with a constant average of last known / next known values
dtNox=1/Channels_Fs(2);
RIPart = 0*Inductance_Thorax;

artranges = 1.0000e-05*[0.001 0.6];
minlengthi = 30/dtNox;

Itemp = Inductance_Thorax<=artranges(1)|Inductance_Thorax>=artranges(2)|Inductance_Abdomen<=artranges(1)|Inductance_Abdomen>=artranges(2);
li = find(diff(Itemp)==1)+1;
ri = find(diff(Itemp)==-1)+1;
[li,ri] = TidyStartEndEventList(li,ri,length(Inductance_Abdomen));

if 1 %remove short artifacts
    lengths = [li(2:end) - ri(1:end-1)]; %time in samples until next artifact
    Itemp = lengths<minlengthi;
    Itemp_ = [0;Itemp];
    li(Itemp_==1)=[];
    ri(Itemp ==1)=[];
    lengths = [li(2:end) - ri(1:end-1)]; %check only
end

for i=1:size(li,1)
    RIPart(li(i):ri(i))=1;
end

ABupper=5*(prctile(Inductance_Abdomen(RIPart==0),95)-median(Inductance_Abdomen(RIPart==0)))+median(Inductance_Abdomen(RIPart==0));
THupper=5*(prctile(Inductance_Thorax(RIPart==0),95)-median(Inductance_Thorax(RIPart==0)))+median(Inductance_Thorax(RIPart==0));

Itemp = Inductance_Thorax<=artranges(1)|Inductance_Thorax>=min([artranges(2) THupper])|Inductance_Abdomen<=artranges(1)|Inductance_Abdomen>=min([artranges(2) ABupper]);
li = find(diff(Itemp)==1)+1;
ri = find(diff(Itemp)==-1)+1;
[li,ri] = TidyStartEndEventList(li,ri,length(Inductance_Abdomen));
if 1 %remove short artifacts
    lengths = [li(2:end) - ri(1:end-1)]; %time in samples until next artifact
    Itemp = lengths<minlengthi;
    Itemp_ = [0;Itemp];
    li(Itemp_==1)=[];
    ri(Itemp ==1)=[];
    lengths = [li(2:end) - ri(1:end-1)]; %check only
end

for i=1:size(li,1)
    i_good=[li(i)-1 ri(i)+1];
    i_good(i_good<1|i_good>length(Inductance_Abdomen))=[];
    Inductance_Thorax(li(i):ri(i))=mean(Inductance_Thorax(i_good));
    Inductance_Abdomen(li(i):ri(i))=mean(Inductance_Abdomen(i_good));
    RIPart(li(i):ri(i))=1;
end

%% Main
Duration = StopTime-StartTime;

N2 = length(Inductance_Thorax);

TimeNoxRIP1 = (StartTime:dtNox:(StartTime+dtNox*(N2-1)))';



%% Import RIP artifact list

directorytemp = PathName;
textfilename=[directorytemp 'RIPart_timeofday.txt'];
artifactsignal=0*Inductance_Thorax;
if exist(textfilename,'file')==2
    ['found RIP artifact file and removed artifact']
    %Position.values_original=Position.values;
    [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
    for i=1:length(col1)
        dttemp=dtNox;
        try
            lefti=round((col1(i)-TimeNoxRIP1(1))/dttemp +1);
        catch me
            lefti=round((col1(i)-0)/dttemp +1);
        end
        if lefti<1, lefti=1; end
        try
            righti=round((col2(i)-TimeNoxRIP1(1))/dttemp +1);
        catch me
            righti=round((col2(i)-0)/dttemp +1);
        end
        if righti>length(TimeNoxRIP1), righti=length(TimeNoxRIP1); end
        ii=[lefti-1 righti+1];
        ii(ii<1|ii>length(TimeNoxRIP1))=[];
        thval = nanmedian(Inductance_Thorax(ii));
        abval = nanmedian(Inductance_Abdomen(ii));
        Inductance_Thorax(lefti:righti)=thval;
        Inductance_Abdomen(lefti:righti)=abval;
        artifactsignal(lefti:righti)=1;
    end
end

%% If flow copy is available in Nox

noisewavRIP = interp1(Time,noisewav,TimeNoxRIP1,'nearest','extrap');
   
   
    if ~asymmetricPnasal
        FlowSpikeSync2 = interp1(Time,(-2*invertflow+1)*eval(desiredchannel),TimeNoxRIP1,'linear','extrap');
    else
        FlowSpikeSync2 = interp1(Time,Flow2,TimeNoxRIP1,'linear','extrap');
    end
    
    if asymmetricPnasal==0&&flowtype==2
        XbaselineMTArs_NoxSync=interp1(Time,XbaselineMTArs,TimeNoxRIP1,'linear','extrap');
        FlowSpikeSync2=FlowSpikeSync2-XbaselineMTArs_NoxSync;
    end
    
    noisewav_NoxSync=interp1(Time,noisewav,TimeNoxRIP1,'nearest','extrap');
    
    FlowSpikeSync2(noisewav_NoxSync==3)=0;
    
    %%
if exist('Channel_7')
    TimeNoxFlow = (StartTime:1/Channels_Fs(1):(StartTime+1/Channels_Fs(1)*(length(Channel_7)-1)))';
    Channel_7_25Hz = interp1(TimeNoxFlow,Channel_7,TimeNoxRIP1,'linear','extrap');
    
    figure(3);
    ax3(1)=subplot(2,1,1); 
        plot(TimeNoxFlow,Channel_7); 
        ylabel('Nox Channel 7');
    ax3(2)=subplot(2,1,2); 
        plot(Time,eval(desiredchannel)); 
        ylabel(['Spike ' desiredchannel]);
    linkaxes(ax3,'x');
    
    %% Downsample Flow to 25 Hz
    
    FlowSpikeSync = interp1(Time,eval(desiredchannel),TimeNoxRIP1,'linear','extrap');
    
%     flipreffactor=1; %for when Pmask was input as flow, or flow was upside-down when input to Nox C1.
%     if exist([directory 'flipref.txt'],'file')==2
%         flipreffactor=-1;
%     end
    
    testflipref=1;
    %% Run crosscorr function to calculate lags between near-identical signals
    [LagaverageExact,flipreffactor]=SyncFlowExact(FlowSpikeSync,Channel_7_25Hz,TimeNoxRIP1,0.9,testflipref);
    
    %note Channel_7_25Hz not used after this
    %% calibrate RIP
    TimeNoxRIP2=TimeNoxRIP1+LagaverageExact; %breakpoint
    Inductance_Thorax_ = interp1(TimeNoxRIP2,Inductance_Thorax,TimeNoxRIP1,'linear','extrap');
    Inductance_Abdomen_ = interp1(TimeNoxRIP2,Inductance_Abdomen,TimeNoxRIP1,'linear','extrap');
    
    [kNox,alphaNox]=calibrateRIPusingFlow(Inductance_Thorax_,Inductance_Abdomen_,TimeNoxRIP1,FlowSpikeSync2,noisewavRIP,flowtype,asymmetricPnasal,dtNox,settings);
    
    %% 
    figure(12)
    ax12(1)=subplot(4,1,1); plot(TimeNoxRIP1,FlowSpikeSync2)
    ax12(2)=subplot(4,1,2); plot(TimeNoxRIP1,Inductance_Thorax_)
    ax12(3)=subplot(4,1,3); plot(TimeNoxRIP1,Inductance_Abdomen_)
    linkaxes(ax12,'x')
    
    LagaverageTotal = LagaverageExact;
else
    
    %% Make Vol Signal from pneumotach or Pnasal needed for RIP sync   
    %filter
    filter_HFcutoff_butter0 = 1/30; filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dtNox/2),'high');
    
    VolNoxSync = filtfilt(B_butter0,A_butter0,cumsum(FlowSpikeSync2)*dtNox);
    tempupper = prctile(VolNoxSync,99.75); templower = prctile(VolNoxSync,0.25);
    VolNoxSync(VolNoxSync>tempupper)=tempupper;  VolNoxSync(VolNoxSync<templower)=templower;   
    
    %% Make RIP vol signal - first assuming k=1
    RIPvol = filtfilt(B_butter0,A_butter0,(Inductance_Abdomen + Inductance_Thorax));
    tempupper = prctile(RIPvol,99.75); templower = prctile(RIPvol,0.25);
    RIPvol(RIPvol>tempupper)=tempupper;  RIPvol(RIPvol<templower)=templower;
    
    %%
    if flowtype==1
        Rthres=0.5;
    else
        Rthres=0.25;1
        
    end
    Lagaverage1=SyncFlowExact(VolNoxSync,RIPvol,TimeNoxRIP1,Rthres);
    TimeNoxRIP2=TimeNoxRIP1+Lagaverage1; %breakpoint
    Inductance_Thorax_ = interp1(TimeNoxRIP2,Inductance_Thorax,TimeNoxRIP1,'linear','extrap');
    Inductance_Abdomen_ = interp1(TimeNoxRIP2,Inductance_Abdomen,TimeNoxRIP1,'linear','extrap');
    
    Fisnan = sum(isnan(Lagaverage1))/length(Lagaverage1)
    
    [kNox,alphaNox]=calibrateRIPusingFlow(Inductance_Thorax_,Inductance_Abdomen_,TimeNoxRIP1,FlowSpikeSync2,noisewavRIP,flowtype,asymmetricPnasal,dtNox,settings);
    
    %% recalculate lag using alpha K
    RIPvol = filtfilt(B_butter0,A_butter0,(Inductance_Abdomen_.*kNox.^0.5 + Inductance_Thorax_./(kNox.^0.5)));
    tempupper = prctile(RIPvol,99.75); templower = prctile(RIPvol,0.25);
    RIPvol(RIPvol>tempupper)=tempupper;  RIPvol(RIPvol<templower)=templower;
    
    Lagaverage2=SyncFlowExact(VolNoxSync,RIPvol,TimeNoxRIP1,Rthres);
    
    Fisnan = sum(isnan(Lagaverage1))/length(Lagaverage1)
    
    LagaverageTotal=Lagaverage2+Lagaverage1;
end

%% Make lag-corrected versions of all Nox signals

resampleto125=1;
ChannelsList = [ChannelsList,'kNox','alphaNox']
Channels_Fs = [Channels_Fs;Channels_Fs(2);Channels_Fs(2)]
ChannelFound = [ChannelFound;1;1];

for j=1:length(ChannelsList)
    if ChannelFound(j)==0
        continue
    end
    NoxTime = (StartTime:1/Channels_Fs(j):(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
    NoxTimeCorrected = NoxTime + interp1(TimeNoxRIP1,LagaverageTotal,NoxTime,'linear','extrap');
    temp = interp1(NoxTimeCorrected,eval(ChannelsList{j}),NoxTime,'linear','extrap');
    eval([ChannelsList{j} '_Corrected=temp;']);
    if resampleto125
        %NoxTime2 = (StartTime:1/125:(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
        temp = interp1(NoxTime,eval([ChannelsList{j} '_Corrected']),Time,'linear','extrap');
        eval([ChannelsList{j} '_125Hz=temp;']);
    end
end

NoxPosition_125Hz(isnan(NoxPosition_125Hz))=0;

%%
if 0
   figure(2); clf(2); 
   plot(TimeNoxRIP1,[LagaverageExact,LagaverageTotal]);
end

%%
figure(2); clf(2); set(gcf,'color',[1 1 1]);

Nplots=6;

ax1(1)=subplot(Nplots,1,1); plot(Time,eval(desiredchannel)); set(gca,'xtick',[]); %,TimeNoxRIP1,FlowSpikeSync

if exist('Channel_7_125Hz')
ax1(2)=subplot(Nplots,1,2);
j=1;
NoxTime = (StartTime:1/Channels_Fs(j):(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
plot(NoxTime,Channel_7_Corrected); set(gca,'xtick',[]);
hold('on');
plot(Time,Channel_7_125Hz); set(gca,'xtick',[]);
end

ax1(3)=subplot(Nplots,1,3);
j=2;
NoxTime = (StartTime:1/Channels_Fs(j):(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
plot(NoxTime,Inductance_Thorax_Corrected); set(gca,'xtick',[]);
hold('on');
plot(Time,Inductance_Thorax_125Hz); set(gca,'xtick',[]);


ax1(4)=subplot(Nplots,1,4);
plot(NoxTime,Inductance_Abdomen_Corrected); set(gca,'xtick',[]);
hold('on');
plot(Time,Inductance_Abdomen_125Hz); set(gca,'xtick',[]);

ax1(5)=subplot(Nplots,1,5);
j=5;
NoxTime = (StartTime:1/Channels_Fs(j):(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
plot(NoxTime,NoxPosition_Corrected); set(gca,'xtick',[]);
hold('on');
plot(Time,NoxPosition_125Hz);

ax1(6)=subplot(Nplots,1,6);
j=4;
NoxTime = (StartTime:1/Channels_Fs(j):(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
plot(NoxTime,Audio_Volume_Corrected); set(gca,'xtick',[]);
hold('on');
plot(Time,Audio_Volume_125Hz);

linkaxes(ax1,'x');



%% Export Spike Data
if 1
    retest=1;
    if retest
        cedpath = [cd '\CEDS64ML2017'];
        addpath(cedpath); % so CEDS64LoadLib.m is found
        CEDS64LoadLib(cedpath); % load ceds64int.dll
        %fhand1 = CEDS64Open(file);
    end
    
    Fs=125;
    ChannelsToSave = {'Inductance_Abdomen_125Hz','Inductance_Thorax_125Hz', ...
        'NoxPosition_125Hz','Audio_Volume_125Hz','kNox_125Hz','alphaNox_125Hz','Channel_7_125Hz',desiredchannel}';
    ChannelsToSaveNames = {'AbNoxRIP','ThNoxRIP', ...
        'NoxPos','NoxAudio','kFlow','alphaFlow',['Nox' desiredchannel],desiredchannel}';
    
    writedirect=0; %not working at present, also risky
    if ~writedirect
        % Create a new Blank .smr file
        file2 = [file(1:max(strfind(file,'.'))-1) '_RIP.smrx'];
        fhand2 = CEDS64Create(file2,32,2);
        tbase = CEDS64TimeBase(fhand1);%CEDS64TimeBase( fhand1 );%%%obtain time divisions from the original file
        CEDS64TimeBase(fhand2,tbase); %% set the time division for new file
    else
        fhand2 = CEDS64Open(file);
        tbase = CEDS64TimeBase(fhand1);%CEDS64TimeBase( fhand1 );%%%obtain time divisions from the original file
    end
    %hardcoded based on other example recordings, works fine
    
    for i=1:length(ChannelsToSave)
        if ~exist(ChannelsToSave{i})
            continue
        end
        if writedirect
            ch=SpikeBlank(i);
        else
            ch=i;
        end
        iOk(i,1)=CEDS64SetWaveChan(fhand2,ch,1/(Fs*tbase),9); %%% set channel 1-- 9 is for channel type "waveform"
        % CEDS64WriteWave(fhand2,ch,Channel_7_125Hz,0);
        eval(['CEDS64WriteWave(fhand2,ch,' ChannelsToSave{i} ',0);']);
        iOk(i,2) = CEDS64ChanTitle( fhand2, ch,ChannelsToSaveNames{i});
    end
    disp(iOk==0)
    CEDS64Close(fhand2);
    
end
disp('Completed Spike export');
% Close Spike
CEDS64CloseAll(); % close all the files
unloadlibrary ceds64int; % unload ceds64int.dll


%% Export MAT Data
%
% save([file(1:max(strfind(file,'.'))-1) '_RIP.mat'],'ThNoxRIP','AbNoxRIP','alphaFlow','kFlow','StarttimeSpike','NoxTimeShift','-v7.3');
% disp('Completed Matlab export');
%


ChannelsList_ = ChannelsList;
ChannelsList_(ChannelFound==0)=[];
for i=1:length(ChannelsList)
    try
        eval([ChannelsList{i} '=' ChannelsList{i} '_125Hz;']);
    catch me
    end
end

TimeNoxRIP=TimeNoxRIP1;
TimeLag=LagaverageTotal;
VarsToSaveMat = [ChannelsList_,'StartTime','TimeNoxRIP','TimeLag'];


filenameMat = [file(1:max(strfind(file,'.'))-1) '_RIP.mat'];
save(filenameMat,VarsToSaveMat{:},'-v7.3')

%Code to sync with Spike time in future:
% NoxTime = (StartTime:1/Channels_Fs(j):(StartTime+1/Channels_Fs(j)*(length(eval(ChannelsList{j}))-1)))';
% NoxTimeCorrected = NoxTime + interp1(TimeNoxRIP,TimeLag,NoxTime,'linear','extrap');
% temp = interp1(NoxTimeCorrected,eval(ChannelsList{j}),NoxTime,'linear','extrap');
% eval([ChannelsList{j} '_Corrected=temp;']);


if 0
    save workspacetemp
end






