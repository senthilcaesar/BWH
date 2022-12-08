clear all


%Using new Nox reader, and broader sync range (15 s) since 7/7/2018

%Before 9/19/17, K was selected based on an average, but this should really have
%been a average of log values. Corrected on 9/19/17. Rerun those prior to
%this date for optimal K calculation.

%User guide:
%0. Run this file (F5)
%1. User selects spike file.
%2. Nox folder(s) inside the spike file directory is automatically found.
%3. New Spike and .mat files with RIP data are output to the same folder.
%4. Later merge spike data, see "HowTo.txt"

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
[FileName,PathName] = uigetfile([default],'Select Spike file');
file=[PathName FileName]

fid = fopen(a_fnamelastpath,'wt');
fprintf(fid, '%s', [PathName '*.*']);
fclose(fid);

%% Load spike library

addpath(genpath('G:\Partners Healthcare Dropbox\SATP Group\NoxSignalsToSpike\NoxSignalsToSpike\'));
cedpath = [cd '\CEDS64ML2017'];

%%
addpath( cedpath ); % so CEDS64LoadLib.m is found
CEDS64LoadLib(cedpath); % load ceds64int.dll


%% channels to load
clear channelnameoptions
channelnameoptions.Flow={'Flow','Vflow'}; %'Vflow'
desiredchannel = char(fieldnames(channelnameoptions));

%% Open File

fhand1 = CEDS64Open(file);
if fhand1==-21 || fhand1==-1
    disp('Spike file is already open, close it and then try again');
end

if (fhand1 <= 0); unloadlibrary ceds64int;
    'failed at opening Spike file, CEDS64Open'
    return;
end

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
% 
% figure(1);
% if flowtype==1
%     plot(Time,Flow);
%     ylabel('Flow');
%     invertflow = 0;
%     asymmetricPnasal=0;
%     
%     Flow(isnan(Flow))=0;
% elseif flowtype==2
%     plot(Time,Pnasal);
%     ylabel('Pnasal');
%     invertflow = 1;
%     asymmetricPnasal=1;
%     %was added to deal with asymmetric flow (Pnasal)
%     %Q: is failure to sync due to the DC offset. set to zero to use fast
%     %baseline search. answer: no, was due more to asymmetry.
%     
%     Pnasal(isnan(Pnasal))=0;
% end
%% Code to edit signal

Flow = -Flow;
%
%Flow((Time-Time(1))>15150)=-Flow((Time-Time(1))>15150);


%% Export Spike Data
if 1
    retest=1;
    if retest
        cedpath = [cd '\CEDS64ML2017'];
        addpath(cedpath); % so CEDS64LoadLib.m is found
        CEDS64LoadLib(cedpath); % load ceds64int.dll
        %fhand1 = CEDS64Open(file);
    end
        
    ChannelsToSave = {desiredchannel}';
    ChannelsToSaveNames = {desiredchannel}';
    FsList = [Fs]';
    
    writedirect=1; %working (just needed to open in write mode), but needs testing as it is inherently risky
    if ~writedirect
        % Create a new Blank .smr file
        file2 = [file(1:max(strfind(file,'.'))-1) '_edits.smrx'];
        fhand2 = CEDS64Create(file2,32,2);
        tbase = CEDS64TimeBase(fhand1);%CEDS64TimeBase( fhand1 );%%%obtain time divisions from the original file
        CEDS64TimeBase(fhand2,tbase); %% set the time division for new file
%         if 1 %note if we simply close and open the file here we then cannot write data to it
%             CEDS64Close(fhand2);
%             fhand2 = CEDS64Open(file2,0);
%         end
    else
        fhand2 = CEDS64Open(file,0);
        tbase = CEDS64TimeBase(fhand2);%CEDS64TimeBase( fhand1 );%%%obtain time divisions from the original file
    end
    %hardcoded based on other example recordings, works fine
    
    clear iOk
    for i=1:length(ChannelsToSave)
        if ~exist(ChannelsToSave{i})
            continue
        end
        if writedirect %working!
            ch=SpikeBlank(i);
        else
            ch=i;
        end
        iOk(i,1)=CEDS64SetWaveChan(fhand2,ch,1/(FsList(i)*tbase),9); %%% set channel 1-- 9 is for channel type 
        %iOk(i,3) = CEDS64WriteWave(fhand2,ch,eval(ChannelsToSave{i}),0);
        eval(['CEDS64WriteWave(fhand2,ch,' ChannelsToSave{i} ',0);']);
        iOk(i,2) = CEDS64ChanTitle( fhand2, ch,ChannelsToSaveNames{i});
    end
    disp(iOk==0)
    CEDS64Close(fhand2);
    
end
disp('Completed Spike export');
%% Close Spike
CEDS64CloseAll(); % close all the files
unloadlibrary ceds64int; % unload ceds64int.dll




