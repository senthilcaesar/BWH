function outdata = NDFreadN(noxdir,SaveData)
%%Signals only, no events data
if ~exist('noxdir','var') || isempty(noxdir)
    noxdir = uigetdir();
end
    
if 0
    LoadAndCompare=1;
    ploton=1;
else
    LoadAndCompare=0;
    ploton=0;
end
outdata = [];

devmode=0;
if ~exist('SaveData','var')
SaveData=0;
end
LoadSavedData=0;

%noxdir = 'E:\Dropbox (Partners HealthCare)\PAtO Studies\ScoredES\DPW2131PatoPhysio\20210726T223651 - fe7d3\'

if ~exist('noxdir','var') || isempty(noxdir)
    noxdir = uigetdir('Select Actual Nox Folder');
end

if noxdir(end)~=filesep
   noxdir = [noxdir filesep];
end

dir_ = dir(noxdir);
signames = struct2table(dir_);
signames.isndf = contains(signames.name,'.ndf');

%Audio wav file is not a standard .ndf signal file
Iaudio = find(signames.name == "audio.ndf");
signames.isndf(Iaudio)=0;

%resp rate is breaking, need to debug 
Iaudio = find(signames.name == "resp rate.ndf");
signames.isndf(Iaudio)=0;

%%
%siglist= {'inductance abdomen.ndf'}
%siglist= {'inductance thorax.ndf'}
siglist = signames.name(signames.isndf)';
%siglist= {'channel 7.ndf'}
%siglist= {'voltage (battery).ndf'}
%siglist= {'audio volume db.ndf'}
%siglist= {'audio volume.ndf'}
%siglist= {'activity.ndf'}
%%
FsFactor = 2.083224969283037e-05; %applied to most, not ExcludeFromFsFixList1 or 2; only used when fsTrue is an integer 
%FsFactor2 = -3.694864854830726e-05; 
FsFactor2 = (99.996380101692400/100)-1; %used in files with NDFread fsTrue==100

if any(siglist=="audio volume.ndf")
    audioinfo = NDFread1([noxdir 'audio volume.ndf']);
    audioFs = audioinfo.fsTrue;
    if audioFs<99.997 %
        FsFactor2 = 99.996214041937170/99.996040564196800-1; %used in files with NDFread fsTrue==99.996040564196800;
    end
end
        
%         ExcludeFromFsFixList2 = {};  %older data
%         %NoxMAT actual Fs = 99.996214041937170
%         %NDFread fsTrue shows up as 99.996040564196800
%         %99.996214041937170/99.996040564196800-1
%     else
        ExcludeFromFsFixList2 = {'Audio_Volume','Voltage__battery_','Voltage__bluetooth_','Voltage__core_'};
%     end
% end
%         
% if ~any(siglist=="audio volume db.ndf")
%     ExcludeFromFsFixList2 = {};  %older data, recognizable from absence of an "audio volume db" signal
% else
%     ExcludeFromFsFixList2 = {'Audio_Volume','Voltage__battery_','Voltage__bluetooth_','Voltage__core_'};
%     %apply FsFactor2
% end
%ExcludeFromFsFixList1 = {'Heart_Rate','K','RIP_Phase'};

%Calculated channels with whole number sample rates are typically correct.
ExcludeFromFsFixList1 = {'Heart_Rate','K','RIP_Phase','Activity','Audio_Volume_dB'}; %updated 8/10/23, Activity is actually exactly 20 Hz not 20.0004---

%Audio_Volume Fs = 99.996276081160800 corrected to 99.992581353899740 in newer data from 20211104T213422 without audio volume db present; works well using FsFactor above
% correction from 99.996040564196800 causes problems when corrected to 99.992345845637770 in 20210712T215000 - 982c9% 
% exact Fs for audio volume from NoxDataTransferTool 99.996380101692400 8/8/2023
%Voltage signals are still not quite correct
%%
for i=1:length(siglist)
    filename =([noxdir siglist{i}]);
    try
    info = NDFread1(filename);
    
    catch
	disp(['error: failed ' siglist{i}]);
    continue
    end
    %Fs fix
    if 1 && ~any(ExcludeFromFsFixList1==string(info.Label)) && ~any(ExcludeFromFsFixList2==string(info.Label)) %except Herat Rate, K, RIPphase, Voltage_Battery, ... also Audio_Volume needs a separate fix
        if round(info.fsTrue)==info.fsTrue
            info.fsTrue = info.fsTrue*(1+FsFactor);
        else
            info.fsTrue = info.fsTrue;
        end
    end
    if any(ExcludeFromFsFixList2==string(info.Label))
        info.fsTrue = info.fsTrue*(1+FsFactor2); 
    end
    %interestingly the Fs from audio volume db.ndf seems correct, but this signal isn't always present
%     if siglist{i}=="audio volume.ndf" && any(signames.name(signames.isndf)'=="audio volume db.ndf")
%        info2 = NDFread1([noxdir 'audio volume db.ndf']);
%        info.fsTrue=info2.fsTrue;
%     end
    
    StartTime = info.StartTime;
    Fs=(info.fsTrue);
    

    
    %handle gaps in sessions here
    if info.Nsessions>1
        
        clear sessionstart sessionend
        for j=1:info.Nsessions
            sessionstart(j,1) = mod(datenum(info.SessionStart(j),'HH:MM:SS.FFF'),1)*86400; % SO 20230328, get start and end of EDF, needed below
            sessionend(j,1) = sessionstart(j)+(info.SessionLength(j)-1)*(1/Fs); %last sample
        end
        Inextday = sessionstart<86400/2;
        sessionend(Inextday)=sessionend(Inextday)+86400;
        sessionstart(Inextday)=sessionstart(Inextday)+86400;
        gapsize = [(sessionstart(2:end) - sessionend(1:end-1)) - 1/Fs; 0];
        
        data = [];
        sessionstarts = 1+[0 cumsum(info.SessionLength(1:end-1))];
        for j=1:info.Nsessions
            sessiondata = info.data(sessionstarts(j)-1+[1:info.SessionLength(j)]);
            gapdata = zeros(round(gapsize(j)*Fs),1);
            data = [data;sessiondata;gapdata];
        end
        info.data = data;
    end
    
    N_timeXHz=length(info.data);
    Time=(StartTime:(1/Fs):StartTime+(N_timeXHz-1)*(1/Fs))'; % This is the time vector associated with the _XHz Flow data.
    info.EndTime = Time(end);
    
    if ploton
        figure(101);clf(101);
        plot(Time,info.data);
        hold on;
    end
    outdata = setfield(outdata,info.Label,info);
    %
    if LoadAndCompare && ploton
        getfilename = dir([noxdir '\Matlab-Files\']);
        getfilename(1:2)=[];
        temp=getfilename.name;
        h=load([noxdir '\Matlab-Files\' temp],info.Label);
        %h=load([noxdir '\Matlab-Files\NoxMat_debug.mat'],info.Label);
        if isempty(fieldnames(h))
            continue
        end
        d=getfield(h,info.Label);
        
        StartTime3=d.StartTime;
        Fs2=d.fsTrue;
        N_timeXHz2=length(d.data);
        Time3=(StartTime3:(1/Fs2):StartTime3+(N_timeXHz2-1)*(1/Fs2))'; % This is the time vector associated with the _XHz Flow data.
        
        plot(Time3,d.data - 2*nanstd(d.data));
        xlim([prctile(Time3,95)+[0 60]]);
        ylabel(info.Label,'interpreter','none');
        pause()
    end
    %%
end

%% Get the EndTimes from each signal
if 0
varlist = {'fsTrue','StartTime','EndTime','N'};
siginfoT = DealNaNIntoTable(length(siglist),varlist);
siginfoT.Label=strings(length(siglist),1);
fieldnames_=fieldnames(outdata);
for i=1:length(fieldnames_)
    temp=getfield(outdata,fieldnames_{i});
    siginfoT.Label(i)=string(getfield(temp,'Label'));
    siginfoT.StartTime(i)=getfield(temp,'StartTime');
    siginfoT.EndTime(i)=getfield(temp,'EndTime');
    siginfoT.fsTrue(i)=getfield(temp,'fsTrue');
    siginfoT.N(i)=getfield(temp,'fsTrue');
    siginfoT.N(i)=length(getfield(temp,'data'));
    
    if round(siginfoT.fsTrue(i))==siginfoT.fsTrue(i)
        siginfoT.fsBetter(i) = siginfoT.fsTrue(i)*(1+FsFactor);
    else
        siginfoT.fsBetter(i) = siginfoT.fsTrue(i);
    end
end
i = find(siginfoT.Label=="cRIP_Flow")
siginfoT.EndTimeDiff = siginfoT.EndTime-siginfoT.EndTime(i);

%siginfoT.fsBetter= 1./((siginfoT.EndTime(i)-siginfoT.StartTime)./(siginfoT.N-1));

for i=1:length(fieldnames_)
    sig=getfield(outdata,fieldnames_{i});
    sig.fsBetter = siginfoT.fsBetter(i);
    outdata=setfield(outdata,fieldnames_{i},sig);
end

end
%%
%temp = 2.083224969283037e-05;
%outdata.Inductance_Abdomen.fsTrue = outdata.Inductance_Abdomen.fsTrue*(1+temp);
%outdata.Inductance_Thorax.fsTrue = outdata.Inductance_Thorax.fsTrue*(1+temp);

%% Hack Fs for Audio Volume, removed and applied above using FsFactor2
% try
%     %this sampling rate is incorrect in the original Nox file, should be 99.9964 not 100 Hz
%     outdata.Audio_Volume.fsTrue = outdata.Audio_Volume_dB.fsTrue;
% end

%% Load Events

    try
        mksqlite('open',[noxdir 'Data.ndb']);
        % Query the database
        tableslist = struct2table(mksqlite(['show tables']));
        EventsT = struct2table(mksqlite(['SELECT * FROM scoring_marker']));
        mksqlite('close');
        
        EventsT = sortrows(EventsT,'starts_at');
        date=floor(min(double(EventsT.starts_at))/10000000/86400);
        EventsT.StartTime = 86400*(double(EventsT.starts_at)/10000000/86400-date);
        if min(EventsT.StartTime)<86400/2
            date=date-1;
            EventsT.StartTime = EventsT.StartTime+86400;
        end
        EventsT.Duration = double(EventsT.ends_at - EventsT.starts_at)/10000000;
        EventsT.StartTimeHHMM = datetime(date + 367 + EventsT.StartTime/86400,'ConvertFrom','datenum');
    end



%% savedata
if SaveData
    disp('Finished importing, started saving, ~15 sec');
    tic
    if ~exist([noxdir 'Matlab-Files\'],'dir')
        mkdir([noxdir 'Matlab-Files\']);
    end
    save([noxdir '\Matlab-Files\' 'NdfMat'],'-struct','outdata','-v7.3');
    toc
    disp('Finished saving');
end

if LoadSavedData
    x=load([noxdir 'NdfMat']);
end

%%
if devmode %development checking mode only
figure(101);clf(101);

plotList = {'Channel_7','C3'} %,'cRIP_Flow','C3' ,'Inductance_Abdomen','Inductance_Thorax'
for i=1:length(plotList)
ax(i) = subplot(length(plotList),1,i);
sig = getfield(outdata,plotList{i});
StartTime = sig.StartTime;
Fs=(sig.fsTrue);
%Fs=(sig.fsBetter);
%temp = (Fs2-Fs)/Fs;
%Fs=Fs2;
N_timeXHz=length(sig.data);
Time=(StartTime:(1/Fs):StartTime+(N_timeXHz-1)*(1/Fs))'; % This is the time vector associated with the _XHz Flow data.
plot(Time,sig.data);
hold on;
TidyPlot()

    if 1
            h=load([noxdir '\Matlab-Files\NoxMat_debug.mat'],plotList{i});
            if isempty(fieldnames(h))
                continue
            end
            d=getfield(h,plotList{i});

            StartTime3=d.StartTime;
            Fs2=d.fsTrue;
            N_timeXHz2=length(d.data);
            Time3=(StartTime3:(1/Fs2):StartTime3+(N_timeXHz2-1)*(1/Fs2))'; % This is the time vector associated with the _XHz Flow data.

            plot(Time3,d.data);


    end

    


end

linkaxes(ax,'x')
end