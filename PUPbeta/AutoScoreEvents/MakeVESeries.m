%% Start
clear all;
clc;
addpath(genpath(pwd));
addpath(genpath('E:\PUPbeta_git\PUPbeta20190629'));

%% Choose dataset

global settings
settings.CurrentCodeVersion = 'PUPbeta20190629';
mydir  = mfilename('fullpath'); %only works under "run" not "evaluate line"
idcs   = strfind(mydir,filesep);
settings.workdir = mydir(1:idcs(end-1));
settings.codedir = ['E:\PUPbeta_git\' settings.CurrentCodeVersion '\'];

addpath(genpath('E:\PUPbeta_git\PUPbeta20190629\'));

savename = 'MrOS_Oct2019';
settings.AMasterdir = 'D:\MrOS\PUPStart\';
settings.AMasterdir = 'D:\MrOS\PUPStart\';
AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
idcs   = strfind(settings.AMasterdir,filesep);
settings.workdir = settings.AMasterdir(1:idcs(end-1));
[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
Npatients = size(patients,1);
clear filedir
for i=1:Npatients
    filedir{i,1} = ['D:\MrOS\Analyzed\' savename '_' num2str(i)];
end
for i=1:Npatients
    filedirC{i,1} = [settings.workdir 'Converted\' patients{i,1}];
end
settings.savename='MrOs2019';
noscoredarinwake=0; %0 is have scoring of arousals in wake


% savename = 'RICCADSA';
% settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\PUPstart\';
% settings.AMasterdir = 'G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\PUPstart\';
% AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
% idcs   = strfind(settings.AMasterdir,filesep);
% settings.workdir = settings.AMasterdir(1:idcs(end-1));
% [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
% Npatients = size(patients,1);
% clear filedir
% for i=1:Npatients
%     filedir{i,1} = ['G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Analyzed\' savename '_' num2str(i)];
% end
% for i=1:Npatients
%     filedirC{i,1} = [settings.workdir 'Converted\' patients{i,1}];
% end
% MultiScorer = xlsread(AMasterSpreadsheet,1,'R4:R10003');
% filedir(MultiScorer>0)=[];
% filedirC(MultiScorer>0)=[];
% Npatients = size(filedir,1)
% settings.savename='RICCADSA2019';
% noscoredarinwake=0; %0 is have scoring of arousals in wake

addpath(genpath(settings.workdir));
addpath(genpath(settings.codedir));
addpath(genpath(pwd));

%% Select subject and load Converted and Analyzed data
n=2

loadpath=[filedir{n}];

loadpathC=[filedirC{n}];
load(loadpathC);

load(loadpath,'BreathDataTable');

%%
[BreathDataTable1,BreathDataTable2]=GetNonOvlappedVE(BreathDataTable);

%find indices at which data starts new windows
temp = find(diff(BreathDataTable1.Time0)~=0);
endi = [temp;length(BreathDataTable1.Time0)];
starti = [1;temp+1];
clear temp

%%
NN = size(DataEventHypnog_Mat,1);

VEseries = NaN*ones(NN,1);
CountSeries = zeros(NN,1);

for i=1:length(starti)
    timestart = BreathDataTable1.Time_start(starti(i):endi(i));
    timestart_ = BreathDataTable1.Time_end(endi(i));
    timestart = [timestart;timestart_];
    VI = [BreathDataTable1.VE(starti(i):endi(i));NaN];
    temp2 = DataEventHypnog_Mat(:,1);
    temp = interp1(timestart,VI,DataEventHypnog_Mat(:,1),'previous'); %upsample to staircase from breath domain
    
    I=find(~isnan(temp));
    %        tempC = CountSeries(I);
    %         if sum(tempC)==0
    win = ones(length(I),1);
    %         end
    
    I2=find(~isnan(temp)&~isnan(VEseries)); %indices of overlap between existing signal (being build) and new window signal
    if ~isempty(I2)
        i1 = I2(1);
        i2 = I2(end)+1;
        m = 1/(i2-i1); c = -m*i1;
        win_ = I2*m + c;
        %             temp__ = win(I2)
        
        tempi= I2-I(1)+1;
        win(tempi)=win_;
    end
    
    VEseriesNew = nansum([VEseries(I).*(1-win) , temp(I).*win],2); %calculate weighted average of two VE signals
    VEseries(I)=VEseriesNew; %copy into existing signal
    
    if 0
        figure(1); clf(1);
        ax(1)=subplot(3,1,1);
        plot(DataEventHypnog_Mat(:,1),temp,'.');
        hold on
        stairs(timestart,VI);
        ax(2)=subplot(3,1,2);
        plot(DataEventHypnog_Mat(:,1),VEseries);
        ax(3)=subplot(3,1,3);
        plot(DataEventHypnog_Mat(:,1),CountSeries);
        linkaxes(ax,'x');
        xlim([min(timestart)-1200 max(timestart)+120]);
    end
    %pause(2)
end

%%
DataEventHypnog_Mat = [DataEventHypnog_Mat VEseries];
ChannelsList = [ChannelsList {'VE'}];
ChannelsFs = [ChannelsFs;ChannelsFs(1)];

%% Calculate Eupnea
filter_HFcutoff_butter0 = 1/(120*2*pi);
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(ChannelsFs(1)/2),'low');

signalfiltered = nanfilter(B_butter0,A_butter0,VEseries,1);

if sum(strcmp(ChannelsList,'VEeupnea'))==0
    DataEventHypnog_Mat = [DataEventHypnog_Mat signalfiltered];
    ChannelsList = [ChannelsList {'VEeupnea'}];
    ChannelsFs = [ChannelsFs;ChannelsFs(1)];
    clear signalfiltered
else
    idx = find(strcmp(ChannelsList,'VEeupnea')==1);
    DataEventHypnog_Mat(:,idx) = signalfiltered;
    clear signalfiltered
end

%%
VE_peupnea = 100*DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'VE')))./DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'VEeupnea')));
%to do:
%make breath-domain data and divie there, then upsample to stairs again,
%will keep each breath value flat/constant for the breath duration

[ApneaHypopneamat2]=FindApneaHypopnea(VE_peupnea,DataEventHypnog_Mat,ChannelsFs(1),ChannelsList);



if 0
figure(11); clf(11);
ax(1)=subplot(3,1,1);
plot(DataEventHypnog_Mat(:,1),ApneaHypopneamat2,'.');
hold on
stairs(timestart,VI);
ax(2)=subplot(3,1,2);
plot(DataEventHypnog_Mat(:,1),VEseries);
ax(3)=subplot(3,1,3);
plot(DataEventHypnog_Mat(:,1),CountSeries);
linkaxes(ax,'x');
xlim([min(timestart)-1200 max(timestart)+120]);
end


%% save the autoscored events
if 0
    if sum(strcmp(ChannelsList,'VEpeupnea'))==0
        DataEventHypnog_Mat = [DataEventHypnog_Mat VE_peupnea];
        ChannelsList = [ChannelsList {'VEpeupnea'}];
        ChannelsFs = [ChannelsFs;ChannelsFs(1)];
        clear VE_peupnea
    else
        idx = find(strcmp(ChannelsList,'VEpeupnea')==1);
        DataEventHypnog_Mat(:,idx) = VE_peupnea;
        clear VE_peupnea
    end
    
    if sum(strcmp(ChannelsList,'EventRespAuto'))==0
        DataEventHypnog_Mat = [DataEventHypnog_Mat ApneaHypopneamat2];
        ChannelsList = [ChannelsList {'EventRespAuto'}];
        ChannelsFs = [ChannelsFs;ChannelsFs(1)];
        clear  ApneaHypopneatemp
    else
        idx = find(strcmp(ChannelsList,'EventRespAuto')==1);
        DataEventHypnog_Mat(:,idx) = ApneaHypopneamat2;
        clear  ApneaHypopneatemp
    end
else
    disp('Warning: Not saving autoscoring to DataEventHypnog_Mat');
end



% position=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Position')==1));


%% for sao2/arousals

IdxEd=find([diff(ApneaHypopnea_)<0; 0]);
IdxSt=find([diff(ApneaHypopnea_)>0; 0]);



Time_new=-100:1:100;
Time_OtherSigs=repmat(Time_new,length(IdxEd),1)+repmat(evtEndtmp,1,length(Time_new)); %table of times +/- 100s surrounding events

SaO2Ve=interp1(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1)),Time_OtherSigs,'nearest'); %table of signal matching the above time array
Sao2Ensmb=nanmean(SaO2Ve,2); %ensemble averaged signal;

ArVe=interp1(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1)),Time_OtherSigs,'nearest'); %table of signal matching the above time array
ArEnsmb=nanmean(ArVe,2); %ensemble averaged signal;

%for each event: type, tstart, tend, duration, state, position, desat, arousalYN, average depth, max depth.






