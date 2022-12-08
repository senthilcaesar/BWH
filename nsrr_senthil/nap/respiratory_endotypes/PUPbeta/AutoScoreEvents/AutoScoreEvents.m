function [SigT,ChannelsList,ChannelsFs,Evts] = AutoScoreEvents(BreathDataTable,SigT,ChannelsList,ChannelsFs,Evts,plotfig,minEventDuration)
global settings

    if ~exist('minEventDuration')==1 || isempty('minEventDuration')==1
        minEventDuration=10;
    end
    
% addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta20190629'))
[BreathDataTable1,~]=GetNonOvlappedVE(BreathDataTable);

%find indices at which data starts new windows
temp = find(diff(BreathDataTable1.Time0)~=0);
endi = [temp;length(BreathDataTable1.Time0)];
starti = [1;temp+1];
clear temp

%%
NN = size(SigT,1);

VEseries = NaN*ones(NN,1);
CountSeries = zeros(NN,1);

for i=1:length(starti)
%     i
    timestart = BreathDataTable1.Time_start(starti(i):endi(i));
    timestart_ = BreathDataTable1.Time_end(endi(i));
    timestart = [timestart;timestart_];
    VI = [BreathDataTable1.VE(starti(i):endi(i));NaN];
    temp2 = SigT{:,1};
    temp = interp1(timestart,VI,SigT.Time,'previous'); %upsample to staircase from breath domain
    
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
        figure(101); clf(101);
        ax(1)=subplot(3,1,1);
        plot(SigT.Time,temp,'.');
        hold on
        stairs(timestart,VI);
        ax(2)=subplot(3,1,2);
        plot(SigT.Time,VEseries);
        ax(3)=subplot(3,1,3);
        plot(SigT.Time,CountSeries);
        linkaxes(ax,'x');
        xlim([min(timestart)-1200 max(timestart)+120]);
    end
    %pause(2)
end

%%
if sum(strcmp(ChannelsList,'VE'))==0
    SigT.VE = VEseries;
    ChannelsList = [ChannelsList {'VE'}];
    ChannelsFs = [ChannelsFs;ChannelsFs(1)];
    clear signalfiltered
else
    idx = find(strcmp(ChannelsList,'VE')==1);
    SigT.VE  = VEseries;
    clear signalfiltered
end

%% Calculate Eupnea
filter_HFcutoff_butter0 = 1/(120*2*pi);
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(ChannelsFs(1)/2),'low');

signalfiltered = nanfilter(B_butter0,A_butter0,VEseries,1);
%should be nanfilter(B,A,Signal,...) If there is an error here there may be
%an older version of nanfilter that has inputs in a different order
%(Signal,B,A). Fix nanfilter, not this code.

if sum(strcmp(ChannelsList,'VEeupnea'))==0
    SigT.VEeupnea  = signalfiltered;
    ChannelsList = [ChannelsList {'VEeupnea'}];
    ChannelsFs = [ChannelsFs;ChannelsFs(1)];
    clear signalfiltered
else
    
    SigT.VEeupnea = signalfiltered;
    clear signalfiltered
end

%%
VE_peupnea = 100*SigT.VE./SigT.VEeupnea;
%to do:
%make breath-domain data and divie there, then upsample to stairs again,
%will keep each breath value flat/constant for the breath duration

[EventsRespAuto]=FindApneaHypopnea(VE_peupnea,SigT,ChannelsFs(1),ChannelsList,minEventDuration);

if 1
figure(102); clf(102);
ax(1)=subplot(3,1,1);
plot(SigT.Time,EventsRespAuto,'.');
hold on
stairs(timestart,VI);
ax(2)=subplot(3,1,2);
plot(SigT.Time,VEseries);
ax(3)=subplot(3,1,3);
plot(SigT.Time,CountSeries);
linkaxes(ax,'x');
xlim([min(timestart)-1200 max(timestart)+120]);

figure(66); clf(66);
yyaxis right
plot(SigT.Time,EventsRespAuto,'Color',[1 0.5 0.6]); % events in sleep only
yyaxis left
plot(SigT.Time,SigT.Flow);
legend('flow','events in sleep');
xlabel('time(s)');
end
%% save the autoscored events
if 1
    if sum(strcmp(ChannelsList,'VEpeupnea'))==0
        SigT.VEpeupnea = VE_peupnea;
        ChannelsList = [ChannelsList {'VEpeupnea'}];
        ChannelsFs = [ChannelsFs;ChannelsFs(1)];
    else
        SigT.VEpeupnea = VE_peupnea;
    end
    
    if sum(strcmp(ChannelsList,'EventsRespAuto'))==0
        SigT.EventsRespAuto = EventsRespAuto;
        ChannelsList = [ChannelsList {'EventsRespAuto'}];
        ChannelsFs = [ChannelsFs;ChannelsFs(1)];
    else
        SigT.EventsRespAuto = EventsRespAuto;
    end
else
    disp('Warning: Not saving autoscoring to SigT');
end
% position=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Position')==1));


%for each event: type, tstart, tend, duration, state, position, desat, arousalYN, average depth, max depth.


end



