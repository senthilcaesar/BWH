function [DataEventHypnog_Mat,ChannelsList,ChannelsFs]=runWSanalysis(DataEventHypnog_Mat,ChannelsList,ChannelsFs)

global settings

%% remove old versions of signals
ChannelsListNew = {'WPr','ArPr','EventsArWS','WPrB','ArPrB','EventsArWSB'};
M = string(ChannelsList)==string(ChannelsListNew');
temp = repmat(1:length(ChannelsList),length(ChannelsListNew),1);
temp = M.*temp;
temp(temp==0)=NaN;
temp = max(temp')
temp(isnan(temp))=[];
I=temp;

ChannelsList(I)=[];
DataEventHypnog_Mat(I,:)=[];
ChannelsFs(I)=[];

%% WS analysis
disp('WS analysis');
    load mdlA
    load mdlNoise
    load RefTable
    load mdlAcc
    
    if ~isfield(settings,'WSArVersion')
        settings.WSArVersion=2; %default to version 2 from 1/7/2020 2:10pm.
    end
    
    if settings.WSArVersion==1
        load mdlAR
    elseif settings.WSArVersion==2
        load mdlARwsbalance
    end
   
    tic
    [WPr,ArPr,WPrB,ArPrB,WSinfo] = RunWS(DataEventHypnog_Mat,ChannelsList,(1-settings.scoredarousalsinwake),mdlA,mdlNoise,RefTable,mdlAcc,mdlAR);
    ArSig = max([WPr ArPr]')' > 0.5;
    %ArSig = WPr > 0.25;
    minarlength=2.5;
    EventsArWS = RemoveShortSegments(ArSig,minarlength,1./settings.Fs,1); % downsampling first would lead to some speed up here.
    ArSigB = max([WPrB ArPrB]')' > 0.5;
    %ArSigB = WPrB > 0.25;
    EventsArWSB = RemoveShortSegments(ArSigB,minarlength,1./settings.Fs,1);
    toc
    
    if 1
    figure(33); clf(33);
    ax1(1)=subplot(3,1,1);
    ch = find(strcmp(ChannelsList,'Epochs')==1);
    plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,ch)/4,'linewidth',2);
    hold on;
    ch = find(strcmp(ChannelsList,'EventsAr')==1);
    plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,ch),'g','linewidth',2);
    plot(DataEventHypnog_Mat(:,1),max([WPr ArPr]')','color',[0.5 0.5 0.5]);
    %plot(DataEventHypnog_Mat(:,1),max([WPrB ArPrB]')');
    plot(DataEventHypnog_Mat(:,1),WPr,'color',[0.5 0.5 0.9]);
    plot(DataEventHypnog_Mat(:,1),EventsArWS,'r','linewidth',1);
    %plot(DataEventHypnog_Mat(:,1),ArSig);
    ax1(2)=subplot(3,1,2);
    ch = find(strcmp(ChannelsList,'EEG1')==1);
    ch2 = find(strcmp(ChannelsList,'EEG3')==1);
    plot(DataEventHypnog_Mat(:,1),[DataEventHypnog_Mat(:,ch) DataEventHypnog_Mat(:,ch2)+nanstd(DataEventHypnog_Mat(:,ch))]);
    ax1(3)=subplot(3,1,3);
    ch = find(strcmp(ChannelsList,'Flow')==1);
    plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,ch));
    linkaxes(ax1,'x')
    
    diff(get(gca,'xlim'));
    temp = get(gca,'xlim');
    set(gca,'xlim',temp);
    end
    
    DataEventHypnog_Mat = [DataEventHypnog_Mat,WPr,ArPr,EventsArWS,WPrB,ArPrB,EventsArWSB];
    ChannelsList = [ChannelsList,'WPr','ArPr','EventsArWS','WPrB','ArPrB','EventsArWSB'];
    ChannelsFs = [ChannelsFs;settings.Fs;settings.Fs;settings.Fs;settings.Fs;settings.Fs;settings.Fs];
    
    WakeSleepInfo.WSinfo = WSinfo;