global settings
settings.WSanalysis=1;
settings.scoredarousalsinwake=1;
settings.PlotWSBalanceFigure=1;
dt = DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1);
F_samp = 1/dt;

%%
if isfield(settings,'WSanalysis') && settings.WSanalysis==1
        disp('WS analysis');
        load mdlA
        load mdlNoise
        load RefTable
        load mdlAcc
    
    settings.WSArVersion=2;
    
        if settings.WSArVersion==1
            load mdlAR
        elseif settings.WSArVersion==2
            load mdlARwsbalance
        end

        tic
        [WPr,ArPr,WPrB,ArPrB,WSinfo] = RunWS(DataEventHypnog_Mat,ChannelsList,(1-settings.scoredarousalsinwake),mdlA,mdlNoise,RefTable,mdlAcc,mdlAR);

        %takes about 1 min for 2xRemoveShortSegments (@100Hz)
        ArSig = max([WPr ArPr]')' > 0.5;
        %ArSig = WPr > 0.25;
        minarlength=2.5;
        EventsArWS = RemoveShortSegments(ArSig,minarlength,1./F_samp,1); % downsampling first would lead to some speed up here.
        ArSigB = max([WPrB ArPrB]')' > 0.5;
        %ArSigB = WPrB > 0.25;
        EventsArWSB = RemoveShortSegments(ArSigB,minarlength,1./F_samp,1);
        toc
    
    ArPr_ = ArPr;
    ArPrB_ = ArPrB;
    WPr_ = WPr;
    EventsArWS_ = EventsArWS;
    
        settings.WSArVersion=1;

        if settings.WSArVersion==1
            load mdlAR
        elseif settings.WSArVersion==2
            load mdlARwsbalance
        end

        tic
        [WPr,ArPr,WPrB,ArPrB,WSinfo] = RunWS(DataEventHypnog_Mat,ChannelsList,(1-settings.scoredarousalsinwake),mdlA,mdlNoise,RefTable,mdlAcc,mdlAR);

        %takes about 1 min for 2xRemoveShortSegments (@100Hz)
        ArSig = max([WPr ArPr]')' > 0.5;
        %ArSig = WPr > 0.25;
        minarlength=2.5;
        EventsArWS = RemoveShortSegments(ArSig,minarlength,1./F_samp,1); % downsampling first would lead to some speed up here.
        ArSigB = max([WPrB ArPrB]')' > 0.5;
        %ArSigB = WPrB > 0.25;
        EventsArWSB = RemoveShortSegments(ArSigB,minarlength,1./F_samp,1);
        toc
    
end

%%
    if settings.PlotWSBalanceFigure==1
    
    figure(133); clf(133);
    ax1(1)=subplot(4,1,1:2);
    ch = find(strcmp(ChannelsList,'Epochs')==1);
    plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,ch)/4,'linewidth',2);
    hold on;
    
    if exist('ArPr_') && exist('WPr_')
    plot(DataEventHypnog_Mat(:,1),2+max([WPr_ ArPr_]')','color',[0.5 0.5 0.5]);    
    end
    plot(DataEventHypnog_Mat(:,1),max([WPr ArPr]')','color',[0.5 0.5 0.5]);
    %plot(DataEventHypnog_Mat(:,1),max([WPrB ArPrB]')');
    
    if exist('WPr_')
        plot(DataEventHypnog_Mat(:,1),2+WPr_,'color',[0.5 0.5 0.9]);    
    end
    plot(DataEventHypnog_Mat(:,1),WPr,'color',[0.5 0.5 0.9]);
    
    
    ch = find(strcmp(ChannelsList,'EventsAr')==1);
    if exist('EventsArWS_')
    plot(DataEventHypnog_Mat(:,1),2+DataEventHypnog_Mat(:,ch),'g','linewidth',2);
    end
    
    plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,ch),'g','linewidth',2);
    
    if exist('EventsArWS_')
    plot(DataEventHypnog_Mat(:,1),2+EventsArWS_*0.9+0.05,'r','linewidth',1); 
    end
    plot(DataEventHypnog_Mat(:,1),EventsArWS*0.9+0.05,'r','linewidth',1);
    
    
    %plot(DataEventHypnog_Mat(:,1),ArSig);
    ax1(2)=subplot(4,1,3);
    ch = find(strcmp(ChannelsList,'EEG1')==1);
    ch2 = find(strcmp(ChannelsList,'EEG3')==1);
    plot(DataEventHypnog_Mat(:,1),[DataEventHypnog_Mat(:,ch) DataEventHypnog_Mat(:,ch2)+nanstd(DataEventHypnog_Mat(:,ch))]);
    ax1(3)=subplot(4,1,4);
    ch = find(strcmp(ChannelsList,'Flow')==1);
    plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,ch));
    linkaxes(ax1,'x')
    
    diff(get(gca,'xlim'));
    temp = get(gca,'xlim');
    set(gca,'xlim',temp);
    end
    
    DataEventHypnog_Mat = [DataEventHypnog_Mat,WPr,ArPr,EventsArWS,WPrB,ArPrB,EventsArWSB];
    ChannelsList = [ChannelsList,'WPr','ArPr','EventsArWS','WPrB','ArPrB','EventsArWSB'];
    
    ChannelsFs = [ChannelsFs;F_samp;F_samp;F_samp;F_samp;F_samp;F_samp];
    WakeSleepInfo.WSinfo = WSinfo;
    %add signals back to DataEventHypnog_Mat ChannelsList ChannelsFs
    

