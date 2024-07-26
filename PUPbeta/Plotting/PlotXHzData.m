%function PlotXHzData()

% run StartHere.m first

addpath(genpath(pwd));
if exist('settings') && isfield(settings,'PlotUsingDateTime') && settings.PlotUsingDateTime==1
    usedatetime=1;
else
    usedatetime=0;
end
%usedatetime=1;
FilterFlow25=1;

try
    ChannelsList=SigT.Properties.VariableNames;
end

if 0
    AddFilterTo = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10'};
else
    AddFilterTo = {};
end

Filter1 = [0.3 45 1];
    if Filter1(2)>0.999*ChannelsFs(1)/2, Filter1(2)=0.999*ChannelsFs(1)/2, end


if ~exist('plotchanneltext')
    plotchanneltext = { ...
        {'Epochs'},{'Flow'},{'Position'},{'EventsResp'},{'EventsAr','EventsArWS'},{'WArPr'},{'SpO2'},{'EEG11'},{'HR'}...
        };
    %     {'Epochs','EventsResp','EventsAr'},{'Flow'},{'SpO2'},{'Edi'},{'GGpmax'},{'TonicREM','PhasicREM'}...
    %     }; ,{'WArPrB','EventsArWSB'}
    
    %     plotchanneltext = { ...
    %         {'Epochs'},{'EventsResp'},{'WArPr','EventsAr'},{'Flow'},{'SpO2'},{'Abdomen'},{'Position'},{'Pmask'}...
    %         };
    % plotchanneltext = { ...
    %     {'EventsAr'},{'WSintensity'},{'PrWake'},{'Pmask'},{'Flow'},{'Edi'},{'FlowDrive'},{'VEseries'},{'SpO2'}...
    %     };
    %     plotchanneltext = { ...
    %       {'PrWake'},{'Flow'},{'SpO2'},{'VI','Vdr_est'}...
    %      };
    %  plotchanneltext = { ...
    %      {'Epochs'},{'EventsAr'},{'WPr','ArPr','EventsArWS'},{'Pbetalogfilt1','Palphalogfilt1','Pthetalogfilt1','Pdeltalogfilt1'},{'EEG1'},{'Flow'}...
    %      }; %'WPr','ArPr','EventsArWS'
    
    % = { ...
    %     {'Epochs'},{'Flow'},{'VE','VEeupnea'},{'VEpeupnea'},{'EventRespAuto'},{'EventsResp'} ...
    %     {'SpO2'}...
    %     };
    %{'EventsAr','WakeSleep'}, Pbetalogfilt
    
    % Basic view
    % plotchanneltext = { ...
    %    {'Epochs','EventsResp','EventsAr'},{'Flow'},{'Thorax'},{'Abdomen'},{'SpO2'}};
end

if ~exist('SigT') % convert DataEventHypnog_Mat to SigT
    SigT=array2table(DataEventHypnog_Mat);
    SigT.Properties.VariableNames = ChannelsList;
    clear DataEventHypnog_Mat
end

try
    I_=find(strcmp(ChannelsList,'WArPr')==1);
    if isempty(I_)
        SigT.WArPr=max([SigT.WPr SigT.ArPr]')';

        logit = @(p) log(p./(1-p));
        WSBalance = logit(SigT.WPr);
        WSBalance(WSBalance>10)=10;
        WSBalance(WSBalance<-10)=-10;
        SigT.WSBalance = WSBalance;
        
        ChannelsList = SigT.Properties.VariableNames;
    end
end

try
    I_=find(strcmp(ChannelsList,'WArPrB')==1);
    if isempty(I_)
        SigT.WArPrB=max([SigT.WPrB SigT.ArPrB]')';

        logit = @(p) log(p./(1-p));
        WSBalance = logit(SigT.WPrB);
        WSBalance(WSBalance>10)=10;
        WSBalance(WSBalance<-10)=-10;
        SigT.WSBalanceB = WSBalance;
        
        ChannelsList = SigT.Properties.VariableNames;
    end
end



%% noise plot
try
    I_=find(strcmp(ChannelsList,'Noise')==1);
    if isempty(I_)
        SigT.Noise = FlowSignalToNoise(SigT.Time,SigT.Flow,0);
        ChannelsList = SigT.Properties.VariableNames;
    end
end

%Find channel numbers; ignore channels that do not exist
plotchannels=[];
ChannelEmpty=zeros(length(plotchanneltext),1);
for i=1:length(plotchanneltext)
    temp = plotchanneltext{i};
    plotchannels{i}=[];
    for j=1:length(temp)
        I=find(strcmp(ChannelsList,temp{j})==1);
        if ~isempty(I)
            plotchannels{i}(j)=I;
        else
            ChannelEmpty(i)=1;
        end
    end
end
plotchannels(ChannelEmpty==1)=[];

global ax2
% ax2=[];

figure(412); clf(412);

if usedatetime
    t = datetime(SigT.Time/86400,'ConvertFrom','excel');
end
for i=1:length(plotchannels)
    try
    ax2(i)=subplot(length(plotchannels),1,i);
    Ydata = SigT{:,plotchannels{i}};
    
    if exist('plotTransformLogit')
        if plotTransformLogit(i)==1
            Ydata = logit(Ydata);
        end
    end
    
    
    if sum(string(ChannelsList(plotchannels{i}))=="Flow") && FilterFlow25
        TotalTime = (height(SigT)-1)*(1/ChannelsFs(1));
        Time_ = [0:(1/ChannelsFs(1)):TotalTime]';
        TimeDS=[0:(1/25):TotalTime]';
        FlowDS = interp1(Time_,Ydata,TimeDS,'linear'); % downsample
        Ydata = interp1(TimeDS,FlowDS,Time_,'linear');  % upsample back to original length.
        clear Time_ TimeDS FlowDS
    end
    
    for j=1:size(plotchannels{i})
        if sum(string(AddFilterTo')==string(ChannelsList(plotchannels{i}(j))))
            %found match
            [B_butter0,A_butter0] = butter(Filter1(3),[Filter1(1) Filter1(2)]/(ChannelsFs(1)/2));
            Ydata(:,j) = nanfilter(B_butter0,A_butter0,Ydata(:,j),1); %filtfilt, otherwise flow signal is right-shifted
        end
    end
    
    if usedatetime
        t = datetime(SigT.Time/86400,'ConvertFrom','excel');
        plot(t,Ydata);
    else
        plot(SigT.Time,Ydata);
    end
    if i<length(plotchannels)
        set(gca,'xtick',[],'xcolor',[1 1 1]);
    end
    set(gca,'box','off','tickdir','out');
    postemp = get(gca,'position');
    postemp2 = postemp;
    postemp2(1)=0.09;
    postemp2(3)=0.90;
    postemp2(4)=1.2*postemp(4);
    set(gca,'position',postemp2);
    ylabel(ChannelsList(plotchannels{i}),'fontsize',8);
    if sum(string(ChannelsList(plotchannels{i}))=="Epochs")
        set(gca,'ytick',[0 1 2 3 4 8],'yticklabels',{'N3','N2','N1','R','W','?'});
    end
    if sum(string(ChannelsList(plotchannels{i}))=="EventsResp")
        set(gca,'ytick',[2 3 4 5 6],'yticklabels',{'OA','CA','OH','MA','CH'});
    end
    end
end
warning('off');
linkaxes(ax2,'x');
clear i plotchannels postemp postemp2 plotchanneltext temp I j
set(gcf,'color',[1 1 1]);
%%
global xvalues yvalues range h1 sliderstep P
xvalues=[];
yvalues=[];
timeminmax = [SigT.Time(1) SigT.Time(end)];
if usedatetime
    timeminmax = datetime(timeminmax/86400,'ConvertFrom','excel');
end
range = diff(timeminmax)*1.0;

plotwithslider(timeminmax);

