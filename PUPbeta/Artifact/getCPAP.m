function [CPAPoff,CPAP] = getCPAP(SigT,settings,ChannelsList,ploton)

if ~exist('ploton')
    try
        ploton = settings.plotfigure;
    catch me
        ploton=0;
    end
end

Time = SigT.Time;
N=length(Time);
if (isempty(find(strcmp(ChannelsList,'CPAP')==1)) && ...
    isempty(find(strcmp(ChannelsList,'Pmask')==1))) || ...
    settings.ignoreCPAPdata
    CPAP = zeros(N,1);
    Pmask = zeros(N,1); %for plot
    CPAPoff = ones(N,1);
else
    %CPAP = DataEventHypnog_Mat(:,ColumnHeads(find(strcmp(ChannelsList,'CPAP')==1))); % not indexed, 'find' stays
    Pmask = SigT.Pmask; % not indexed, 'find' stays
    Pmask(isnan(Pmask))=0; % DLM replaced NaN's with zeros, because can't filtfilt on NaN data
    dt = 1/settings.Fs;
    
    %% Find when CPAP is off
    filter_HFcutoff_butter0 = 1/10; %4Hz=250ms
    filter_order0 = 4;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    CPAP = filtfilt(B_butter0,A_butter0,Pmask);
    CPAP(CPAP<-20)=-20; CPAP(CPAP>20)=20;
    CPAPrange = [-settings.minabsPmaskforCPAPoff settings.minabsPmaskforCPAPoff];
    CPAPoff = CPAP<CPAPrange(2)&CPAP>CPAPrange(1);
    
    %plot(downsample(Time,dsf),2*downsample(CPAPoff,dsf),'g');
    minlength=180; minlengthi = (minlength/dt);
    diffCPAPoff1=diff(CPAPoff);
    CPAPswitch = find(diffCPAPoff1==1|diffCPAPoff1==-1);
    CPAPswitchI = [1;CPAPswitch;length(CPAP)];
    difflengths = diff(CPAPswitchI);
    CPAPoff_ = mod(1+CPAPoff(1)+(1:length(difflengths)),2)'; %
    temp = [difflengths,CPAPoff_];
    Iremove = find(difflengths<minlengthi&CPAPoff_==1);
    for i=1:length(Iremove)
        CPAPoff(CPAPswitchI(Iremove(i)):CPAPswitchI(Iremove(i)+1))=0;
    end
    
end
%% Plot
if ploton
    dsf=100;
    %nicefig = @() set(gca,'box','off','tickdir','out','xcolor',[1 1 1],'xticklabel',[],'fontsize',8,'fontname','arial narrow');
    nicefiglast = @() set(gca,'box','off','tickdir','out','fontsize',8,'fontname','arial narrow');
    
    figure(101);
    ax1(1)=subplot(4,1,1); plot(downsample(Time,dsf),downsample(Pmask,dsf));  ylabel('CPAP');
    hold('on');
    plot(downsample(Time,dsf),downsample(CPAP,dsf),'r');
    
    plot(downsample(Time,dsf),5*downsample(1-CPAPoff,dsf),'g');
    hold('off');
    nicefiglast();
end
