function [EMGdi,Err]=EdiAnalysisSpike(f,EdiOpt,SignalFormat)

% This code uses EditDist to match up this signals and channels
% These are the names of the signals I am going to analyze
sigs = {'Flow','Edi1','Edi2','Edi3','Edi4','Edi5'}; %,'Pes'
%index into sigs using {}, not ()

if strcmp(SignalFormat,'Luo')
    numSigs = length(sigs);
    
    numChans = length(f.channel_specs);  %number of channels collected by Luo
    %numChans = length(f);
    chan = cell(1,numChans);
    
    for n = 1:numChans
        chan{n} = f.channel_specs(1,n).name;
    end
    
    sigName = cell(1,numSigs);  %the names of the signals based on what Luo named them
    
    % find the "distance" between two strings. If this hits a bug, then just
    % do it by hand below
    distBwStr = NaN(numSigs,numChans);
    for n = 1:numSigs
        for m = 1:numChans
            distBwStr(n,m) = EditDist(sigs{n},chan{m}); %external function on fileexchange
        end
        [~,idx] = min(distBwStr(n,:));
        sigName{n} = chan{idx};
    end
    
    for i=1:length(sigs)
        eval([sigs{i} '=importLabChartData(f,sigName{i},f.n_records);']);
    end
    
elseif strcmp(SignalFormat,'Spike')
    %read directly from spike channel output and assing to channel labels
    for i=1:length(sigs)
        eval([sigs{i}  '=[f.(sigs{i}).values];']);
    end
end

clear chan sigName numsigs numChans n m idx distBwStr

FlowScaleFactor=round(length(Edi1)/length(Flow)); %work out ratio of EMG/flow channel sampling

Flow=-Flow/60;

N = length(Edi1);%(Flow) %changed as Flow different sample rate

if strcmp(SignalFormat,'Luo')
    dt = f.channel_specs(1,1).dt(1);
elseif strcmp(SignalFormat, 'Spike')
    dt=f.Edi1.interval;%changed to read from Spike EMG output
end

Fs = 1/dt;
Time = (0:dt:(N-1)*dt)'; %Time now for EMG signals only

EdiSum = Edi1 + Edi2 + Edi3 + Edi4 + Edi5;
Temp = EdiSum;

% Find noise periods where the data is invalid
if strcmp(SignalFormat,'Luo')
    Vflow_ds = downsample(Flow,20); %need to downsample the signal or it's too slow
    time_ds = downsample(Time,20);
elseif strcmp(SignalFormat, 'Spike')
    Vflow_ds =Flow; %For nox, flow already downsampled
    time_ds = downsample(Time,FlowScaleFactor);
    if length(time_ds) ~= length(Vflow_ds)
        dtFlow = time_ds(2) - time_ds(1);
        time_ds = (0:dtFlow:(length(Flow)-1)*dtFlow)';
    end
end



noisewav_ds = FlowSignalToNoise(time_ds,Vflow_ds,0);
noisewav = interp1(time_ds,noisewav_ds,Time,'nearest','extrap');



%%
Nwin=ceil(N*dt/60);

%EdiOpt = 6;
switch EdiOpt
    case 1
        ratios = [1 0 0 0 0];
    case 2
        ratios = [0 1 0 0 0];
    case 3
        ratios = [0 0 1 0 0];
    case 4
        ratios = [0 0 0 1 0];
    case 5
        ratios = [0 0 0 0 1];
    case 6
        ratios = [-1 1 2 1 -1];
end
% ratios = [-1 0 1 2 1];
% ratios = [1 2 1 0 -1];
% ratios = [-1 1 2 1 -1];
ratios = ratios./sum(abs(ratios));
signalT = ratios(1)*Edi1+ratios(2)*Edi2+ratios(3)*Edi3+ratios(4)*Edi4+ratios(5)*Edi5;

clear Edi1 Edi2 Edi3 Edi4 Edi5 Pes
%%
indrows=(1:N)';
wint = 60;
M = round(wint/dt);
ioverlap = round(0.2*M); %20% overlap
indarray = buffer(indrows,M,ioverlap,'nodelay');
clear indrows;
Time_ = buffer(Time,M,ioverlap,'nodelay');
Nwin = size(indarray,2);

% indarray(indarray==0)=NaN;

windowwidth=0.025;
Foverlap=0.8; %sec, check Foverlap*windowwidth/dt is an integer
temp = Time(indarray(:,1));
tempout = buffer(temp,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
Timerms = tempout(round(size(tempout,1)/2),:);
N2 = length(Timerms);

EdiRMS_ = NaN*indarray;


%%

plotonfly=0;
for nn=1:(Nwin-1)
    
    nn
            
    temp = Time(indarray(:,nn));
    tempout = buffer(temp,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
    Timerms = tempout(round(size(tempout,1)/2),:);
    
    I = indarray(:,nn);
    I_ds=indarray((1:FlowScaleFactor:end),nn); %added downsampled array for Flow timepoints
    signal = signalT(I);
    
    % condition to ignore Edi data if it's all 0 (happens when there is a
    % break in data collection)
    if sum(EdiSum(I)) == 0
        continue
    end
    
    if sum(noisewav(I)>=1)>0 %>0.9*length(signal) %
        %[num2str(nn) ' is noise']
        continue
    end
    
    if plotonfly
        figure(1); clf(1);
        Nsubs=3;
        ax(1)=subplot(Nsubs,1,1); plot(Time(I_ds),Flow(round(I_ds./FlowScaleFactor))); set(gca,'xtick',[]); %Adjusted flow timepoints
        ax(2)=subplot(Nsubs,1,2); plot(Time(I),signal); set(gca,'xtick',[]); hold('on');
        ax(2)=subplot(Nsubs,1,2); plot(Time(I),EdiSum(I)-std(EdiSum(I))); set(gca,'xtick',[]); hold('on');
        subploti = 2;
    end
    
    ECG_peak_i = EKGpeakdetection(EdiSum(I),Time(I),dt,1,0);
    %storedECG__peak(nn
    
    RRnormal=median(diff(ECG_peak_i)*dt);
    if plotonfly
        ax(Nsubs)=subplot(Nsubs,1,subploti); plot(Time(I(ECG_peak_i)),signal((ECG_peak_i)),'r.');
        linkaxes(ax,'x');
    end
    
    EKGkeepFr = [0.2 0.8];
    
    temp = diff(ECG_peak_i);
    li=ECG_peak_i(1:end-1)+round(EKGkeepFr(1)*temp);
    ri=ECG_peak_i(1:end-1)+round(EKGkeepFr(2)*temp);
    EKGfree = 0*EdiSum(I);
    for i=1:length(li)
        EKGfree(li(i):ri(i))=1;
    end
    temp = signal;
    temp(EKGfree==0)=NaN;
    
    if plotonfly
        ax(subploti)=subplot(Nsubs,1,subploti); plot(Time(I),temp); set(gca,'xtick',[]);
    end
    cutoff = 40;
    filter_order = 4;
    [B_butter,A_butter] = butter(filter_order,[cutoff]/(Fs/2),'high');
    signalfilt = filtfilt(B_butter,A_butter,signal);
    
    lefttemplatetime = RRnormal*0.2; righttemplatetime = RRnormal*0.2;
    leftdi=round(1/dt*lefttemplatetime); rightdi=round(1/dt*righttemplatetime);
    polyorder=3;
    if plotonfly
        templateonly=1;
        [~,mediantemplate] = crosscontaminated(signalfilt,ECG_peak_i,leftdi,rightdi,templateonly,polyorder);
    end
    [signalfiltdecontaminated,mediantemplate] = crosscontaminated(signalfilt,ECG_peak_i,leftdi,rightdi,0,polyorder);
    
    lefttemplatetime = RRnormal*0.33; righttemplatetime = RRnormal*0.0;
    leftdi=round(1/dt*lefttemplatetime); rightdi=round(1/dt*righttemplatetime);
    templateonly=1;
    [~,mediantemplate] = crosscontaminated(signalfiltdecontaminated,ECG_peak_i,leftdi,rightdi,templateonly,polyorder);
    if plotonfly
        [signalfiltdecontaminated,mediantemplate] = crosscontaminated(signalfiltdecontaminated,ECG_peak_i,leftdi,rightdi,0,polyorder);
    end
    
    if 1
        lefttemplatetime = RRnormal*0.5; righttemplatetime = RRnormal*0.5;
        leftdi=round(1/dt*lefttemplatetime); rightdi=round(1/dt*righttemplatetime);
        if plotonfly
            templateonly=1;
            [~,mediantemplate] = crosscontaminated(signalfiltdecontaminated,ECG_peak_i,leftdi,rightdi,templateonly,polyorder);
        end
        [signalfiltdecontaminated,mediantemplate] = crosscontaminated(signalfiltdecontaminated,ECG_peak_i,leftdi,rightdi,0,polyorder);
    end
    
    if 1
        lefttemplatetime = RRnormal*0.0; righttemplatetime = RRnormal*0.75;
        leftdi=round(1/dt*lefttemplatetime); rightdi=round(1/dt*righttemplatetime);
        if plotonfly
            templateonly=1;
            [~,mediantemplate] = crosscontaminated(signalfiltdecontaminated,ECG_peak_i,leftdi,rightdi,templateonly,polyorder);
        end
        [signalfiltdecontaminated,mediantemplate] = crosscontaminated(signalfiltdecontaminated,ECG_peak_i,leftdi,rightdi,0,polyorder);
    end
    
    if plotonfly
        subploti = 3;
        
        figure(1)
        ax(subploti)=subplot(Nsubs,1,subploti);
        hold('on')
        plot(Time(I),signalfiltdecontaminated); %set(gca,'xtick',[]); hold('off')
    end
    
    
    temp = signalfiltdecontaminated;
    tempout = buffer(temp,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
    temprms = 2*rms(tempout);
    
    if plotonfly
        ax(subploti)=subplot(Nsubs,1,subploti);
        hold('on')
        plot(Timerms,temprms,'r'); set(gca,'xtick',[]); hold('off')
    end
    
    % Throws an error - Timerms and temp different lengths - investigate?
    %     if plotonfly
    %         ax(subploti)=subplot(Nsubs,1,subploti);  hold('on');
    %         plot(Timerms,temp,'k'); set(gca,'xtick',[]); hold('off');
    %     end
    
    dT = Timerms(2)-Timerms(1);
    ECG_peak_i_rms = round(ECG_peak_i*dt/dT);
    
    lefttemplatetime = RRnormal*0.2; righttemplatetime = RRnormal*0.2;
    leftdi=round(1/dT*lefttemplatetime); rightdi=round(1/dT*righttemplatetime);
    polyorder=2;
    if plotonfly
        templateonly=1;
        [~,mediantemplate] = crosscontaminated(temprms,ECG_peak_i_rms,leftdi,rightdi,templateonly,polyorder);
    end
    [temprmsdecontaminated,mediantemplate] = crosscontaminated(temprms,ECG_peak_i_rms,leftdi,rightdi,0,polyorder);
    
    lefttemplatetime = RRnormal*0.33; righttemplatetime = RRnormal*0.0;
    leftdi=round(1/dT*lefttemplatetime); rightdi=round(1/dT*righttemplatetime);
    polyorder=2;
    if plotonfly
        templateonly=1;
        [~,mediantemplate] = crosscontaminated(temprmsdecontaminated,ECG_peak_i_rms,leftdi,rightdi,templateonly,polyorder);
    end
    [temprmsdecontaminated,mediantemplate] = crosscontaminated(temprmsdecontaminated,ECG_peak_i_rms,leftdi,rightdi,0,polyorder);
    
    lefttemplatetime = RRnormal*0.6; righttemplatetime = RRnormal*0.4;
    leftdi=round(1/dT*lefttemplatetime); rightdi=round(1/dT*righttemplatetime);
    polyorder=2;
    if plotonfly
        templateonly=1;
        [~,mediantemplate] = crosscontaminated(temprmsdecontaminated,ECG_peak_i_rms,leftdi,rightdi,templateonly,polyorder);
    end
    [temprmsdecontaminated,mediantemplate] = crosscontaminated(temprmsdecontaminated,ECG_peak_i_rms,leftdi,rightdi,0,polyorder);
    
    if plotonfly
        figure(1)
        ax(subploti)=subplot(Nsubs,1,subploti);  hold('on');
        plot(Timerms,temprmsdecontaminated,'b','linewidth',0.5); set(gca,'xtick',[]); hold('off')
    end
    
    EdiRMS_(:,nn) = interp1(Timerms,temprmsdecontaminated,Time(I),'linear','extrap');
end

%% Merge windows

Noverlap=find(Time_(:,2)>=Time_(end,1),1);
Noverlap=round(Noverlap*2/3);
removerows = round(Noverlap/2);
NoverlapRemaining = Noverlap - removerows;

N3 = size(EdiRMS_,1);
Foverlap2 = ones(N3,1);
temp = linspace(0,NoverlapRemaining,NoverlapRemaining)'/NoverlapRemaining;
Foverlap2(1:removerows) = 0;
Foverlap2((end-removerows+1):end) = 0;
Foverlap2((removerows+1):(removerows+NoverlapRemaining)) = temp;
Foverlap2((end-NoverlapRemaining-removerows+1):(end-removerows)) = flipud(temp);

SignalAllrs = 0*Time;

for nn=1:(size(EdiRMS_,2)-1)
    SignalAllrs(indarray(:,nn))=SignalAllrs(indarray(:,nn))+Foverlap2.*EdiRMS_(:,nn); %actual
end

%% Clear memory

clear EdiRMS_ EdiRMSlow_ Time_ EdiSum indrows noisewav Temp EdiRMS_t


%% Artifacts

windowwidth=1;
Foverlap=0.95;
tempout = buffer(SignalAllrs,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
tempsignal = nanmedian(tempout)';
tempout = buffer(Time,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
temptime = tempout(round(size(tempout,1)/2),:);
SignalMedian1 = interp1(temptime,tempsignal,Time,'linear','extrap');
clear signalfiltlow tempout temptime

thres1=3;
SignalAllrs_ = SignalAllrs;
SignalAllrs_(SignalAllrs>thres1*SignalMedian1)=thres1*SignalMedian1(SignalAllrs>thres1*SignalMedian1);

if 1
    mergeminarti=round(0.2/dt);
    Err = (SignalAllrs-SignalAllrs_)>0;
    I1=find(diff(Err)==1);
    I2=find(diff(Err)==-1);
    [I1_,I2_] = TidyStartEndEventList(I2,I1,N); %interdistance
    Warts_ = (I2_-I1_);
    rem = find(Warts_>mergeminarti);
    Warts_(rem)=[];
    I1_(rem)=[];
    I2_(rem)=[];
    for i=1:length(I1_)
        Err(I1_(i):I2_(i))=1;
    end
    
    I1=find(diff(Err)==1);
    I2=find(diff(Err)==-1);
    [I1,I2] = TidyStartEndEventList(I1,I2,N);
    Warts = (I2-I1);
    for i=1:length(I1)
        SignalAllrs_(I1(i):I2(i))=NaN;
    end
end

%% Interpolate through artifacts

I1=find(diff(Err)==1);
I2=find(diff(Err)==-1);
[I1,I2] = TidyStartEndEventList(I1,I2,N);
Warts = (I2-I1);
rem = find(Warts<mergeminarti);
Warts(rem)=[];
I1(rem)=[];
I2(rem)=[];

Err=0*Err; %replaced to keep just long artifacts
ploton=0;
if 1
    for i=1:length(I1)
        li1 = I1(i)-round(Warts(i));
        ri1 = I2(i)+round(Warts(i));
        temp = SignalAllrs_(li1:ri1);
        timetemp = Time(li1:ri1);
        temp2 = temp;
        timetemp2 = timetemp;
        isnans = isnan(temp)|isnan(timetemp);
        temp(isnans)=[];
        timetemp(isnans)=[];
        
        s=fitglm(timetemp,temp);
        ypred=predict(s,timetemp);
        I = ~isoutlier(ypred-temp);
        s=fitglm(timetemp(I),temp(I));
        ypred=predict(s,timetemp2);
        tempnhannw = tukeywin(length(li1:ri1),0.67);
        ypred2 = (tempnhannw).*ypred+(1-tempnhannw).*SignalAllrs(li1:ri1);
        
        SignalAllrs_(li1:ri1)=ypred2;
        Err(li1:ri1)=1;
        
        if ploton
            figure(2); clf(2);
            diplot=round(5/dt);
            Iplot = (li1-diplot):(ri1+diplot);
            plot(Time(Iplot),SignalAllrs(Iplot));  hold('on');
            %         plot(timetemp,temp);
            plot(timetemp2,ypred2,'k');
            hold('off');
            
            pause(0.1)
        end
    end
end

%% recalculate SignalMedian moving time average
%500ms moving time average

windowwidth=0.5;
Foverlap=0.95;
tempout = buffer(SignalAllrs_,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
tempsignal = nanmedian(tempout)';
tempout = buffer(Time,round(windowwidth/dt),round(Foverlap*windowwidth/dt),'nodelay');
temptime = tempout(round(size(tempout,1)/2),:);
SignalMedian = interp1(temptime,tempsignal,Time,'linear','extrap');
clear signalfiltlow tempout temptime

%% Plot
finalploton=1;
if finalploton
    figure(4); clf(4);
    set(gcf,'color',[1 1 1]);
    dsfactor=10;
    ax(1)=subplot(2,1,1);
    hx(1)=plot(downsample(Time,dsfactor),downsample(signalT,dsfactor),'color',[0.9 0.9 0.9]);
    hold('on')
    hx(2)=plot(downsample(Time,dsfactor),downsample(SignalAllrs,dsfactor),'color',[0.8 0.8 1]);
    hx(3)=plot(downsample(Time,dsfactor),downsample(SignalAllrs_,dsfactor),'color',[1 0.9 0.7]);
    hx(5)=plot(downsample(Time,dsfactor),downsample(Err*0.2,dsfactor),'color',[1 0.8 0.8]); %delete(hx)
    hx(4)=plot(downsample(Time,dsfactor),downsample(SignalMedian,dsfactor),'k'); %delete(hx)
    ylim([0 0.3]);
    box('off');
    %ylim([-5 5]);
    ax(2)=subplot(2,1,2);
    %plot(downsample(Time,dsfactor),downsample(Flow,dsfactor),'k'); %delete(hx)
    plot(time_ds,Vflow_ds,'k'); %delete(hx) %Flow already downsampled, keep orig if Spike, if Luo use downsampled (though sf 20 currently)
    box('off');
    linkaxes(ax,'x')
    %% page up
    if 0
        xlims=get(gca,'xlim'); xlims = xlims + diff(xlims); xlim(xlims);
    end
end
%%

EMGdi=interp1(Time,SignalMedian,time_ds);
EMGdi(EMGdi>100) = nan;
Err=interp1(Time,Err,time_ds);
Time125=time_ds;

% EMGdi=resample(SignalMedian,1,16);
% Err=resample(Err,1,16);
% Time125=resample(Time,1,16);

