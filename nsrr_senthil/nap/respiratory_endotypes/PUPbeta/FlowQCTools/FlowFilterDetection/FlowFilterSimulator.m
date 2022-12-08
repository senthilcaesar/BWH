function [log10Prel,log10PrelH]=FlowFilterSimulator(flow,Time,LFmethod,LFinputtest,HFmethod,HFinputtest,ploton,plotdims)

if ~exist('plotdims')
    plotdims = [8 6];
end

%load the _XHz data then run this cell
%ploton=0;
predictlowpass=0;
%LFinputtest = 5;
freqsweep = [1:15];
%LFmethod=4

%
% flow = DataEventHypnog_Mat(:,2);
dt = Time(2,1)-Time(1,1);
%    

if ploton
figure(10);
set(gcf,'color',[1 1 1]);
%clf(10);
end
%
predrange = [1:15]; %just for the prediction model

   
    Ierr = isnan(flow)| flow==Inf | flow==-Inf;
    flow(Ierr)=nanmedian(flow);
    
    
    %or skip the above and provide "flow" and "dt" then proceed below
    %note: need to load "FilterSmoothCutoffDetector" to get mdlSmooth
switch LFmethod
    case 0
        flowfiltered = flow;
    case 1  %simulate a 1st order lowpass filter
        %clear log10Prel
        filter_HFcutoff_butter1 = LFinputtest;
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter1]/(1/dt/2),'low');
        flowfiltered = filter(B_butter0,A_butter0,flow);
    case 2  %simulate a 1st order lowpass filter, order 2
        %clear log10Prel
        filter_HFcutoff_butter1 = LFinputtest;
        filter_order0 = 2;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter1]/(1/dt/2),'low');
        flowfiltered = filter(B_butter0,A_butter0,flow);
    case 3  %simulate downsampling, and reverse upsampling
        downsamplefs=2*LFinputtest;
        TotalTime = (length(flow)-1)*(dt);
        Time_ = [0:dt:TotalTime]';
        TimeDS=[0:(1/downsamplefs):TotalTime]';
        FlowDS = interp1(Time_,flow,TimeDS,'linear','extrap'); % downsample
        flowfiltered = interp1(TimeDS,FlowDS,Time_,'linear','extrap');  % upsample back to original length.
        clear Time_ TimeDS FlowDS
    case 4  %simulate downsampling, and reverse upsampling
        downsamplefs=2*LFinputtest;
        TotalTime = (length(flow)-1)*(dt);
        Time_ = [0:dt:TotalTime]';
        TimeDS=[0:(1/downsamplefs):TotalTime]';
        FlowDS = interp1(Time_,flow,TimeDS,'pchip','extrap'); % downsample
        flowfiltered = interp1(TimeDS,FlowDS,Time_,'pchip','extrap');  % upsample back to original length.
        clear Time_ TimeDS FlowDS 
end

switch HFmethod
    case 0
        %do nothing
    case 1  %simulate a 1st order lowpass filter
        if HFinputtest~=0 %simulate a highpass filter, order 1
            %HFinputtest = 0.1;
            %clear log10Prel
            filter_HFcutoff_butter1 = HFinputtest;
            filter_order0 = 1;
            [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter1]/(1/dt/2),'high');
            flowfiltered = filter(B_butter0,A_butter0,flowfiltered);
        end
    case 2
        if HFinputtest~=0 %simulate a highpass filter, order 1
            %HFinputtest = 0.1;
            %clear log10Prel
            filter_HFcutoff_butter1 = HFinputtest;
            filter_order0 = 2;
            [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter1]/(1/dt/2),'high');
            flowfiltered = filter(B_butter0,A_butter0,flowfiltered);
        end
end


    if 0 %simulate clipping
        clipprctiles = 0;
        %clear log10Prel
        flowfiltered(flowfiltered>prctile(flowfiltered,100-clipprctiles))=prctile(flowfiltered,100-clipprctiles);
        flowfiltered(flowfiltered<prctile(flowfiltered,clipprctiles))=prctile(flowfiltered,clipprctiles);
    end


    freqsweep = freqsweep(predrange);
    nfft = 200*round(1/dt); %200 seconds per fft
    [Flowpsd,F] = pwelch(flowfiltered-mean(flowfiltered),hann(nfft),round(nfft/2),nfft,1/dt);
    df = F(2)-F(1);
    Fmax = 30;
    Fnyquist = 1/dt/2;
    I = F>Fmax;
    Flowpsd(I)=[];
    F(I)=[];
    
    Fl = 0.1; Fr = 12.5; Fli = round(Fl/df) + 1; Fri = round(Fr/df) ; %not inclusive
    % DLM. Fri was > length Flowpsd, and again lower down, so here we set Fri = length(Flowpsd)
    if Fri>length(Flowpsd); Fri=length(Flowpsd); end % DLM added this line for when Fri > length(Flowpsd)
    P1to12 = sum(Flowpsd(Fli:Fri))*df;
    
    Fr1=1;
    Fr1i= round(Fr1/df);
    P1 = sum(Flowpsd(Fli:Fr1i))*df;
    P1sd = mean(Flowpsd(Fli:Fr1i));
    
    
    if ploton
        figure(10);% clf(10);
        fig = gcf;
        fig.Units = 'Inches';
        currentPos = fig.Position;
        fig.Position = [currentPos(1), currentPos(2), plotdims(1), plotdims(2)];
        %plot fft:
        Nsubplots=7;
        subplot(Nsubplots,3,3*[1:Nsubplots]); loglog(F,Flowpsd/P1sd); 
        xlim([0.01 30]);
        ylim([10^-10 10^2]);
        grid on
        
        if 0
        hold on
        end
        
        set(gca,'box','off'); xlabel('Frequency, Hz');  ylabel('Flow Spectral Density')
        
        %plot example1:
        durationeg = 120;
        Fstart = [0.125:0.125:0.25+(Nsubplots-1)*0.125];
        for i=1:Nsubplots
            ileft = round((size(flow,1)*Fstart(i)));
            irange=ileft:ileft+round(durationeg/dt);
            tempx = Time(irange);
            temp=flowfiltered(irange);
            temp = (temp-nanmean(temp))/nanstd(temp);
            subplot(Nsubplots,3,(i-1)*3+[1 2]); 
            plot(tempx,temp); set(gca,'box','off','ylim',[-4 4],'xlim',[min(tempx) max(tempx)],'xtick',[],'xcolor',[1 1 1],'ytick',[],'ycolor',[1 1 1]);
            if i==1
               h=ylabel('Example signals','fontweight','normal','color',[0 0 0],'position',[9.5170e+04 -15.8327 -1.0000]);
            end
        end
        
    end
    
    
    range_ = freqsweep';
    clear Prel Prel2 log10PrelTest
    for i=1:size(range_,1)
        Fl = range_(i); Fr = range_(i)+1; Fli = round(Fl/df) + 1; Fri = round(Fr/df) ; %not inclusive
        if Fri>length(Flowpsd); Fri=length(Flowpsd); end % DLM added this line for when Fri > length(Flowpsd)
        Prel(i,1) = sum(Flowpsd(Fli:Fri))*df/P1to12;
        Prel2(i,1) = sum(Flowpsd(Fli:Fri))*df/P1;
    end
    log10PrelTest = log10(Prel)';
    log10Prel = log10(Prel2)';
    
    clear PrelH
    rangeH_ = [0.005 0.01:0.01:0.1]';
    for i=1:size(rangeH_)-1
        Fl = rangeH_(i); Fr = rangeH_(i+1); Fli = round(Fl/df) ; Fri = round(Fr/df) -1 ; %not inclusive
        if Fri>length(Flowpsd); Fri=length(Flowpsd); end % DLM added this line for when Fri > length(Flowpsd)
        PrelH(i,1) = sum(Flowpsd(Fli:Fri))*df/P1;
    end
    log10PrelTest = log10(Prel)';
    log10PrelH = log10(PrelH)';
    
    %
    % SDratiodiagnosis = interp1(freqsweep,log10Prel,filter_HFcutoff_butter1)
    if predictlowpass
    [FSmoothEst,FSmoothEst95CI]=predict(mdlSmooth,log10PrelTest);
    
    disp(['FSmoothEst = ' num2str(FSmoothEst)])
    end
    
    