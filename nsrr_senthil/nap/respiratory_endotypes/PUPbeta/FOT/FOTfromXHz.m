function [SigT,ChannelsFs,FOTinfo]=FOTfromXHz(SigT,ChannelsFs,ploton)
ChannelsList=SigT.Properties.VariableNames;
%FOT in spike:
%1. 0.03s smoothing on pO2.
%2. try to add procedure to correct for delays in CO2/O2.


fs2=128; %
compFreq = 8;
timeatend = 0; %rather than middle of window
clear Vol

delayvalue=0.0160; %157 cm
%Note: phase lag equivalent to 0.0160 s for a 157cm (long) Pmask tube.

%Pmask tube in Rm 16 is 92.5 cm (3/2015), equivalent to 0.0094 s. Note
%0.0094 is the EXACT delay value calculated in 1681N1pm data.

%delay down airway should be approx 0.15m/(300m/s)=0.0005 s;


%% Time arrays


time = SigT.Time;
N = length(time);
dt = time(2)-time(1);
fs = 1/dt;


%% Naming of variables

ch_P = find(strcmp(ChannelsList,'PmaskOrig'));
ch_V = find(strcmp(ChannelsList,'FlowOrig'));
if isempty(ch_P)
    ch_P = find(strcmp(ChannelsList,'Pmask'));
    ch_V = find(strcmp(ChannelsList,'Flow'));
end
Pmask = SigT{:,ch_P};
Vflow = SigT{:,ch_V};

try 
    Pepi = SigT.Pepi; 
end


%% examine phase shift between Pmask and Pepi (should be zero)
if exist('Pepi') && ~isempty(Pepi)
    filtOrd = 4;
    wcHigh = (compFreq*0.5)/(fs/2);
    wcLow = (compFreq*1.5)/(fs/2);
    [numHigh,denomHigh] = butter(filtOrd,wcHigh,'high');
    [numLow,denomLow] = butter(filtOrd,wcLow,'low');
    
    % bandpass filter the input signals
    pHigh1 = filtfilt(numHigh,denomHigh,Pmask);
    Pmask1_high = filtfilt(numLow,denomLow,pHigh1); %mask
    
    peHigh1 = filtfilt(numHigh,denomHigh,Pepi);
    Pepi1_high = filtfilt(numLow,denomLow,peHigh1); %mask
    
    nlags=round(1*fs/20); %range of lags = 1 sec either side.
    [XCF,Lags,Bounds] = crosscorr(Pmask1_high,Pepi1_high,nlags);
    Lags_us = -nlags:0.01:nlags;
    tempXCF = resample(XCF,100,1);
    templags = resample(Lags,100,1);
    [maxXCF_us,maxXCF_i]=max(tempXCF);
    
    delay_PepitoPmask=-templags(maxXCF_i)*dt; %%%%%%%%%% OUTPUT INFO 
    
    figure(130), plot(Lags*dt,XCF,Lags_us(maxXCF_i)*dt,maxXCF_us,'.');
    
    %XCF_us = interp1(Lags,XCF,Lags_us,'spline');
    %[maxXCF_us,maxXCF_i]=max(XCF_us);
    %delay_PepitoPmask=-Lags_us(maxXCF_i)*dt
    if ploton==2
        figure(130), plot(Lags*dt,XCF,Lags_us(maxXCF_i)*dt,maxXCF_us,'.');
        xlabel('Lag of expired Pepi with respect to Pmask (s)');
        ylabel('R value');
    end
end

%% Circshift Pmask
delay_correction_timeshift=0;
if delay_correction_timeshift
    Pmask = lagmatrix(Pmask,-round(delayvalue*fs))'; % i.e. just 1 element.
    actual_timeshift=round(delayvalue*fs)/fs;
    Pmask(isnan(Pmask))=0;
end

%% Check Spectra


df = 1/(time(end)-time(1));
[Pxx1,f] = pwelch(Vflow,N,0,N,fs);
[Pyy1,f] = pwelch(Pmask,N,0,N,fs);
if ploton==2
figure(2)
subplot(2,1,1); semilogy(f,Pxx1); xlim([0 50]);
subplot(2,1,2); semilogy(f,Pyy1); xlim([0 50]);
end
%% Raw data plot
if ploton==2
figure(2)
ax2(2)=subplot(4,1,1), plot(time,Vflow)
ax2(3)=subplot(4,1,3), plot(time,Pmask)
linkaxes(ax2,'x')
end

%% Interpolate to 128 Hz for simplicity (FOT signal of 8Hz)

dt2=1/fs2;
N2=floor(N*fs2/fs);
time2 = time(1):dt2:((N2-1)*dt2+time(1));
Vflow2 = interp1(time,Vflow,time2,'spline');
Pmask2 = interp1(time,Pmask,time2,'spline');
if 0
figure(200)
ax200(1)=subplot(2,1,1), plot(time,Vflow,time2,Vflow2)
ax200(2)=subplot(2,1,2), plot(time,Pmask,time2,Pmask2)
end

%% FOT analysis 
% define the bandpass filter (cascaded low and high pass)
filtOrd = 4;
wcHigh = (compFreq*0.5)/(fs2/2);
wcLow = (compFreq*1.5)/(fs2/2);
[numHigh,denomHigh] = butter(filtOrd,wcHigh,'high');
[numLow,denomLow] = butter(filtOrd,wcLow,'low');

% bandpass filter the input signals
qHigh = filtfilt(numHigh,denomHigh,Vflow2);
pHigh = filtfilt(numHigh,denomHigh,Pmask2);
Vflow2_high = filtfilt(numLow,denomLow,qHigh);        
Pmask2_high = filtfilt(numLow,denomLow,pHigh); %mask

% filter off the FOT from the flow/pressure signals; performed at each of the first 3 harmonics
filtOrd = 1;
bandwidthnotch=1;
[Bstop,Astop] = butter(filtOrd,(compFreq+bandwidthnotch*[-1 1])/(fs2/2)','stop');
Vflow2_notch = filtfilt(Bstop,Astop,Vflow2);
Pmask2_notch = filtfilt(Bstop,Astop,Pmask2);
[Bstop,Astop] = butter(filtOrd,(2*compFreq+bandwidthnotch*[-1 1])/(fs2/2)','stop');
Vflow2_notch = filtfilt(Bstop,Astop,Vflow2_notch);
Pmask2_notch = filtfilt(Bstop,Astop,Pmask2_notch);
[Bstop,Astop] = butter(filtOrd,(3*compFreq+bandwidthnotch*[-1 1])/(fs2/2)','stop');
Vflow2_notch = filtfilt(Bstop,Astop,Vflow2_notch);
Pmask2_notch = filtfilt(Bstop,Astop,Pmask2_notch);
if ploton==2
    figure(200);
    ax2(1)=subplot(2,1,1), plot(time2,Vflow2_high,time2,Vflow2_notch);
    ax2(2)=subplot(2,1,2), plot(time2,Pmask2_high,time2,Pmask2_notch);
    linkaxes(ax2,'x');
end
%% Calculate impedances (slow section)

NcyclesPerWindowYFast=1;
fftLen = NcyclesPerWindowYFast*fs2/compFreq; %1 cycle per window
stepSize = 1; %single sample step size!
numFFTs = 1 + floor((N2 - fftLen)/stepSize);
compIndex = compFreq*fftLen/fs2 + 1;
Pxy = zeros(numFFTs,1); 
Pxx = zeros(numFFTs,1);
Pyy = zeros(numFFTs,1);
timeFOT1= zeros(numFFTs,1);
timeFOT1(1)=time2(1);

%Pmask delay:
delaycorrection1 = delayvalue; %a positive value 'undoes' the delay [delay: exp(-j.omega.T)]
phasecorrection1 = 2*pi*delaycorrection1*compFreq;
phasecorrectionfactor1 = exp(j*phasecorrection1);

for k = 1:numFFTs
  begIndex = (k - 1)*stepSize + 1;
  endIndex = begIndex + fftLen - 1;
  if ~timeatend
    timeFOT1(k)=time2(round((begIndex+endIndex)/2));
  else
    timeFOT1(k)=time2(endIndex);  
  end
    
  Y = fft(Pmask2_high(begIndex:endIndex)); 
  X = fft(Vflow2_high(begIndex:endIndex))*phasecorrectionfactor1; %%% mult by delaycorrection1
  Pxy(k) = conj(X(compIndex)).*Y(compIndex); %*phasecorrectionfactor1
  Pxx(k) = conj(X(compIndex)).*X(compIndex);
  Pyy(k) = conj(Y(compIndex)).*Y(compIndex);
end
RealZfast = real(Pxy./Pxx); %previously means
ImagZfast = imag(Pxy./Pxx); %previously means
Zabsfast = abs(Pxy./Pxx);
Yabsfast = 1./abs(Pxy./Pxx);


%% Coh: longer time scale
NwavesForCoh=3;
cohwindowwidth=1/compFreq*NwavesForCoh; %denotes 3 FOT waves; %cohwindowwidth is seconds of data used
winLen = cohwindowwidth*fs2; %Nsamples per window @128 Hz
stepSize = 8; %steps by 8 samples, half a cycle 
numAvgs = 1 + floor((numFFTs - winLen)/stepSize);
Coh = zeros(numAvgs,1);
Pxx_ = zeros(numAvgs,1);
Pxy_ = zeros(numAvgs,1);
timeFOT2 = zeros(numAvgs,1);
for k = 1:numAvgs
    begIndex = (k - 1)*stepSize + 1;
    endIndex = begIndex + winLen - 1;
    if ~timeatend
       timeFOT2(k)=mean(timeFOT1((begIndex:endIndex)));
    else
       timeFOT2(k)=timeFOT1(endIndex);  
    end
    Coh(k) =    abs(mean(Pxy(begIndex:endIndex)))^2 / (mean(Pxx(begIndex:endIndex))) / mean(Pyy(begIndex:endIndex));
    Pxx_(k) =   median(Pxx(begIndex:endIndex)); %unused
    Pxy_(k) =   median(Pxy(begIndex:endIndex)); %unused
end

%RealZfast_ = real(Pxy_./Pxx_); %unused
%ImagZfast_ = imag(Pxy_./Pxx_); %unused
%Zabsfast_ = abs(Pxy_./Pxx_); %unused
Yabsslow = 1./abs(Pxy_./Pxx_); %unused


%% OUTPUTS
FlowOrig = Vflow;
PmaskOrig = Pmask;
FlowClean = interp1(time2,Vflow2_notch,time);
PmaskClean = interp1(time2,Pmask2_notch,time);
FlowHigh = interp1(time2,Vflow2_high,time);
PmaskHigh = interp1(time2,Pmask2_high,time);
PpmaskFOT = interp1(timeFOT1,abs(Pyy),time); %Y is Pres
PflowFOT = interp1(timeFOT1,abs(Pxx),time); %X is Flow
YFOT = interp1(timeFOT1,Yabsfast,time);
CohFOT = interp1(timeFOT2,Coh,time);

%% CleanY
%Fcoh = 1/(1-0.9)*(CohFOT-0.9); %(1->1, 0.9->0; slope = 1/(1-0.9), xint=0.9)
%Fcoh(Fcoh<0)=0;
FOTmincohthres = 0.5;

YabsfastCoh = YFOT;
YabsfastCoh(CohFOT<FOTmincohthres) = NaN;

%%Also noise
minPp_Flog10Pctile95 = 3;
thres = prctile(log10(PpmaskFOT),[95])-minPp_Flog10Pctile95; %100 times smaller than larger values are noise; ratio Y=Flow/Pres(~0) is meaningless
YabsfastCoh(log10(PpmaskFOT)<thres) = NaN;

minPp = 0.33; %100 times smaller than larger values are noise; ratio Y=Flow/Pres(~0) is meaningless
YabsfastCoh(log10(PpmaskFOT)<log10(minPp)) = NaN;

if 0
thres = prctile(log10(PflowFOT),[95])-2; %100 times smaller than larger values are noise.
YabsfastCoh(log10(PflowFOT)<thres) = NaN;

thres = prctile(log10(PpmaskFOT)+log10(PflowFOT),[95])-3; %10 times smaller than larger values are noise; implies leak or dropout.
YabsfastCoh((log10(PpmaskFOT)+log10(PflowFOT))<thres) = NaN;
end

% upper outlier, remove completely
thres = 5*diff(prctile(YabsfastCoh,[50 90])) + prctile(YabsfastCoh,[90]);
YabsfastCoh(YabsfastCoh>thres) = NaN;

if 0 %slow and may remove too many segments of good data
    Out = RemoveShortSegments(~isnan(YabsfastCoh),0.5,1/fs,0); %keep good segments no shorter than 0.5 s
    YabsfastCoh(Out==0)=NaN;
end

Fisnan = sum(isnan(YabsfastCoh))/length(YabsfastCoh);
disp(['Fisnan = ' num2str(Fisnan,2)]);

%% Details
FOTinfo.FOTmincohthres=FOTmincohthres;
FOTinfo.Fisnan=Fisnan;
FOTinfo.minPp_Flog10Pctile95=minPp_Flog10Pctile95;
FOTinfo.minPp=minPp;
FOTinfo.NwavesForCoh=NwavesForCoh;
FOTinfo.NcyclesPerWindowYFast=NcyclesPerWindowYFast;

%% Main Plots if desired
if ploton==1
    figure(10); clf(10);
    
    ax3(1)=subplot(5,1,1);
    plot(time,Vflow,time,2+Pmask/nanstd(PmaskClean)*nanstd(FlowClean)); ylabel('Flow/Pmask clean');
    
    hold on
    plot(time,FlowClean,'k',time,2+PmaskClean/nanstd(PmaskClean)*nanstd(FlowClean),'k'); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    
    ax3(2)=subplot(5,1,2); plot(time,FlowHigh,time,0.5+PmaskHigh/nanstd(PmaskHigh)*nanstd(FlowHigh)); ylabel('Flow/Pmask HF'); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    ax3(3)=subplot(5,1,3); semilogy(time,[(abs(PflowFOT)) (abs(PpmaskFOT))]); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    ylabel('Pxx, Pyy');
    ax3(4)=subplot(5,1,4); plot(time,YFOT); ylim([0 1.5]); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    hold on
    plot(time,[YabsfastCoh],'k','linewidth',2.5);
    hold off
    ylabel('|Y=V/P|');
    %ax3(3)=subplot(4,1,3); plot(time2,volume1,time2,detrend(volume1));
    ax3(5)=subplot(5,1,5); plot(time,CohFOT); set(gca,'box','off');
    ylabel('Coh');
    linkaxes(ax3,'x')
end
dt_FOT2 = median(diff(timeFOT2));

%% Add to Matrix
%SigT = array2table(DataEventHypnog_Mat);
%DataEventHypnog_MatT.Properties.VariableNames = ChannelsList;

%overwrite
SigT.Flow = FlowClean;
SigT.Pmask = PmaskClean;

%add channels to end
SigT.YFOT = YFOT;
SigT.YFOTclean = YabsfastCoh; 
SigT.CohFOT = CohFOT;
SigT.PflowFOT = PflowFOT;
SigT.PpmaskFOT = PpmaskFOT;
SigT.FlowOrig = FlowOrig;
SigT.PmaskOrig = PmaskOrig;

%% Replace
ChannelsList = SigT.Properties.VariableNames;
%DataEventHypnog_Mat = DataEventHypnog_MatT{:,:};
ChannelsFs(end:length(ChannelsList))=ChannelsFs(end);






