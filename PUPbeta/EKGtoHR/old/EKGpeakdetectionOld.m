function [ECG_peak_i,HR] = EKGpeakdetection(EKG,Time,dt,BreakEarly,plotfig)

RemoveNoisyQRSforOutput=1;
%settings for EDIfilter
% RemoveNoisyQRSforOutput=1; %irrelevant if BreakEarly
% BreakEarly=1;
% plotfig=1; 

%settings for normal HR use 
% BreakEarly=0;
% plotfig=0; 

EKG = EKG(:);
Time = Time(:);
HR=NaN;


%% Find RR
fs = round(1/dt);
N = length(Time);
%Filter EKG to aid peak identification. 
filter_LFcutoff_butter = 10; %highpass, DC remove
filter_HFcutoff_butter = 45; %lowpass, smooth

if filter_HFcutoff_butter>fs/2
    filter_HFcutoff_butter=fs/2-0.01;
end

filter_order = 2; 
[B_butter,A_butter] = butter(filter_order,[filter_LFcutoff_butter filter_HFcutoff_butter]/(fs/2));
ECG_filtered = filtfilt(B_butter,A_butter,EKG);

%gently filter raw EKG.
filter_LFcutoff_butter = 1; %highpass, DC remove
filter_HFcutoff_butter = min([70 0.999*fs/2]); %lowpass, smooth
filter_order = 2; 
[B_butter,A_butter] = butter(filter_order,[filter_LFcutoff_butter filter_HFcutoff_butter]/(fs/2));
EKG_ = filtfilt(B_butter,A_butter,EKG);

clear EKG

% ECG_flip_on=1;
% FlipECG = ((mean(ECG_filtered)-prctile(ECG_filtered,0.5))>(prctile(ECG_filtered,99.5)-mean(ECG_filtered)));
% 
    Fs = 1/(Time(2)-Time(1));
    
    % the windowduration and secslide values were hardcoded here
    % they did not use the values from the options worksheet
    secslide=120;
    WindowDuration=180;
    
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    %FlipEKGvotes = NaN*ones(numwind,1);
    FlipEKGmag = NaN*ones(numwind,1);
    for winNum=1:numwind
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        if ri>length(Time)
            ri=length(Time);
        end
        I=round((li:ri));
        %FlipEKGvotes(winNum) = ((mean(ECG_filtered(I))-prctile(ECG_filtered(I),0.5))>(prctile(ECG_filtered(I),99.5)-mean(ECG_filtered(I))));
        FlipEKGmag(winNum) = ((mean(ECG_filtered(I))-prctile(ECG_filtered(I),0.5))-(prctile(ECG_filtered(I),99.5)-mean(ECG_filtered(I))));
    end
    FlipECG = nansum(FlipEKGmag)>0;
    
    
        if FlipECG
            ECG_filtered=-ECG_filtered;
            EKG_=-EKG_;
            disp('flipped EKG');
        end

       % plot(Time,ECG_filtered)
    
%% window-by-window thresholding
ECG_filtered_threshold_fraction=0.6;
if BreakEarly
    ECG_filtered_threshold_fraction=0.6;
end

PrctileLower=50;
PrctileUpper=99;

if Time(end)-Time(1)>1200 %if shorter window then don't use time-varying threshold
    Fs = 1/(Time(2)-Time(1));
    
    % the windowduration and secslide values were hardcoded here
    % they did not use the values from the options worksheet
    secslide=1;
    WindowDuration=6;
    
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    ECG_filtered_threshold_ = NaN*ones(numwind,1);
    ECG_filtered_threshold_t = NaN*ones(numwind,1);
    for winNum=1:numwind
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        if ri>length(Time)
            ri=length(Time);
        end
        I=round((li:ri));
        ECG_filtered_threshold_(winNum) = ECG_filtered_threshold_fraction*[prctile(ECG_filtered(I),PrctileUpper)-prctile(ECG_filtered(I),PrctileLower)]+prctile(ECG_filtered(I),PrctileLower);
        ECG_filtered_threshold_t(winNum) = nanmean(Time(I));
    end
    Itemp = isnan(ECG_filtered_threshold_);
    ECG_filtered_threshold_t(Itemp)=[];
    ECG_filtered_threshold_(Itemp)=[];
    %figure(),plot(ECG_filtered_threshold_t,ECG_filtered_threshold_,'.')
    ECG_filtered_threshold__ = interp1(ECG_filtered_threshold_t,ECG_filtered_threshold_,Time,'nearest','extrap');
else %single value threshold
    ECG_filtered_threshold__ = Time*0 + ECG_filtered_threshold_fraction*[prctile(ECG_filtered,99)-prctile(ECG_filtered,1)]+prctile(ECG_filtered,1);
end

%% Peak Detect, Initial

count=1;
clear ECG_peak ECG_peak_i ECG_peak_t;
for i=1:N %faster than findpeaks function
    if (i>2)&&(i<N)
        if ((ECG_filtered(i)>ECG_filtered_threshold__(i)&&(ECG_filtered(i)>ECG_filtered(i-1))&&(ECG_filtered(i)>ECG_filtered(i+1)) || (ECG_filtered(i)>ECG_filtered(i-2))&&(ECG_filtered(i)>ECG_filtered(i+1))&&(ECG_filtered(i)==ECG_filtered(i-1))))
            ECG_peak_i(count)=i;
            count=count+1;
        end       
    end
end
%Store data
ECG_peak_i_step1=ECG_peak_i;
%ECG_peak_i(EKG(ECG_peak_i)<0)=[]; 

if 0
RR = [diff(Time(ECG_peak_i_step1))' NaN];
HR = zoh(Time(ECG_peak_i_step1),60./RR',Time); %very fast version
plot(Time,HR)
end

%% Remove array elements (RAE) [fast]
    %Removes EKG peaks that are too close to each other
    %The smaller of the two peaks are removed. 
  
%RAE Setup
ECG_peak_t_threshold_percentofmedian = [0.05 15 20 25 30 33 33]; %specify

Nrejected = [];
for X=1:length(ECG_peak_t_threshold_percentofmedian) %runs several times. 2 cases usually enough.
RRa = [diff(ECG_peak_i)*dt NaN];
RRb = [NaN RRa(1:end-1)];
RRc = [NaN NaN RRa(1:end-2)];
data = [RRa;RRb;RRc];

covRR = std(data)./mean(data); %local cov

%find normal RR
upper = 4;
lower = 0.25; %hard lower value
covmax = 1/6; %just for median estimate

if BreakEarly
    lower=0.2;
    if X==4
        break
    end
end

ECG_peak_t_delta_median = nanmedian(RRa(RRa<upper&RRa>lower));
if X==length(ECG_peak_t_threshold_percentofmedian)
    ECG_peak_t_delta_median = nanmedian(RRa(RRa<upper&RRa>lower&covRR<covmax));
end
ECG_peak_t_threshold_min = ECG_peak_t_delta_median*ECG_peak_t_threshold_percentofmedian(X)/100;
if ECG_peak_t_threshold_min<lower
    ECG_peak_t_threshold_min=lower;
end
%find heights
temp = [ECG_filtered(ECG_peak_i(1:end))];
EKGheight = [temp [temp(2:end);NaN]];
leftpeaksmallest = EKGheight(:,1)<EKGheight(:,2);

%find closely-spaced EKG to reject

Inext = RRa<ECG_peak_t_threshold_min&(leftpeaksmallest'==0);
    Inext = [NaN Inext(1:end-1)];
Icurrent = RRa<ECG_peak_t_threshold_min&(leftpeaksmallest'==1);
I = Inext==1 | Icurrent==1;

Nrejected(X)=sum(I==1);
ECG_peak_i(I) = [];
end

ECG_peak_t_threshold_percentofmedian=ECG_peak_t_threshold_percentofmedian(end);

%%
if BreakEarly
    return
end

%% Remove ECG_peak_i if removing beat greatly reduces local cov(RR)
COVratioThres = 5; %2

Nloops=2; %300

N_removed = [];
for i=1:Nloops
    RRa = [diff(ECG_peak_i)*dt NaN];
    RRb = [NaN RRa(1:end-1)];
    RRc = [NaN NaN RRa(1:end-2)];
    RRz = [RRa(2:end) NaN];
    data = [RRz;RRa;RRb;RRc]; %5 beats (4 RRs), centered
    covRR = nanstd(data)./nanmean(data); %local cov
    RRab = RRa + RRb;
    dataExcl = [RRz;RRab;RRc]; %5 beats (4 RRs), centered
    covRRExcl = nanstd(dataExcl)./nanmean(dataExcl); %local cov
    sumisnans = sum(isnan(dataExcl));
    COVratio = covRR./covRRExcl;
    criteria = COVratio>COVratioThres & sumisnans<=1 & max(data)<upper;
    N_removed(i)=sum(criteria==1);
    ECG_peak_i(criteria==1)=[];
    if i>1 && N_removed(i)==0 && N_removed(i-1)==0
        break
    end
end

N_removed = [];
for i=1:Nloops
    RRa = [diff(ECG_peak_i)*dt NaN];
    RRb = [NaN RRa(1:end-1)];
    RRc = [NaN NaN RRa(1:end-2)];
    RRd = [NaN NaN NaN RRa(1:end-3)];
    data = [RRa;RRb;RRc;RRd]; %5 beats (4 RRs), centered
    covRR = nanstd(data)./nanmean(data); %local cov
    RRab = RRa + RRb;
    dataExcl = [RRab;RRc;RRd]; %5 beats (4 RRs), centered
    covRRExcl = nanstd(dataExcl)./nanmean(dataExcl); %local cov
    sumisnans = sum(isnan(dataExcl));
    COVratio = covRR./covRRExcl;
    criteria = COVratio>COVratioThres & sumisnans<=1 & max(data)<upper;
    N_removed(i)=sum(criteria==1);
    ECG_peak_i(criteria==1)=[];
    if i>1 && N_removed(i)==0 && N_removed(i-1)==0
        break
    end
end

N_removed = [];
for i=1:Nloops
    RRa = [diff(ECG_peak_i)*dt NaN];
    RRb = [NaN RRa(1:end-1)];
    RRy = [RRa(3:end) NaN NaN];
    RRz = [RRa(2:end) NaN];
    data = [RRy;RRz;RRa;RRb]; %5 beats (4 RRs), centered
    covRR = nanstd(data)./nanmean(data); %local cov
    RRab = RRa + RRb;
    dataExcl = [RRy;RRz;RRab]; %5 beats (4 RRs), centered
    covRRExcl = nanstd(dataExcl)./nanmean(dataExcl); %local cov
    sumisnans = sum(isnan(dataExcl));
    COVratio = covRR./covRRExcl;
    criteria = COVratio>COVratioThres & sumisnans<=1 & max(data)<upper;
    N_removed(i)=sum(criteria==1);
    ECG_peak_i(criteria==1)=[];
    if i>1 && N_removed(i)==0 && N_removed(i-1)==0
        break
    end
end

%% Detect single missing ECG_peak_i 
%backup = ECG_peak_i;
N_to_add2 = [];
N_to_add3 = [];
newbeats = 0*ECG_peak_i;
for i=1:1
    RRa = [diff(ECG_peak_i)*dt NaN];
    RRb = [NaN RRa(1:end-1)];
    RRc = [NaN NaN RRa(1:end-2)];
    RRz = [RRa(2:end) NaN];
    data = [RRz;RRa;RRb;RRc]; %5 beats (4 RRs), centered
    covRR = nanstd(data)./nanmean(data); %local cov
    
    dataExcl = [RRz;RRa/2;RRa/2;RRb;RRc]; %5 beats (4 RRs), centered
    covRRExcl = nanstd(dataExcl)./nanmean(dataExcl); %local cov
    
    dataExcl3 = [RRz;RRa/3;RRa/3;RRa/3;RRb;RRc]; %5 beats (4 RRs), centered
    covRRExcl3 = nanstd(dataExcl3)./nanmean(dataExcl3); %local cov
    
    sumisnans = sum(isnan(data));
    
    COVratio = covRR./covRRExcl;
    COVratio3 = covRR./covRRExcl3;
    
    criteria = COVratio>COVratioThres & COVratio>COVratio3 & sumisnans<=1 & max(data)<upper & RRa>1.1*ECG_peak_t_delta_median & RRa<2.9*ECG_peak_t_delta_median;
    
    N_to_add2(i)=sum(criteria==1);
    Newdata_ECG_peak_i = ECG_peak_i(criteria) + round((RRa(criteria)/2)/dt);
    newbeats = [newbeats 1+0*Newdata_ECG_peak_i];
    ECG_peak_i = [ECG_peak_i Newdata_ECG_peak_i];
    
    criteria3 = COVratio3>COVratioThres & COVratio3>COVratio & sumisnans<=1 & RRa>1.9*ECG_peak_t_delta_median & RRa<3.9*ECG_peak_t_delta_median;
    
    N_to_add3(i)=sum(criteria3==1);
    Newdata_ECG_peak_i = [ECG_peak_i(criteria3) + 1*round((RRa(criteria3)/3)/dt) ...
                          ECG_peak_i(criteria3) + 2*round((RRa(criteria3)/3)/dt)];
    
    newbeats = [newbeats 1+0*Newdata_ECG_peak_i];
    ECG_peak_i = [ECG_peak_i Newdata_ECG_peak_i];
    
    temp = sortrows([ECG_peak_i(:) newbeats(:)]);
    ECG_peak_i = temp(:,1)';
    newbeats = temp(:,2)';
end


%% Detect Noise based on very long RR OR highly variable RR
%do not run recursively as edge effects destroy good data. 
covRRthres=0.6;
maxRRthres=upper;

RRa = [diff(ECG_peak_i)*dt NaN];
RRb = [NaN RRa(1:end-1)];
RRc = [NaN NaN RRa(1:end-2)];
RRz = [RRa(2:end) NaN];
data = [RRz;RRa;RRb;RRc]; %5 beats (4 RRs), centered
covRR = nanstd(data)./nanmean(data); %local cov


minabsdiff = min(abs(diff(data)));
temp = diff(data);

covRRseries = zoh(Time(ECG_peak_i),covRR',Time); %very fast version
RRseries = zoh(Time(ECG_peak_i),RRa',Time); %very fast version
Peakheights = zoh(Time(ECG_peak_i),ECG_filtered(ECG_peak_i),Time); %very fast version
minabsdiffseries = zoh(Time(ECG_peak_i),minabsdiff',Time); %very fast version

% figure()
% plot(Time,RRseries)

Noise = 0+0*Time;
    secslide=5; %was 5
    Fs = 1/(Time(2)-Time(1));
    WindowDuration=15; %was 15
    
if Time(end)-Time(1)>WindowDuration %if shorter window then don't use time-varying threshold
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    %SpikinessW = zeros(1,numwind);
    %timetemp = zeros(1,numwind);
    Minabsdiffseries  = 0*Time;
    Spikiness = 0*Time;
    for winNum=1:numwind
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        if ri>length(Time)
            ri=length(Time);
        end
        I=round((li:ri));
        %timetemp(winNum)=mean(Time(I));
%         if winNum==78
%             1
%         end
        try
            Spikiness(I) = nanmean(Peakheights(I))/prctile((ECG_filtered(I)),80)<6;
        catch me
            Spikiness(I) = 1; %done because I may extend too far past peakheights data
            Noise(I)=1;
        end
        try
            Minabsdiffseries(I) = nanmean(minabsdiffseries(I))>0.5;
        catch me
            Minabsdiffseries(I) = 1; %done because I may extend too far past minabsdiffseries data
            Noise(I)=1;
        end
        try
            if nanmean(covRRseries(I))>covRRthres||(max(RRseries(I))>maxRRthres)%||(sum(isnan(RRseries(I)))>0)
                Noise(I)=1;
            end
        catch me
            Noise(I)=1;
        end
    end
    %Spikiness = interp1(timetemp,SpikinessW,Time,'nearest','extrap');
else 
    %not coded for shorter windows yet
end

if 0
    figure(44)
    plot(Time,[Noise,Spikiness+1.1,Minabsdiffseries+2.2])
end

Noise = Noise==1|Spikiness==1;%|Minabsdiffseries==1;
clear minabsdiffseries Minabsdiffseries

%% NearPeak

NotNearQRSpeak = 1+0*Time;
NotNearQRSpeakithres = round(0.05*Fs);

for i=1:length(ECG_peak_i)
    li=round(ECG_peak_i(i)-NotNearQRSpeakithres);
    ri=round(ECG_peak_i(i)+NotNearQRSpeakithres);
    if li<1
        li=1;
    end
    if ri>length(Time)
        ri=length(Time);
    end
    NotNearQRSpeak(li:ri)=0;
end

%% FLocalTimeAboveThres

FLocalTimeAboveThres = 1+0*ECG_peak_i';

minsearchi = round(Fs*(ECG_peak_t_threshold_min-0.05));
maxsearchi = round(Fs*(ECG_peak_t_delta_median));

for i=1:length(ECG_peak_i)
    diplus=minsearchi;
    diminus=minsearchi;
    if i<length(ECG_peak_i)
        diplus = (ECG_peak_i(i+1)-ECG_peak_i(i))/2;
    end
    if diplus>maxsearchi
        diplus=maxsearchi;
    end
    if i>1
        diminus = (ECG_peak_i(i)-ECG_peak_i(i-1))/2;
    end
    if diminus>maxsearchi
        diminus=maxsearchi;
    end
    li=round(ECG_peak_i(i)-diminus);
    ri=round(ECG_peak_i(i)+diplus);
    if li<1
        li=1;
    end
    if ri>length(Time)
        ri=length(Time);
    end
    I=li:ri;
    FLocalTimeAboveThres(i)=sum(ECG_filtered(I)>ECG_filtered_threshold__(I)&NotNearQRSpeak(I)>0)/sum(NotNearQRSpeak(I)>0);
end

%% RMS
    secslide=nanmedian(RRa)*0.5;
    WindowDuration=secslide*2;
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    %FlipEKGvotes = NaN*ones(numwind,1);
    %RMS = NaN*ones(numwind,1);
    RMSraw = NaN*ones(numwind,1);
    %RMSoverthres = NaN*ones(numwind,1);
    RMS_t = NaN*ones(numwind,1);
    for winNum=1:numwind
        li=round((winNum-1)*secslide*Fs+1);
        ri=round((winNum-1)*Fs*secslide+WindowDuration*Fs);
        if li<1
            li=1;
        end
        if ri>length(Time)
            ri=length(Time);
        end
        I=li:ri;
        
        RMS_t(winNum) = mean(Time(I));
        %RMS(winNum) = rms(ECG_filtered(I).*NotNearQRSpeak(I));
        RMSraw(winNum) = rms(EKG_(I).*NotNearQRSpeak(I));
        I2 = ECG_filtered(I)>ECG_filtered_threshold__(I);
        temp = ECG_filtered(I(I2)).*NotNearQRSpeak(I(I2))-ECG_filtered_threshold__(I(I2)).*NotNearQRSpeak(I(I2));
        %RMSoverthres(winNum) = rms(temp);
    end
     %I=find(isnan(RMS_t)|isnan(RMS)|isnan(RMSraw));
     I=find(isnan(RMSraw));
     %RMS(I)=[];
     RMSraw(I)=[];
    %RMSoverthres(I)=[];
     RMS_t(I)=[];
    %RMS__ = interp1(RMS_t,RMS,Time,'nearest','extrap');
    RMSraw__ = interp1(RMS_t,RMSraw,Time,'nearest','extrap');
    %RMSoverthres__ = interp1(RMS_t,RMSoverthres,Time,'nearest','extrap');
    
%% Make short lengths of good data => noise
if 0
minlength=10;

I1=find(diff(1-Noise)==1); 
I2=find(diff(1-Noise)==-1); 
[I1,I2] = TidyStartEndEventList(I1,I2,length(Noise));
lengths = (I2-I1)*dt;

for i=1:length(I1)
    if lengths(i)<minlength
        Noise(I1(i):I2(i))=1;
    end
end
end

%% Move index to peak of unfiltered EKG
delta=max([2 round(Fs*0.05)]);

for i=1:length(ECG_peak_i)
    rangei = ECG_peak_i(i)-delta:ECG_peak_i(i)+delta;
    
    rangei(rangei<1)=[];
    rangei(rangei>length(EKG_))=[];
    
    [~,temp] = max(EKG_(rangei));
    ECG_peak_i(i)=rangei(temp);
end

%% Calculate RR (quadratic fit to peak to improve RR resolution from 1/fs)

ECG_peak_i = unique(ECG_peak_i);

if 1
    ECG_peak_quadfit_t = Time(ECG_peak_i);
    EditedQuadraticFit = 0*ECG_peak_quadfit_t;
    for i=1:length(ECG_peak_i)
        if ECG_peak_i(i)<2||ECG_peak_i(i)>(N-1)
            continue
        end
        I = ECG_peak_i(i)-1:ECG_peak_i(i)+1;
        if EKG_(I(2))>EKG_(I(1))&&EKG_(I(2))>EKG_(I(3))
            ECG_peak_quadfit_t(i)=PeakFitQuadratic(Time(I),EKG_(I));
            EditedQuadraticFit(i)=1;
        end
    end
    RR = [diff(ECG_peak_quadfit_t)' NaN];
else
    RR = [diff(Time(ECG_peak_i))' NaN];
end

HR = zoh(Time(ECG_peak_i),60./RR',Time); %very fast version

RMSRatioThres=0.75;
EKGpeakheightRatio = [ECG_filtered(ECG_peak_i)./RMSraw__(ECG_peak_i)];


NoiseRMS = zoh(Time(ECG_peak_i),EKGpeakheightRatio,Time)<RMSRatioThres; %very fast version
NoiseRMSleft = zoh(Time(ECG_peak_i(1:end-1)),EKGpeakheightRatio(2:(end)),Time)<RMSRatioThres; %very fast version

FLocalTimeAboveThresThres=0.005;
NoiseFlocal = zoh(Time(ECG_peak_i),FLocalTimeAboveThres,Time)>FLocalTimeAboveThresThres; %very fast version
NoiseFlocalleft = zoh(Time(ECG_peak_i(1:(end-1))),FLocalTimeAboveThres(2:(end)),Time)>FLocalTimeAboveThresThres; %very fast version
NoiseFlocalleft2 = zoh(Time(ECG_peak_i(1:(end-2))),FLocalTimeAboveThres(3:(end-1)),Time)>FLocalTimeAboveThresThres; %very fast version

NoiseRMS = NoiseRMS|NoiseRMSleft|NoiseFlocalleft2;
NoiseFlocal = NoiseFlocal|NoiseFlocalleft;
if 1 %remove Noise
    %HR(Noise==1)=NaN;
    HR(NoiseRMS==1)=NaN;
    HR(NoiseFlocal==1)=NaN;
    HR(HR<0)=NaN;
    HR(HR>4/ECG_peak_t_delta_median*60)=NaN;
end
if 1 %remove segments with lots of noise
   In = 1*isnan(HR);
   dsfactor = round(Fs/4);
   In2 = downsample(In,dsfactor); %run at 4Hz for speed (e.g. 125 Hz is slow!)
   Timeds=downsample(Time,dsfactor);
   Out2 = RemoveShortSegments(In2,30,dt*dsfactor);
   Out = interp1(Timeds,Out2,Time,'nearest'); %vq = interp1(x,v,xq)
   HR(Out==1)=NaN; 
end
if 1 %find NaN edges and remove beats with excessive deltas
   isnansignal = isnan(HR); 
   InanOn=find(diff(isnansignal)==1);
   for i=1:length(InanOn)
       endHR = HR(InanOn(i));
       x = find(HR(InanOn(i):-1:1)~=endHR,1,'first');
       if isempty(x)
           break
       end
       prevHR = HR(InanOn(i)-x+1);
       if abs(endHR-prevHR)>15 %%%%%%%%%%%%%%%%%%%%%% threshold jump in HR at startNaN to remove
           HR(InanOn(i)-x+2:InanOn(i))=NaN;
       end
   end
   InanOff=find(diff(isnansignal)==-1);
   for i=1:length(InanOff)
       startHR = HR(InanOff(i)+1);
       x = find(HR(InanOff(i)+1:end)~=startHR,1,'first');
       if isempty(x)
           break
       end
       nextHR = HR(InanOff(i)+x);
       if abs(startHR-nextHR)>15 %%%%%%%%%%%%%%%%%%%%%% threshold jump in HR at startNaN to remove
           HR(InanOff(i)+1:InanOff(i)+x-1)=NaN;
       end
   end
end

%% Plot
                     % 'plotfig' switched to zero by DLM
if plotfig % ||Review_data      % 'or Review' commented out by DLM      
% Plot Peak detection
%figure('Name',['RR detect'],'color',[1 1 1]);
figure(5); clf(5); set(gcf,'color',[1 1 1]);
axx2(1)=subplot(3,1,1); 
plot(Time,EKG_,'k',Time(ECG_peak_i),EKG_(ECG_peak_i),'r.');
box('off'); set(gca,'xtick',[]);
ylabel('EKG');

axx2(2)=subplot(3,1,2); 
plot(Time,ECG_filtered,Time(ECG_peak_i_step1),ECG_filtered(ECG_peak_i_step1),'r.',Time(ECG_peak_i),ECG_filtered(ECG_peak_i),'k.',Time,ECG_filtered_threshold__,'r');
box('off'); set(gca,'xtick',[]); 
ylabel('EKG Filt.');
scaling = 4*nanmedian(ECG_filtered_threshold__);
hold('on')
stairs(Time,covRRseries*scaling,'g');
plot(Time(ECG_peak_i),covRRthres*scaling+0*ECG_peak_i,'g:')
%plot(Time,Noise*scaling,'color',[1 0.3 0.2]);
% plot(Time,RMS__,'r--');

plot(Time(ECG_peak_i),4*scaling*FLocalTimeAboveThres,'b');
plot(Time,RMSraw__,'k');

ylimguess = nanmedian(ECG_filtered(ECG_peak_i));
ylim([-1*ylimguess ylimguess*5]);

hold('off')

axx2(3)=subplot(3,1,3); plot(Time,HR);
ylabel('Heart rate (b/min)');
xlabel('Time (s)');
box('off');

linkaxes(axx2,'x');

end

% xlims=get(gca,'xlim')
% try
% set(gca,'xlim',xlims)
% catch me
% end

%% Remove noisy QRS for output
if RemoveNoisyQRSforOutput
    ECG_peak_i(Noise(ECG_peak_i)==1) = [];
end



