function [ECG_peak_i,HR] = EKGpeakdetection(EKG,Time,dt)

%% Find RR
fs = round(1/dt);
N = length(Time);
%Filter EKG to aid peak identification. 
filter_LFcutoff_butter = 5; %highpass, DC remove
filter_HFcutoff_butter = 15; %lowpass, smooth
filter_order = 2; 
[B_butter,A_butter] = butter(filter_order,[filter_LFcutoff_butter filter_HFcutoff_butter]/(fs/2));
ECG_filtered = filtfilt(B_butter,A_butter,EKG);

% ECG_flip_on=1;
% FlipECG = ((mean(ECG_filtered)-prctile(ECG_filtered,0.5))>(prctile(ECG_filtered,99.5)-mean(ECG_filtered)));
% 
    Fs = 1/(Time(2)-Time(1));
    
    % the windowduration and secslide values were hardcoded here
    % they did not use the values from the options worksheet
    secslide=120;
    WindowDuration=180;
    
    skipCPAPon = 1;
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    %FlipEKGvotes = NaN*ones(numwind,1);
    FlipEKGmag = NaN*ones(numwind,1);
    for winNum=1:numwind
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        I=round((li:ri));
        %FlipEKGvotes(winNum) = ((mean(ECG_filtered(I))-prctile(ECG_filtered(I),0.5))>(prctile(ECG_filtered(I),99.5)-mean(ECG_filtered(I))));
        FlipEKGmag(winNum) = ((mean(ECG_filtered(I))-prctile(ECG_filtered(I),0.5))-(prctile(ECG_filtered(I),99.5)-mean(ECG_filtered(I))));
    end
    FlipECG = nansum(FlipEKGmag)>0;
    
    
        if FlipECG
            ECG_filtered=-ECG_filtered;
            EKG=-EKG;
            disp('flipped EKG');
        end
    
%% window-by-window thresholding
ECG_filtered_threshold_fraction=0.7;

if Time(end)-Time(1)>1200 %if shorter window then don't use time-varying threshold
    Fs = 1/(Time(2)-Time(1));
    
    % the windowduration and secslide values were hardcoded here
    % they did not use the values from the options worksheet
    secslide=120;
    WindowDuration=180;
    
    skipCPAPon = 1;
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    ECG_filtered_threshold_ = NaN*ones(numwind,1);
    ECG_filtered_threshold_t = NaN*ones(numwind,1);
    for winNum=1:numwind
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        I=round((li:ri));
        ECG_filtered_threshold_(winNum) = ECG_filtered_threshold_fraction*[prctile(ECG_filtered(I),99)-prctile(ECG_filtered(I),1)]+prctile(ECG_filtered(I),1);
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

%%

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


%% Remove array elements (RAE) [fast]
    %Removes EKG peaks that are too close to each other
    %The smaller of the two peaks are removed. 
  
%RAE Setup
ECG_peak_t_threshold_percentofmedian = 33; %specify

Nrejected = [];
for X=1:3 %runs several times. 2 cases usually enough.
RRa = [diff(ECG_peak_i)*dt NaN];
RRb = [NaN RRa(1:end-1)];
RRc = [NaN NaN RRa(1:end-2)];
data = [RRa;RRb;RRc];

covRR = std(data)./mean(data); %local cov

%find normal RR
upper = 4;
lower = 0.25;
covmax = 1/6;
ECG_peak_t_delta_median = nanmedian(RRa(RRa<upper&RRa>lower&RRa>lower&covRR<covmax));
ECG_peak_t_threshold_min = ECG_peak_t_delta_median*ECG_peak_t_threshold_percentofmedian/100;
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
    secslide=5;
    Fs = 1/(Time(2)-Time(1));
    WindowDuration=15;
    
if Time(end)-Time(1)>WindowDuration %if shorter window then don't use time-varying threshold
    numwind=floor((length(ECG_filtered)-WindowDuration*fs)/(secslide*fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    %SpikinessW = zeros(1,numwind);
    %timetemp = zeros(1,numwind);
    Minabsdiffseries  = 0*Time;
    Spikiness = 0*Time;
    for winNum=1:numwind
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        I=round((li:ri));
        %timetemp(winNum)=mean(Time(I));
        Spikiness(li:ri) = nanmean(Peakheights(I))/prctile((ECG_filtered(I)),80)<6;
        Minabsdiffseries(li:ri) = nanmean(minabsdiffseries(I))>0.5;
        if nanmean(covRRseries(I))>covRRthres||(max(RRseries(I))>maxRRthres)%||(sum(isnan(RRseries(I)))>0)
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

%% Make short lengths of good data => noise

minlength=60;

I1=find(diff(1-Noise)==1); 
I2=find(diff(1-Noise)==-1); 
[I1,I2] = TidyStartEndEventList(I1,I2,length(Noise));
lengths = (I2-I1)*dt;

for i=1:length(I1)
    if lengths(i)<minlength
        Noise(I1(i):I2(i))=1;
    end
end

%% Move index to peak of unfiltered EKG
delta=4;
for i=1:length(ECG_peak_i)
    rangei = ECG_peak_i(i)-delta:ECG_peak_i(i)+delta;
    
    rangei(rangei<1)=[];
    rangei(rangei>length(EKG))=[];
    
    [~,temp] = max(EKG(rangei));
    ECG_peak_i(i)=ECG_peak_i(i)+temp-delta-1;
end

%% Calculate RR (quadratic fit to peak to improve RR resolution from 1/fs)

ECG_peak_i = unique(ECG_peak_i);

ECG_peak_quadfit_t = NaN*ECG_peak_i;
for i=1:length(ECG_peak_i)
    I = ECG_peak_i(i)-1:ECG_peak_i(i)+1;
    ECG_peak_quadfit_t(i)=PeakFitQuadratic(Time(I),EKG(I));
end

RR = [diff(ECG_peak_quadfit_t) NaN];
HR = zoh(Time(ECG_peak_i),60./RR',Time); %very fast version
HR(Noise==1)=NaN;

%% Plot
plotfig=0;                      % 'plotfig' switched to zero by DLM
if plotfig % ||Review_data      % 'or Review' commented out by DLM
% Plot Peak detection
%figure('Name',['RR detect'],'color',[1 1 1]);
figure(5); clf(5); set(gcf,'color',[1 1 1]);
axx2(1)=subplot(3,1,1); 
plot(Time,EKG,'k',Time(ECG_peak_i),EKG(ECG_peak_i),'r.');
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
plot(Time,Noise*scaling,'color',[1 0.3 0.2]);
hold('off')

axx2(3)=subplot(3,1,3); plot(Time,HR);
ylabel('Heart rate (b/min)');
xlabel('Time (s)');
box('off');

linkaxes(axx2,'x');

end


%% Remove noisy QRS for output
ECG_peak_i(Noise(ECG_peak_i)==1) = [];



