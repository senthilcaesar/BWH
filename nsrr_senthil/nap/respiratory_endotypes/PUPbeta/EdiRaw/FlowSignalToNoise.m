function noisewav = FlowSignalToNoise(timewav,respwav,showfigures,AllowNoisierFlowFactor)
% modified by DLM.
%   added "round" to lefti and righti
%   added check on nan in flow signal
%   added check on P1scale, if NaN try setting lower, only switched on if AllowNosierFlowFactor


CheckP1Scale = 0; % default
if ~exist('AllowNoisierFlowFactor')
    AllowNoisierFlowFactor=1; %default 0.25
else
    CheckP1Scale = 1;
end
FlowSignalNoiseFrequencyShift=1;

Fs = 1./(timewav(2)-timewav(1));

%signal to noise per window
secslide2=10;
windur2 = 120; %note no savings in speed by keeping nfft = 2^X (actually lost speed).
% nfft is usually set as next larger power of 2, as:
% nfft = 2^nextpow2(N); % next larger power of 2
nwin=floor((timewav(end)-timewav(1)-windur2)/secslide2 + 1);
StoN_x=zeros(1,nwin);
P1_x=zeros(1,nwin);
Time_x=zeros(1,nwin);
for i=1:nwin%winNum=winnumrange
    lefti=round(1+(i-1)*secslide2*Fs);
    righti=round(lefti+windur2*Fs-1);
    Flow=respwav(lefti:righti);
    Time=timewav(lefti:righti);
    Time_x(i) = Time(1);
    
    if i==1
        df = 1/(Time(end)-Time(1));
        F = 0:df:(length(Time)-1)*df; Fmax = 20;
        Fmaxi = round(Fmax/df+1);
        ranges = [0.1 1 1 10]*FlowSignalNoiseFrequencyShift;
        rangesi = round(ranges/df)-1;
    end
    
    % DLM added the ~any(isnan(flow)) line so that we only do fft
    % and psd on Flow if it doesn't contain NaN. I added this because
    % an fft/psd on NaN data returns complex results, which breaks
    % everything below
    if ~any(isnan(Flow))
        Flowfft = fft(Flow-mean(Flow))*2/length(Flow);
        Flowpsd = Flowfft.*conj(Flowfft)*2/length(Flow);
        F((Fmaxi+1):end)=[]; Flowpsd((Fmaxi+1):end)=[]; Flowpsd(1)=0;
        % should we do a conj at this point, to get absolute values?
        P1 = sum(Flowpsd(rangesi(1):rangesi(2))); %note not divided by bandwidth
        P2 = sum(Flowpsd(rangesi(3):rangesi(4))); %note not divided by bandwidth
        StoN = P1/P2;
        StoN_x(i)=StoN;
        P1_x(i) = P1;
    end
    
end

P1scale = 10^mean(log10(P1_x(StoN_x>10))); % original line, uses fixed value of 10

% DLM also added the following if CheckP1Scale block. 
% the P1Scale line above could return NaN in 'reasonable' data, which then breaks everything below
% so this check is enabled if we explicitly call the function with AllowNoisierFlowFactor set to something
if CheckP1Scale
    if isnan(P1scale)
        upper95C = prctile(StoN_x, 95);
        disp(['WARNING: SNR did not reach >=10 for P1scale, using lower threshold setting of ', num2str(upper95C)]);
        P1scale = 10^mean(log10(P1_x(StoN_x>upper95C))); % modified line, uses 95th centile of data
    end
end

criteria = P1_x/P1scale.*StoN_x;
noisewav=0*timewav;
for i=1:nwin
    if criteria(i)<(1/4/AllowNoisierFlowFactor)
        lefti=round(1+(i-1)*secslide2*Fs);
        righti=round(lefti+windur2*Fs-1);
        noisewav(lefti:righti)=1;
    end
end
for i=1:nwin
    if criteria(i)<(0.5/4/AllowNoisierFlowFactor)
        lefti=round(1+(i-1)*secslide2*Fs);
        righti=round(lefti+windur2*Fs-1);
        noisewav(lefti:righti)=2;
    end
end
for i=1:nwin
    if criteria(i)<(0.25/4/AllowNoisierFlowFactor)||isnan(criteria(i))
        lefti=round(1+(i-1)*secslide2*Fs);
        righti=round(lefti+windur2*Fs-1);
        noisewav(lefti:righti)=3;
    end
end
StoNData.Fnoiseovernight = [sum(noisewav>=1) sum(noisewav>=2) sum(noisewav>=3)]/length(noisewav);
StoNData.StoNovernight = [sum(StoN_x<20) sum(StoN_x<10) sum(StoN_x<5)]/length(noisewav);
if showfigures %StoN figure
    figure(101);
    set(gcf,'color',[1 1 1]);
    ax30(1)=subplot(4,1,3); stairs(Time_x,P1_x/P1scale);
    set(gca,'yscale','log');
    hold('on');
    stairs(Time_x,StoN_x,'k'); box('off');
    stairs(Time_x,P1_x/P1scale.*StoN_x,'r'); box('off');
    stairs(timewav,10.^noisewav,'g'); box('off'); hold('off');
    legend('Power','StoN','PxStoN','lowQ');
    ax30(2)=subplot(4,1,4); plot(timewav,respwav); ylabel('Flow'); box('off');
    %ax30(3)=subplot(3,1,3); plot([1:nwin],NoiseData); ylabel('Noise'); box('off');
    linkaxes(ax30(1:2),'x');
end
