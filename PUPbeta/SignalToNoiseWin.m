%BreathDataTable = SignalToNoiseWin(BreathDataTable,Flow-leak, Time)

function BreathDataTable = SignalToNoiseWin(BreathDataTable,SigIn,Time,DSflag)

%% For use as an addon:
%get the 7 min win of Flow from DataEventHypnog_Mat (converted), called SigIn
%pass in BreathDataTable{i}
%function will append cols to BreathDataTable{i}

%% Define and load subject data
global settings

%% What signals
dt = 1/settings.Fs;

if 0
    Flow = SigIn*(1-2*settings.Pnasaldownisinsp); %does this matter
end

%edit to make this work for any sampling rate
if exist('DSflag') && DSflag==1 %when operated in LGfromFlowBeta this is already done; %see DSflag=1 for use on raw (converted) data
    [p,q] = rat(settings.downsampledFs/(1/dt));
    Flow_ds = resample(SigIn,p,q);
    Flow_ds2 = resample(Flow_ds,q,p);
    Flow=Flow_ds2;
end

SNR_breath = 0*BreathDataTable.BB_i_start;
SNR_window = 0*BreathDataTable.BB_i_start;


%% Window SNR
% WinList = unique(BreathDataTableLong.Time0);
% % Am_breaths = 0*WinList;
% % for p = 1:length(WinList)
% %     Am_breaths(p) = length(find(BreathDataTableLong.Time0 == WinList(p)));
% % end
% %  %SignalBreathdB = zeros(max(Am_breaths),length(WinList));
SignalBreathdB = 0*BreathDataTable.BB_i_start;

li = round((BreathDataTable.Time_start(1)-Time(1))/dt + 1);
ri = round((BreathDataTable.Time_end(end)-Time(1))/dt + 1);
I = li:ri;


% exponent = settings.scalingexponent;
% leak = BreathDataTableLong.leak(jj(1));
% leak2 = BreathDataTableLong.leak2(jj);
% IEratio = BreathDataTableLong.IEratio(jj(1));

% figure(99); clf(99);
% ax99(1)=subplot(3,1,1);
% plot(Time(I),Flow(I)-leak);

% ax99(2)=subplot(3,1,2);
% plot(Time(I),Vflow_out2);
% linkaxes(ax99,'x')

if 0 %remove data near zero crossings (obsolete)
    Izero = [BreathDataTable.BB_i_start(:); BreathDataTable.BB_i_mid(:) ; BreathDataTable.BB_i_end((end))];
    Izero = Izero - BreathDataTable.BB_i_start(1)+1;
    remdT = 0.1;
    remdTi = round(remdT/dt/2);
    remI = repmat([-remdTi:remdTi],length(Izero),1) + Izero;
    remI = remI(:);

    remI(remI<1)=[];
    remI(remI>length(I))=[];
end
%
%
%log10(48)/log10(2)
Nfft=round(0.33/dt); %tailored to sampling rate
Noverlap=round(Nfft*0.75);

filter_HFcutoff_butter0 = 4;
filter_order0 = 1;
[B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'high');
SigInF = nanfilter(B_butter0,A_butter0,SigIn-nanmean(SigIn),1);
%SigInF = filtfilt(B_butter0,A_butter0,SigIn-nanmean(SigIn)); %filtfilt, otherwise flow signal is right-shifted

filter_HFcutoff_butter0 = 4;
filter_order0 = 1;
[B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
SigInFlow = nanfilter(B_butter0,A_butter0,SigIn,1);
%SigInFlow = filtfilt(B_butter0,A_butter0,SigIn); %filtfilt, otherwise flow signal is right-shifted

if 0  %remove data near zero crossings (obsolete)
    SigInF2 = SigInF;
    SigInF2(remI)=NaN; 
end

% figure(99); clf(99);
% ax99(1)=subplot(3,1,1);
% plot(Time(I),[SigIn SigInF SigInF2])
% 
% ax99(2)=subplot(3,1,2);
% plot(Time(I),abs(gradient(SigInFlow)))
% linkaxes(ax99,'x')
% Pxx
t=(1:length(SigInF))';

FlowBuffer = buffer(SigInF,Nfft,Noverlap,'nodelay');
%FlowBuffer=detrend(FlowBuffer);
%FlowBuffer2 = buffer(Flow(I)-leak,Nfft,Noverlap,'nodelay');
IndexBuffer = buffer(t,Nfft,Noverlap,'nodelay');
if 0
    hannwin = hann(Nfft); hannwin = hannwin./rms(hannwin);
    FlowBuffer = FlowBuffer.*hannwin;
end
Frange_low = [5 12];
T01 = Nfft*dt; 
df1 = 1/T01;
Frangei1 = round(Frange_low/df1 + 1); %frequency bins of interest based on frequency range
trueF = (Frangei1-1)*df1 + df1*[-0.5 0.5]';

X1 = fft(FlowBuffer)/Nfft*2;
Pxx1 = abs(X1.^2); %same as Pxx1 = conj(X1).*X1; 
Pow1 = sum(Pxx1(Frangei1(1):Frangei1(2),:))*df1; %Pow1 is power in 5-12 Hz, moving signal

% figure; plot(1/df1:1/df1:length(Pxx1')/df1, Pxx1')

% make this faster in future by making it directly from Time signal
TimeBuffer = buffer(Time(I),Nfft,Noverlap,'nodelay');
%TimePxx = TimeBuffer(1,:) + T01/2;
TimePxx = nanmedian(TimeBuffer);
clear TimeBuffer

NoiseThreshold = log10(prctile(Pow1,25));

% figure(99)
% ax99(2)=subplot(3,1,3);
% plot(TimePxx,[log10(Pow1);NoiseThreshold+0*TimePxx]);
% linkaxes(ax99,'x')
% hold on

% DLM comments:
% Turns out, it's not to do with matlab version at all, but rather that detrend is an overloaded function. 
% on one hand, it subtracts the mean or best-fit line from a timeseries object,
% on another hand, it removes a polynomial trend.
%
% the mean subtract function supports up to three parameters, being:
% (1) the data,
% (2) the method [linear or constant], and 
% (3) row or column indices.
% 
% the polynomial trend function supports lots of parameters:
% (1) the data,
% (2) nth-degree polynomial trend, (not optional)
% (3) breakpoints, (optional)
% (4) nanflag (e.g. 'omitnan')
% (...) name,value pairs.
% 
% so, we have to add a 0 to the polynomial detrend to do mean subtraction
try
    SignalWin = 2*log10(nanmean(detrend(SigInFlow,0,'omitnan').^2).^0.5); %nanrms(detrend()) DLM added 0 for mean detrend
catch 
    SignalWin = 2*log10(nanmean((SigInFlow-nanmean(SigInFlow)).^2).^0.5); % DLM removed omitnan, was attempting to apply to 'detrend', but possibly was menat for 'mean', however not reqd with nanmean.
end

StoNWindB = 10*(SignalWin-NoiseThreshold);

for k=1:length(BreathDataTable.BB_i_start)
    Ii = [BreathDataTable.BB_i_start((k)):BreathDataTable.BB_i_end((k))]-BreathDataTable.BB_i_start((1))+1;
    SignalBreathdB((k)) = 2*log10(rms(detrend(SigInFlow(Ii))));
end

StoNWinBreathdB = 10*(SignalBreathdB-NoiseThreshold);
%temp = 10*(SignalBreathdB-log10(nanmean(10.^SignalBreathdB)));
%amp = (10.^(temp/10)).^0.5;
% figure(10);
% plot(Time(I),SigIn)
% hold on
% stairs(BreathDataTableLong.Time_start(jj),StoNWinBreathdB/max(StoNWinBreathdB)*max(SigIn))


%SNR_window = mean(StoNWinBreathdB);

BreathDataTable.SNRwindow = StoNWindB*ones(length(BreathDataTable.Time0),1);
BreathDataTable.SNRbreath = StoNWinBreathdB(:);

%FbreathOverWin =
%(10.^([BreathDataTable.SNRbreath-BreathDataTable.SNRwindow]/10)).^0.5;

end
