% clear all
% close all

if 0
    pathpiece = 'E:\Dropbox (Partners HealthCare)\PAtO\Studies\Scored';
     [spikefile,path]=uigetfile('*.*','Select original spike file',pathpiece);
    load([path spikefile]);

end

settings.windowwidth = 1200; %set section length for determining best Edi
settings.FTestContribution = [1 0.67 0.1 0.2]; %0.2 %LF noise, HF noise, signal-to-tonic, relative signal mean
settings.Fjumptoconstantcal=0.9;
Nsigs=5

% %%
% addpath([pwd '\ConvertedEdi'])
% filename = 'subject4.mat';
% load(filename)

%%
dt=0.008;
N=length(EMGdi1);
Time = [0:dt:(dt*(N-1))]';

figure(1); clf(1);
ds=5;
gap=0.05;
for i=1:Nsigs
plot(downsample(Time,ds),downsample(eval(['EMGdi' num2str(i)]),ds)-gap*(i-1));
hold('on')
end


%% Tool to examine factors in selected window from Fig 1

if 0
figure(1)
xlims=get(gca,'xlim')
II = round(xlims/dt+1)
if II(2)>length(EMGdi1)
    II(2)=length(EMGdi1)
end
I = (II(1):II(2))';
dataset = [EMGdi1(I) EMGdi2(I) EMGdi3(I) EMGdi4(I) EMGdi5(I)];
prctile(dataset,20)
time = Time(I);
    df = 1/(time(end)-time(1));
    F = 0:df:(length(time)-1)*df; Fmax = 20;
    Fmaxi = round(Fmax/df+1);
    ranges = [0.01 0.1 1 10];
    rangesi = round(ranges/df)-1;  
        rangesi(rangesi==0)=1;
    clear P1 P2 P3 StoN
for i=1:Nsigs
    sig = dataset(:,i);
    sigfft = fft(sig-mean(sig))*2/length(sig);
    sigpsd = sigfft.*conj(sigfft)*2/length(sig);
    F((Fmaxi+1):end)=[]; 
    sigpsd((Fmaxi+1):end)=[]; 
    sigpsd(1)=0;
    P1(i) = sum(sigpsd(rangesi(2):rangesi(3))); %note not divided by bandwidth
    P2(i) = sum(sigpsd(rangesi(3):rangesi(4))); %note not divided by bandwidth
    P3(i) = sum(sigpsd(rangesi(1):rangesi(2))); %note not divided by bandwidth  
end
tests=[P1./P3;P1./P2]
tests=tests./mean(tests,2);
0.5*tests(1,:)+0.5*tests(2,:)

figure(2)
scatter(EMGdi3(I),EMGdi4(I),10,[0 0 0],'filled','markerfacealpha',0.1)
prctile(EMGdi3(I)./EMGdi4(I),90)

end

%% Find best Edi
wwi = settings.windowwidth;
wwi = round(wwi/dt);

%smooth transition window
    df = 1/(Time(wwi)-Time(1));
    F = 0:df:(wwi-1)*df; Fmax = 20;
    Fmaxi = round(Fmax/df+1);
    F((Fmaxi+1):end)=[]; 
    ranges = [0.001 0.1 1 10]; %was0.005
    rangesi = round(ranges/df)-1;  
    rangesi(rangesi==0)=1;

    TimeBuffer = buffer(Time,wwi,0,'nodelay');
    
    clear P1 P2 P3 Tonic
for i=1:Nsigs
    sig = buffer(eval(['EMGdi' num2str(i)]),wwi,0,'nodelay');
    %sig(isnan(sig))=0;
    sigfft = fft(sig-nanmean(sig))*2/wwi;
    sigpsd = sigfft.*conj(sigfft)*2/wwi;
    sigpsd((Fmaxi+1):end,:)=[]; 
    sigpsd(1,:)=0;
    P1(:,i) = abs(sum(sigpsd(rangesi(2):rangesi(3),:))).^0.5; %note not divided by bandwidth
    P2(:,i) = abs(sum(sigpsd(rangesi(3):rangesi(4),:))).^0.5; %note not divided by bandwidth
    P3(:,i) = abs(sum(sigpsd(rangesi(1):rangesi(2),:))).^0.5; %note not divided by bandwidth  
    %eval(['EMGdi' num2str(i) 'Buffer=sig;']);
    Tonic(:,i)=prctile(sig,20);
end
P1onP3=[P1./P3]; %sig (P1) over LF noise
P1onP2=[P1./P2]; %sig (P1) over HF noise
P1onTonic=[P1./Tonic]; %sig (P1) over tonic noise
P1alone_ = P1./mean(P1,2); %sig (P1) compared to mean of others in absolute units

if 0
P1onP3_=P1onP3./mean(P1onP3,2);
P1onP2_=P1onP2./mean(P1onP2,2);
P1onTonic_=P1onTonic./mean(P1onTonic,2);
else
P1onP3_=P1onP3./2.07;
P1onP2_=P1onP2./9.11;
P1onTonic_=P1onTonic./0.011;    
end

if 1
Test = settings.FTestContribution(1)*P1onP3_ + ...
    settings.FTestContribution(2)*P1onP2_ + ...
    settings.FTestContribution(3)*P1onTonic_ + ...
    settings.FTestContribution(4)*P1alone_;
else
    Test = (settings.FTestContribution(1)*P1onP3_ + ...
    settings.FTestContribution(2)*P1onP2_ + ...
    settings.FTestContribution(3)*P1onTonic_).*P1alone_;
end

[~,BestI] = max(Test');
BestI=BestI(:);
TestNaN=isnan(sum(Test,2));
BestI(TestNaN==1)=NaN;

Time_ = TimeBuffer(1,:)';
figure(1);clf
subplot(6,1,1)
h2 = stairs(Time_,BestI)
%h2 = stairs(Time_,-(BestI-1)*gap)
for j=1:5
subplot(6,1,j+1)
hold on
h2 = stairs(Time_,P1onP3_(:,j))
h2 = stairs(Time_,P1onP2_(:,j))
h2 = stairs(Time_,P1onTonic_(:,j))
%h2 = stairs(Time_,P1alone_(:,j)*gap)
ylim([0 2])
end
%delete(h)

L = 2*wwi;
win = hann(L);
Left = round(L/4);
win = [zeros(Left,1);win(1:end-Left)];
win((end-Left+1):end)=1; win(end)=1;

%% Try to update BestI with neighboring values through the unknowns
if 1
BestI = interp1(Time_(~isnan(BestI)),BestI(~isnan(BestI)),Time_,'nearest'); %remove extrap because trouble in Fadjust if there are NaN in the actual signal
end

%% Make a single signal EMGdiX
wincurrent = BestI(find(~isnan(BestI),1,'first'));
winnext = wincurrent; %default, in case of NaN in win 2
currentF = 1;
EMGdiX = NaN*EMGdi1;

BestSingleEdi = mode(BestI)
currentFbest = [];
h=plot(Time,EMGdiX-gap*Nsigs,'k');
EMGdiX(1:wwi,1)=eval(['EMGdi' num2str(wincurrent) '(1:wwi)']);
adjustF=1;
clear currentFbackup
for i=2:size(TimeBuffer,2)
adjustF(i)=1;
irange = ((i-1)*wwi)+1:(i*wwi);
irangeincllast = ((i-2)*wwi)+1:(i*wwi);
irange(irange>length(EMGdiX))=[];
irangeincllast(irangeincllast>length(EMGdiX))=[];
if ~isnan(BestI(i-1))
    wincurrent = BestI(i-1);
end %else do not update
if ~isnan(BestI(i)) 
    winnext = BestI(i);
end %else do not update


EMGdiX(irange) = currentF*eval(['EMGdi' num2str(wincurrent) '(irange)']);   

% if i==17   
%     pause()
% end
if i==size(TimeBuffer,2) || isnan(BestI(i-1))%|| winnext==wincurrent
    continue
end

datalast = eval(['EMGdi' num2str(wincurrent) '(irangeincllast)']);   
datanext = eval(['EMGdi' num2str(winnext) '(irangeincllast)']);  

prct=75;
I=datalast>prctile(datalast,prct)&datanext>prctile(datanext,prct);
adjustF(i) = nanmedian(datalast(I)./datanext(I));

currentF=currentF*adjustF(i);

%will prevent any progressive drift in signal amplitude, so long as best signal is chosen at multiple time points, and relies on best signal holding
%absolute calibration:
if isempty(currentFbest)&&winnext==BestSingleEdi
    currentFbest=currentF;
end

if ~isempty(currentFbest)&&winnext==BestSingleEdi
    currentF=settings.Fjumptoconstantcal*currentFbest+(1-settings.Fjumptoconstantcal)*currentF;
end

newdatanext = currentF*datanext;

EMGdiX(irangeincllast) = EMGdiX(irangeincllast).*(1-win) + newdatanext.*(win); 

end 
%     delete(h);
% h=plot(Time,EMGdiX-gap*5,'k');
prct=75
EMGdiBest = eval(['EMGdi' num2str(BestSingleEdi)]);
I=EMGdiX>prctile(EMGdiX,prct)&EMGdiBest>prctile(EMGdiBest,prct);
adjustFfinal = median(EMGdiBest(I)./EMGdiX(I));

EMGdiX = EMGdiX*adjustFfinal;
h=plot(Time,EMGdiX-gap*(BestSingleEdi-1),'k');

%%
BestEdiInfo.BestSingleEdi = BestSingleEdi;
BestEdiInfo.BestI = BestI;
BestEdiInfo.Time = Time_;
%EMGdiCh=interp1(Time_,BestI,Time,'previous'); %channel selected for EMGdiX w time
EMGdiCh=interp1(Time_,BestI,Time,'nearest'); %channel selected for EMGdiX w time
%%
%save(filename,'EMGdiX','EMGdiBest','settings','BestEdiInfo','-append');
  %save([path 'FiltEMGdi' '.mat'],'EMGdiX','EMGdiCh','EMGdiBest','settings','BestEdiInfo','-append','-v7.3');
 savevarlist={'EMGdiX','EMGdiCh','EMGdiBest','settings','BestEdiInfo'}
   savevarlist2={'EMGdi1','EMGdi2','EMGdi3','EMGdi4','EMGdi5'}
          save([path 'FiltEMGdi' '.mat'], savevarlist2{:} , savevarlist{:},'-v7.3');
%       
%    
%         if exist([path 'FiltEMGdi' '.mat'],'file')==0
%         %    save(['subjectA' num2str(subject) '.mat'],savevarlist2{:},'-v7.3');          
%             save([path 'FiltEMGdi' '.mat'], savevarlist2{:} , savevarlist{:},'-v7.3');
%         else
%             save([path 'FiltEMGdi' '.mat'],savevarlist{:},'settings','BestEdiInfo','-append','-v7.3');
%         end
        %%
        figure(10);clf
        ax(1)=subplot(2,1,1)
        plot(EMGdiCh);hold on
        ax(2)=subplot(2,1,2)
        plot(EMGdiX);hold on
        
        load([path spikefile]);
        ax(1)=subplot(2,1,1)
        plot(EMGdiCh);hold on
        ax(2)=subplot(2,1,2)
        plot(EMGdiX);hold on
        linkaxes(ax,'x')
        
  %% Debug plots
  if 0  
      
N=length(EMGdiX);
dt = 0.008;
Time = [0:dt:(N-1)*dt]';
EdiT = table(Time,EMGdi1,EMGdi2,EMGdi3,EMGdi4,EMGdi5,EMGdiCh,EMGdiX);

signallist={'EdiT.EMGdi1','EdiT.EMGdi2','EdiT.EMGdi3','EdiT.EMGdi4','EdiT.EMGdi5','EdiT.EMGdiCh','EdiT.EMGdiX'};

global xvalues yvalues range ax2
range=10000;
PlotAndSelectData(signallist,Time,1);

  end