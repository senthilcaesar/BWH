function [Data,DataCI,DataN,varlist,AHItotal,Fstates] = SummaryAnalysisOnePes(subj,settings)

%% Hardcoded directory info
directory=settings.directory;%'H:\MESA\Analyzed3';
savename =settings.savename;%'MESA2Sept2017';

scoredarousalsinwake=1;

%% Default for figures
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');

%% Load PUPdata

%subj=1;
getCIs=settings.getCIs;
varlist = {'LG1','LGn','delay','VRA1','ArThres','Vpassive','Vactive','Vcomp'};
Nperwindowvars=5;

loadpath=[directory '\' savename '_' num2str(subj)]
load(loadpath,'AHIdata','SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual')

%M=1;

%%
AHItotal=AHIdata{1}(58);

%%
state=settings.selectstate; %1=nrem1, 4 nrem, 5 rem, 8 ignore state

%% Analyse LG, ArThres

clear LG1 LG2 LGn temp temp2 FVAL VRA1 VRA2 ArThres MeanEx

% events, longestwake, position, FremMAX
minNevents = 1; %1
maxwakethres = 30;

%defaults
minFnrem1=-Inf;
minFnrem2=-Inf;
minFnrem3=-Inf;

if state==4
    maxFREM=0;
    minFrem=-Inf;
    AHIstate=AHIdata{1}(88);
    DurationInState=AHIdata{1}(88-7);
elseif state==1
    maxFREM=0;
    minFrem=-Inf;
    minFnrem1=0.5;
    AHIstate=AHIdata{1}(88); %not correct
    DurationInState=AHIdata{1}(88-7);
elseif state==2
    maxFREM=0;
    minFrem=-Inf;
    minFnrem2=0.5;
    AHIstate=AHIdata{1}(88); %not correct
    DurationInState=AHIdata{1}(88-7);
elseif state==3
    maxFREM=0;
    minFrem=-Inf;
    minFnrem3=0.5;
    AHIstate=AHIdata{1}(88); %not correct
    DurationInState=AHIdata{1}(88-7);
elseif state==5
    maxFREM=Inf;
    minFrem=0.5;
    AHIstate=AHIdata{1}(96);
    DurationInState=AHIdata{1}(96-7);
elseif state==8
    maxFREM=Inf;
    minFrem=-Inf;
    AHIstate=AHIdata{1}(80);
    DurationInState=AHIdata{1}(80-7);
end

%% Find first and last windows containing sleep
TimeOfNightCriteria = ones(length(SleepData{1}(:,1)),1);
if settings.selecttimeofnight
    XTiles=settings.selecttimeofnight_XTiles;
    NthXTile=settings.selecttimeofnight_NthXTile;
    
    FwakePerWindow=SleepData{1}(:,1);
    SleepWin1=find(FwakePerWindow<1,1,'first');
    SleepWinN=find(FwakePerWindow<1,1,'last');
    
    rangeWin=round((SleepWinN-SleepWin1+1)/XTiles);
    lowerWinNs=round(SleepWin1+((1:XTiles)-1)*rangeWin);
    upperWinNs=round(lowerWinNs+rangeWin-1);
    
    TimeOfNightCriteria=0*TimeOfNightCriteria;
    TimeOfNightCriteria(lowerWinNs(NthXTile):upperWinNs(NthXTile))=1;
end

%%
usemediannotmeanLG1=1;
i=1;
clear Tn LG1
tempLG1=LGplusinfo{i}(:,7);
FisNaN(i)=sum(isnan(tempLG1))/length(tempLG1);

temp=LGplusinfo{i};

tempLG1=LGplusinfo{i}(:,7); tempLGn=LGplusinfo{i}(:,5);  tempLG2=LGplusinfo{i}(:,8); tempTn=LGplusinfo{i}(:,6);
tempLG0=LGplusinfo{i}(:,2); temptau=LGplusinfo{i}(:,3);  tempdelay=LGplusinfo{i}(:,9);
tempVRA1=LGplusinfo{i}(:,10); tempVRA2=LGplusinfo{i}(:,11); tempArThres=LGplusinfo{i}(:,12);

temppositiondata=LG_QualityInfo{i}(:,5);
templongestwakedata=SleepData{i}(:,7);
tempFREM=SleepData{i}(:,6);

Pos = LG_QualityInfo{1}(:,5);
if isempty(settings.selectposition)
    Poscriteria = ones(length(Pos),1);
else
    Poscriteria = sum(Pos == settings.selectposition,2)>0;
end

Fnrem1=SleepData{1}(:,3);
Fnrem2=SleepData{1}(:,4);
Fnrem3=SleepData{1}(:,5);
nremXcriteria=Fnrem1>=minFnrem1&Fnrem2>=minFnrem2&Fnrem2>=minFnrem2;

templongestwakedata(length(tempLG1)+1:end)=[];

tempFREM(length(tempLG1)+1:end)=[];
%SleepData(winNum+1,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];

alpha=(tempLG2./tempLG1).^2;
beta=(1-alpha)./(4*alpha-1);
tempLG3min=tempLG1.*((1+beta)./(1+beta*(1/3)^2)).^0.5; %error, needs fixing...
tempLG6min=tempLG1.*((1+beta)./(1+beta*(1/6)^2)).^0.5; %error, needs fixing...
tempLG90s=tempLG1.*((1+beta)./(1+beta*(2/3)^2)).^0.5; %error, needs fixing...

%GetRsquared:
%temp3=fitQual{i};
clear OneMinusRsq
for ii=1:size(fitQual{i},2)
    temp3pt1=fitQual{i}{ii};
    if length(temp3pt1)>1
        OneMinusRsq(ii,:)=1-temp3pt1(2); %Element #2 of FitQual
    else
        OneMinusRsq(ii,:)=NaN;
    end
end
% Rsq=temp3(:,1);
FVAL=LGplusinfo{i}(:,13);

MeanE=LG_QualityInfo{i}(:,3); %is actual scored events. LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position EXITFLAG];
SStot=FVAL./OneMinusRsq;

N_events=LG_QualityInfo{i}(:,2); tempN_arousals=LG_QualityInfo{i}(:,1);

criteria=Poscriteria==1&nremXcriteria==1&TimeOfNightCriteria==1&templongestwakedata<=maxwakethres&N_events>=minNevents&tempFREM<=maxFREM&tempFREM>=minFrem; %(temppositiondata==0|temppositiondata==x(4))&


% New:
tempVRA1(tempN_arousals<1)=NaN;
tempVRA2(tempN_arousals<1)=NaN;


LG1(i,1)=prctile(tempLG1(~isnan(tempLG1)&criteria),50);
LG1values=tempLG1(~isnan(tempLG1)&criteria);

LG1_SD(i,1)=std(tempLG1(~isnan(tempLG1)&criteria));
LG1_N(i,1)=sum(~isnan(tempLG1)&criteria); % no_of_LGmeasures(i,1)=sum(~isnan(tempLG1)&criteria);
LG1_SE(i,1)=LG1_SD(i,1)/LG1_N(i,1)^0.5;
LG1_N_nocriteria(i,1)=sum(~isnan(tempLG1)); % no_of_LGmeasures(i,1)=sum(~isnan(tempLG1)&criteria);

meanSE=mean(LG1_SE);
F_incl=LG1_N./LG1_N_nocriteria;
meanFincl=mean(F_incl);
minLGN=min(LG1_N);

%minmeanmedianNmeasures=[min(no_of_LGmeasures) mean(no_of_LGmeasures) median(no_of_LGmeasures)]
% LGn_=tempLGn(criteria);
% VRA1_=tempVRA1(criteria);
% VRA2_=tempVRA2(criteria);
MeanEx(i,1)=prctile(MeanE(~isnan(MeanE)&criteria),50);
% LG1(i,1)=prctile(tempLG1(criteria),50)
LG2(i,1)=prctile(tempLG2(~isnan(tempLG2)&criteria),50);
LGn(i,1)=prctile(tempLGn(~isnan(tempLGn)&criteria),50);
LGnvalues=tempLGn(~isnan(tempLGn)&criteria);
Tn(i,1)=prctile(tempTn(~isnan(tempTn)&criteria),50);
VRA1(i,1)=100*prctile(tempVRA1(~isnan(tempVRA1)&criteria),50);
VRA1values=100*tempVRA1(~isnan(tempVRA1)&criteria);
VRA2(i,1)=100*prctile(tempVRA2(~isnan(tempVRA1)&criteria),50);
ArThres(i,1)=100*prctile(tempArThres(~isnan(tempArThres)&criteria),50);
ArThresvalues=100*tempArThres(~isnan(tempArThres)&criteria);
% Rsq_median(i,1)=prctile(Rsq(criteria),50)
LG0direct(i,1)=prctile(tempLG0(~isnan(tempLG0)&criteria),50);
tau(i,1)=prctile(temptau(~isnan(temptau)&criteria),50);
delay(i,1)=prctile(tempdelay(~isnan(tempdelay)&criteria),50);
delayvalues=tempdelay(~isnan(tempdelay)&criteria);
LG3min(i,1)=prctile(tempLG3min(~isnan(tempLG3min)&criteria),50);
LG6min(i,1)=prctile(tempLG6min(~isnan(tempLG6min)&criteria),50);
LG90s(i,1)=prctile(tempLG90s(~isnan(tempLG90s)&criteria),50);

%limitVRA1=1.59;
%VRA1_boundaryN(1,i,1)=sum(tempVRA1(~isnan(tempVRA1)&criteria)>limitVRA1);

%%
Fstates1=[mean(Fnrem1(criteria)) mean(Fnrem2(criteria)) mean(Fnrem3(criteria)) mean(tempFREM(criteria))]';
Fstates2=[mean(Fnrem1(criteria&tempN_arousals>0)) mean(Fnrem2(criteria&tempN_arousals>0)) mean(Fnrem3(criteria&tempN_arousals>0)) mean(tempFREM(criteria&tempN_arousals>0))]';


%% Plots LG, ArThres
if settings.plotfigs
    figure(1);
    
    subplot(2,6,3);
    xedges=0:0.1:1.2;
    h=histogram(LG1values,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('LG1')
    
    subplot(2,6,4);
    xedges=0:0.1:1.2;
    h=histogram(LGnvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('LGn');
    
    subplot(2,6,5);
    xedges=0:2:30;
    h=histogram(delayvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('delay')
    
    subplot(2,6,6);
    xedges=0:10:100;
    h=histogram(VRA1values,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('VRA')
    
    subplot(2,6,9);
    xedges=(0.8:0.1:3)*100;
    h=histogram(ArThresvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('ArThres')
end

%get DataN
for i=1:Nperwindowvars
    eval([varlist{i} 'N=length(' varlist{i} 'values);']);
end

if getCIs
    Nbootstrap=50;
    %default
    for i=1:length(varlist)
        eval([varlist{i} 'CI=[-Inf;Inf];']);
    end
    for i=1:Nperwindowvars
        try
            eval([varlist{i} 'CI=bootci(Nbootstrap,@median,' varlist{i} 'values);']);
        catch me
        end
    end
    
    if settings.plotfigs
        i=1;
        subplot(2,6,3);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=2;
        subplot(2,6,4);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=3;
        subplot(2,6,5);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=4;
        subplot(2,6,6);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=5;
        subplot(2,6,9);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    end
else
    Nbootstrap=0;
end
%% Summary
%ArThres(ArThres==0)=NaN;
ArThresPSG1=ArThres;
%% UA phenotype using model drive

unnormalize = 0;
subplot(1,3,1);
fontsize_=10;

minNevents=1;
maxwakethres=300;
Nciles=10;
n=1;

N_events = LG_QualityInfo{n}(:,2);
FREM = SleepData{n}(:,6);
minwake = SleepData{n}(:,7);
criteria = Poscriteria==1 & nremXcriteria==1 & TimeOfNightCriteria==1 & N_events>=minNevents & minwake<=maxwakethres & FREM<=maxFREM & FREM>=minFrem;
NwindowsForUA=sum(criteria);
for i=(Nperwindowvars+1):length(varlist)
    eval([varlist{i} 'N=NwindowsForUA;']);
end

Fstates3=[mean(Fnrem1(criteria)) mean(Fnrem2(criteria)) mean(Fnrem3(criteria)) mean(tempFREM(criteria))]';
% 
% if NwindowsForUA>0
%     if settings.plotfigs
%         plot(1,1,'marker','none'); hold('on');
%     end
%     hold('on');
%     [t x y a win nota veup hyp ttot] = VEVdriveArray(DataOut,n,criteria);
%     
%     x_=1;
%     if unnormalize
%         y = y.*veup;
%         x = x.*veup;
%         ymean = mean(y); %re-normalize with mean over whole night
%         y = y/ymean;
%         x = x/ymean;
%     end
%     
%     arthres(n) = ArThresPSG1(n)/100; % taken from window summary data, median of means
%     
%     if 1 %set minimum for arthres PER SE
%         if arthres(n)<1.00
%             arthres(n)=1.00;
%         end
%     end
%     arthresPSG2(n) = arthres(n);
%     if arthresPSG2(n)<1.05
%         arthresPSG2(n)=1.05;
%     end
%     
%     if scoredarousalsinwake
%         x = x(nota==1);
%         y = y(nota==1);
%     else
%         x = x(nota==1&hyp~=4);
%         y = y(nota==1&hyp~=4);
%     end
%     
%     
%     ploton=settings.plotfigs;
%     if settings.plotfigs
%         plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
%         hold('on');
%         plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
%     end
%     [Vpassive,Vactive,~,~,VpassiveCI,VactiveCI]=VEVdriveAnalysis(100*x,100*y,100*arthresPSG2(n),ploton,100,Nciles,Nbootstrap);
%     VpassiveCI=VpassiveCI(:);
%     VactiveCI=VactiveCI(:);
%     if 0&settings.plotfigs
%         set(gca,'fontsize',fontsize_)
%         xlim([0 max([prctile(100*x,95) 300])]);
%         ylim([0 max([prctile(100*y,95) 100])]);
%         ylabel('Ventilation, %eupnea'); xlabel('Vdrive model, %eupnea');
%         hold('off');
%         title(['AHI=',num2str(round(AHIstate),'%u')]);
%         
%         subplot(2,6,10);
%         xedges=0:0.1:1.4;
%         dx=0.1;
%         h=histogram(100*y(x>(1-dx)&x<(1+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
%         box('off');
%         xlim(100*[min(xedges) max(xedges)]);
%         %ylim([0 max([0.2;h.Values(:)])]);
%         ylabel('Vpassive')
%         
%         subplot(2,6,11);
%         h=histogram(100*y(x>(arthresPSG2-dx)&x<(arthresPSG2+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
%         box('off');
%         xlim(100*[min(xedges) max(xedges)]);
%         %ylim([0 max([0.2;h.Values(:)])]);
%         ylabel('Vactive')
%     end
%     if 0&settings.plotfigs
%         i=6;
%         subplot(2,6,10);
%         ylims=get(gca,'YLim');
%         hold('on');
%         plot(eval(varlist{i})*[1 1],ylims,'k');
%         plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         
%         i=7;
%         subplot(2,6,11);
%         ylims=get(gca,'YLim');
%         hold('on');
%         plot(eval(varlist{i})*[1 1],ylims,'k');
%         plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%     end
%     Vcomp = Vactive-Vpassive;
%     VcompCI = [VactiveCI(1)-VpassiveCI(2);VactiveCI(2)-VpassiveCI(1)];
%     VcompCI = sort(VcompCI);
% else
%     Vpassive=NaN;
%     Vactive=NaN;
%     Vcomp=NaN;
%     VpassiveCI=[-Inf;Inf];
%     VactiveCI=[-Inf;Inf];
%     VcompCI=[-Inf;Inf];
% end
%%
if NwindowsForUA>0
    if settings.plotfigs
        plot(1,1,'marker','none'); hold('on');
    end
    hold('on');
drivecol=16;   %%[16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak 24:FlowEdi, 23:FlowPes]
dividedrivebyttot=0;
normalizeyusingveupwindow=1;
normalizeyusingknownVeupnea=0;
Gequals1=0;
Nrows=8;

Fvol=0; Fpcw=0;


N_events = LG_QualityInfo{n}(:,2);
Pos = LG_QualityInfo{n}(:,5);
Fwake = SleepData{n}(:,1);
FREM = SleepData{n}(:,6);
minwake = SleepData{n}(:,7);

%criteria = N_events>=minNevents & minwake<maxwakethres2 & FREM<=maxFREM ; %& (~isnan(Pos)&supine==1)
Nwindows(n)=sum(criteria);

[t,x,y,a,win,nota,veup,hyp,pcw,t2,e] = VEPesArray(DataOut,n,criteria,drivecol);

y_ = nanmedian(veup);
        y = y./y_*100;
y_n(n)=y_;

%get wake data next (27,23,17)
Wcriteria = minwake>=0; %(~isnan(Pos)&supine==1)
[t,x,y,a,win,nota,veup,hyp,pcw,t2] = VEPesArray(DataOut,n,Wcriteria,drivecol);
if dividedrivebyttot
    x=x./t2;
end
[ad]=howfarawayfromsleep(a,win); %replace a with (a==1)&(hyp==4) if arousals not scored in wake...
%find breaths during wakefulness and arousals
    minNwakebreaths = 50;
    hh=hist(ad,[1:11]); hh(end)=NaN; th=find(hh>minNwakebreaths,1,'last');
    threshold(n) = min([4 th]);
a1 = ad>threshold(n);
g = nanmedian(y(a1)./x(a1));
x_wake(n)=median(x(a1));
y_wake(n)=median(y(a1));
G_N(n) = length(y(a1)./x(a1));

%GET BREATH-BY-BREATH SLEEP DATA
[t,x,y,a,win,nota,veup,hyp,pcw,t2,e] = VEPesArray(DataOut,n,criteria,drivecol);
if dividedrivebyttot
    x=x./t2;
end
%normalize data
y = y./veup; %note Y is unnormalized VE from LGfromflow
y_ = nanmedian(veup);
%use y_ from first section
x_ = y_n(n)/g; %xeup
x_n(n)=x_;

%Arousal threshold

%First find arthres for Vactive/Vcomp:
settings2=[3 0 1 0 1 0 1 0.2];
%1:ignorefirstXbreathsaftersleeponset [4]
%2:swapbreathsforhigherdrive [N]
%3:Nbreathsattributabletoarousal [N]
%4:increasingdriveonly [N]
%5:deleteifbeloweupnea [Y]
%6:setaseupneaifbeloweupnea [N]
[ArThres(n),arthres_N(n),~,ArThresvalues] = ArThresNew(a',x',x_*1.05,win,settings2);
arthres2=ArThres;
ArThresCI=bootci(Nbootstrap,@median,ArThresvalues)/x_*100;   
    

if arthres2(n)<(1.05*x_),arthres2(n)=1.05*x_; end
if ArThres(n)<(1.00*x_),ArThres(n)=1.00*x_; end


%Remove wakefulness breaths for further analysis
if scoredarousalsinwake
    x = x(nota==1);
    y = y(nota==1);
    veup = veup(nota==1);
else
    x = x(nota==1&hyp~=4);
    y = y(nota==1&hyp~=4);
    veup = veup(nota==1&hyp~=4);
end

x=x/x_;
ArThres=ArThres/x_*100;
arthres2=arthres2/x_;

veonvdrive=y./x*100;
veonvwake=y.*veup./y_wake*100;
driveondrivewake=x.*x_/x_wake*100;


ploton=settings.plotfigs;
if settings.plotfigs
    plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
    hold('on');
    plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
end
[Vpassive,Vactive,~,~,VpassiveCI,VactiveCI]=VEVdriveAnalysis(100*x,100*y,100*arthres2(n),ploton,100,Nciles,Nbootstrap);
VpassiveCI=VpassiveCI(:);
VactiveCI=VactiveCI(:);
if settings.plotfigs
    set(gca,'fontsize',fontsize_)
    xlim([0 max([prctile(100*x,95) 300])]);
    ylim([0 max([prctile(100*y,95) 100])]);
    ylabel('Ventilation, %eupnea'); xlabel('Vdrive model, %eupnea');
    hold('off');
    title(['AHI=',num2str(round(AHIstate),'%u')]);
    
    subplot(2,6,9);
    xedges3=0.9:0.1:3;
    
    dx=0.1;
    h=histogram(100*ArThresvalues/x_,100*xedges3,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim(100*[min(xedges3) max(xedges3)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('ArThres')
    
    subplot(2,6,10);
    xedges=0:0.1:1.4;
    xedges2=0:0.1:4.0;
    dx=0.1;
    h=histogram(100*y(x>(1-dx)&x<(1+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim(100*[min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Vpassive')
    
    subplot(2,6,11);
    h=histogram(100*y(x>(arthres2-dx)&x<(arthres2+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim(100*[min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Vactive')
    
    
end
if settings.plotfigs
    
    i=5;
    subplot(2,6,9);
    ylims=get(gca,'YLim');
    hold('on');
    plot(eval(varlist{i})*[1 1],ylims,'k');
    plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    
    i=6;
    subplot(2,6,10);
    ylims=get(gca,'YLim');
    hold('on');
    plot(eval(varlist{i})*[1 1],ylims,'k');
    plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    
    i=7;
    subplot(2,6,11);
    ylims=get(gca,'YLim');
    hold('on');
    plot(eval(varlist{i})*[1 1],ylims,'k');
    plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    
%     subplot(2,6,12);
%     ylims=get(gca,'YLim');
%     hold('on');
%     plot(eval(varlist{i})*[1 1],ylims,'k');
%     plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%     plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
end
Vcomp = Vactive-Vpassive;
VcompCI = [VactiveCI(1)-VpassiveCI(2);VactiveCI(2)-VpassiveCI(1)];
VcompCI = sort(VcompCI);

    figure(2); clf(2);
    varlist2 = {'veonvdrive','veonvwake','driveondrivewake'};
    xlimupper = [140 140 300];
    for i=1:length(varlist2)
        h=histogram(eval(varlist2{i}),0:10:xlimupper(i),'normalization','probability','EdgeAlpha',0);
        hold('on');
        eval([varlist2{i} '(isnan(' varlist2{i} '))=[];']);
        eval([varlist2{i} 'CI=bootci(Nbootstrap,@median,' varlist2{i} ');']);
    end
    box('off');
    xlim([0 max(xlimupper)]);
    
    xlabel('Value, %');
    ylabel('Proportion of Sleep');
    
    ylims=get(gca,'YLim');
    i=1;
	plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0 0.45 0.74]);
    i=2;
    plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0.85 0.33 0.1]);
    i=3;
    plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0.93 0.69 0.13]);
    
    h=legend('Actual Airflow / Intended Airflow','Actual Airflow / Wake Baseline','Intended Airflow / Wake Baseline')
    
else
    Vpassive=NaN;
    Vactive=NaN;
    Vcomp=NaN;
    VpassiveCI=[-Inf;Inf];
    VactiveCI=[-Inf;Inf];
    VcompCI=[-Inf;Inf];
end
%%
disp([x_ y_ x_wake y_wake])

%% Save
criteriaset=[1 1 1 2 2 3 3 3];
clear Data DataCI Fstates
for i=1:length(varlist)
    Data(i) = eval(varlist{i});
    if settings.getCIs
        DataCI(:,i) = eval([varlist{i} 'CI']);
    else
        DataCI(:,i) = NaN;
    end
    DataN(i) = eval([varlist{i} 'N']);
    Fstates(:,i) = eval(['Fstates' num2str(criteriaset(i))]);
end

