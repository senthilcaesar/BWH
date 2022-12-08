function [Data,DataCI,DataN,varlist,AHItotal,Fstates] = SummaryAnalysisOne_DLM(subj,settings)

%% Hardcoded directory info
%directory=settings.directory;%'H:\MESA\Analyzed3';
savename =settings.savename;%'MESA2Sept2017';

%% added by DLM
directory=settings.OutputDataDirectory;
settings.getCIs=1;
settings.plotfigs=1;
settings.selecttimeofnight=0;
settings.selecttimeofnight_XTiles=2;
settings.selecttimeofnight_NthXTile=2;
settings.selectstate=4; %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
settings.selectposition=[]; %empty for unused, %1=supine in profusion

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
    ylabel('LG1')
    
    subplot(2,6,4);
    xedges=0:0.1:1.2;
    h=histogram(LGnvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    ylabel('LGn');
    
    subplot(2,6,5);
    xedges=0:2:30;
    h=histogram(delayvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    ylabel('delay')
    
    subplot(2,6,6);
    xedges=0:10:100;
    h=histogram(VRA1values,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    ylabel('VRA')
    
    subplot(2,6,9);
    xedges=(0.8:0.1:3)*100;
    h=histogram(ArThresvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    ylabel('ArThres')
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
Exclude=0;

unnormalize = 0;
subplot(1,3,1);
fontsize_=10;

minNwindows = 3; %unused
% if state==1
%     maxFREM=0;
%     minFrem=-Inf;
% elseif state==2
%     maxFREM=Inf;
%     minFrem=0.5;
% end

minNevents=1;
maxwakethres=300;
Nciles=10;
n=1;
%for n=1:M
%     if Exclude(n)==1||LG1_N(n)<=minNwindows
%         continue
%     end
    
    %title(filename_{n});
    
    N_events = LG_QualityInfo{n}(:,2);
    
    FREM = SleepData{n}(:,6);
    minwake = SleepData{n}(:,7);
    %criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=maxwakethres;
    criteria = Poscriteria==1 & nremXcriteria==1 & TimeOfNightCriteria==1 & N_events>=minNevents & minwake<=maxwakethres & FREM<=maxFREM & FREM>=minFrem;
    NwindowsForUA=sum(criteria);
    for i=(Nperwindowvars+1):length(varlist)
        eval([varlist{i} 'N=NwindowsForUA;']);
    end
    Fstates3=[mean(Fnrem1(criteria)) mean(Fnrem2(criteria)) mean(Fnrem3(criteria)) mean(tempFREM(criteria))]';

%     if sum(criteria)==0

%         continue
%     end
%get DataN
if NwindowsForUA>0
    
    if settings.plotfigs
    plot(1,1,'marker','none'); hold('on');
    end
    %title(filename_{n});
    hold('on');
    [t x y a win nota veup hyp ttot] = VEVdriveArray(DataOut,n,criteria);
    
    
    x_=1;
    if unnormalize
        y = y.*veup;
        x = x.*veup;
        ymean = mean(y); %re-normalize with mean over whole night
        y = y/ymean;
        x = x/ymean;
    end
    %y=y*100;
    
    %         if kk==1
    %             data_B{n}{iii}=[win t y ttot];
    %         else
    %             data_B_flow{n}=[win t y ttot];
    %         end
    
    %arthres(n) = ArThresFromPSG(1,x,win);
    %increasingdriveonly, deleteifbeloweupnea(0),
    %setaseupneaifbeloweupnea(1), usemediannotROC(1), falsepositiverate(0.2)
    if 0
        settings=[1 0 1 0 1 0 1 0.2];
        %ignorefirstXbreathsaftersleeponset, swapbreathsforhigherdrive, Nbreathsattributabletoarousal, ...
        %increasingdriveonly, deleteifbeloweupnea, setaseupneaifbeloweupnea
        %[arthres(n)] = ArThresNew(a',x',x_,win,settings);
        [arthres(n)] = ArThresNew(a',[x(2:end);NaN]',x_,win,settings);
    elseif 1
        %arthres(n) = ArThresFromPSG(1,x,win);
        arthres(n) = ArThresPSG1(n)/100; % taken from window summary data, median of means
    elseif 0
        ArOnset = [NaN;diff(a)];
        arthres(n)= median(x(ArOnset==1));
    else %window by window average as per LGfromFlow (near exact)
        temparthres=NaN*ones(max(win),1);
        ArOnset = [NaN;diff(a)];
        for i=1:max(win)
            temparthres(i)=nanmean(x(ArOnset==1&win==i));
        end
        arthres(n)=nanmedian(temparthres);
    end
    if 0
        figure(34)
        plot(ArThresPSG1/100,arthres,'.');
    end
    if 1 %set minimum for arthres PER SE
        if arthres(n)<1.00
            arthres(n)=1.00;
        end
    end
    arthresPSG2(n) = arthres(n);
    if arthresPSG2(n)<1.05
        arthresPSG2(n)=1.05;
    end
    
    
    %arthres(n) = 1.6;
    if 0
        x = x(nota==1);
        y = y(nota==1);
    else
        x = x(nota==1&hyp~=4);
        y = y(nota==1&hyp~=4);
    end
    
    
    medianV(n) = median(y);
    
    ploton=settings.plotfigs;
     if settings.plotfigs
    plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
    hold('on');
    plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
     end
    [Vpassive,Vactive,~,~,VpassiveCI,VactiveCI]=VEVdriveAnalysis(100*x,100*y,100*arthresPSG2(n),ploton,100,Nciles,Nbootstrap);
    VpassiveCI=VpassiveCI(:);
    VactiveCI=VactiveCI(:);
    if settings.plotfigs
    set(gca,'fontsize',fontsize_)
    xlim([0 max([prctile(100*x,95) 300])]);
    ylim([0 max([prctile(100*y,95) 100])]);
    ylabel('Ventilation, %eupnea'); xlabel('Vdrive model, %eupnea');
    hold('off');
    title(['AHI=',num2str(round(AHIstate),'%u')]);
    
    subplot(2,6,10);
    xedges=0:0.1:1.4;
    dx=0.1;
    h=histogram(100*y(x>(1-dx)&x<(1+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim(100*[min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    ylabel('Vpassive')
    
    subplot(2,6,11);
    h=histogram(100*y(x>(arthresPSG2-dx)&x<(arthresPSG2+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim(100*[min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    ylabel('Vactive')
    end
    if settings.plotfigs
        
        
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
        
    end
    
%end

Vcomp = Vactive-Vpassive;
VcompCI = [VactiveCI(1)-VpassiveCI(2);VactiveCI(2)-VpassiveCI(1)];
    VcompCI = sort(VcompCI);

else
    Vpassive=NaN;
    Vactive=NaN;
    Vcomp=NaN;
    VpassiveCI=[-Inf;Inf];
    VactiveCI=[-Inf;Inf];
    VcompCI=[-Inf;Inf];
end
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

