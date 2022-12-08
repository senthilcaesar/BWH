% Produce group level endotype plots
% run stats for centiles (turned off)

%clear all;
fixedxupper=260;
%note 1st value contains 0-5th centile, centered on 2.5th centile.

%% Load
[file,dir] = uigetfile();
filedir = [dir, file];
matObj = matfile(filedir);
varlist = who(matObj);
load(filedir);


%% Select data comparison
% 1 - A vs B
% 2 - A vs B1
% 3 - A vs B2
% 4 - B1 vs B2
% 5 - All subjects rising/falling
% 6 - A rising/falling
% 7 - B rising/falling
% 8 - B1 rising/falling
% 9 - B2 rising/falling
%10 - n1 n2 n3
%-1 - single subject
%-2 - all subjects (set color2 blank)

comparisonMode = 5
N=19
holdfigs=10
redOption=[ 0.95 0.2 0.05]
blueOption= [ 0.1 0.4 0.7]
purpleOption= [0.3 0.1 0.4]
blank = [1 1 1]

color1 =  redOption
color2 = blueOption %blueOption %compare group
% 
% color1 =  blank
% color2 = purpleOption %blu

displayColorPlot=1;plotColorBar=0
ignoreGGdotsandlines=0;
ignoreVentdotsandlines=1;
%%
[~,~,raw]=xlsread('E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018\PUPstart\AMasterSpreadsheet','BI4:BI63');
raw = string(raw);
P = NaN*ones(60,1);
P(raw=="A")=1;
P(raw=="B1")=0;
P(raw=="B2")=-1;
P(raw=="B")=-2;

if comparisonMode == 5
[file,dir] = uigetfile();
filedir = [dir, file];
matObj = matfile(filedir);
varlist = who(matObj);
newvar = load(filedir,'DataNArray', 'DataNArray_compare')
DataNArray2=newvar.DataNArray;
DataNArray_compare2=newvar.DataNArray_compare;
end
%%
usefilter121=0; %
if usefilter121>0
    for i=1:size(VEdeciles,2)
        VEdeciles(:,i) = filter121(VEdeciles(:,i),usefilter121);
    end
    for i=1:size(GGpdeciles,2)
        GGpdeciles(:,i) = filter121(GGpdeciles(:,i),usefilter121);
    end
    for i=1:size(GGtdeciles,2)
        GGtdeciles(:,i) = filter121(GGtdeciles(:,i),usefilter121);
    end
    
    try
        for i=1:size(VEdeciles_compare,2)
            VEdeciles_compare(:,i) = filter121(VEdeciles_compare(:,i),usefilter121);
        end
        for i=1:size(GGpdeciles_compare,2)
            GGpdeciles_compare(:,i) = filter121(GGpdeciles_compare(:,i),usefilter121);
        end
        for i=1:size(GGtdeciles_compare,2)
            GGtdeciles_compare(:,i) = filter121(GGtdeciles_compare(:,i),usefilter121);
        end
    catch me
    end
end

DataNArray=DataNArray_compare;
    DataArray=DataArray_compare;
    GGdata=GGdata_compare;
    GGpdeciles=GGpdeciles_compare;
    GGtdeciles=GGtdeciles_compare;
    Vdrivedeciles=Vdrivedeciles_compare;
    VEdeciles=VEdeciles_compare;
    Veupnea=Veupnea_compare;
    EdecilesMean=EdecilesMean_compare;


if comparisonMode<5 
    DataNArray_compare=DataNArray;
    DataArray_compare=DataArray;
    GGdata_compare=GGdata;
    GGpdeciles_compare=GGpdeciles;
    GGtdeciles_compare=GGtdeciles;
    Vdrivedeciles_compare=Vdrivedeciles;
    VEdeciles_compare=VEdeciles;
    Veupnea_compare=Veupnea;
    EdecilesMean_compare=EdecilesMean;
else
    %load('C:\Users\lg373\Dropbox (Partners HealthCare)\PhenotypeDrive2018\Summary\RisingFalling\SummaryAnalysis.mat')
end

settings.runcompare=1;
minNthres=1;
setpositionflowdrive = [880   346   363   351];
setpositionGG = [786   474   247   266];



%% FlowVsDrive PhenoPlot Average
if settings.runcompare
    % A - exclude others (P~=1)
    % B - exclude isnan and A: (P==1) | isnan(P)
    % B1 - GG nonresponsive: exclude others (P~= 0)
    % B2 - GG responsive / ineffective: (P~=-1)
    % Both B => exclude (P~=-1) & (P~= 0)
    % Compare=blue
    %Iexclude_compare = DataNArray_compare(:,6)<2 | DataNArray(:,6)<2;
    %Iexclude = DataNArray_compare(:,6)<2 | DataNArray(:,6)<2;
    
    %Initialise data to A
    Iexclude = DataNArray(:,1)<minNthres | (P~=1)
    
    switch comparisonMode
        
        case 1  %A vs B all:
            Iexclude =  DataNArray(:,1)<minNthres |(P~=1) %A (red)
            Iexclude_compare = DataNArray_compare(:,1)<minNthres |(P==1) | isnan(P) %B (blue)
            
        case 2 %A vs B1 nonresponsive GG
            Iexclude = DataNArray(:,1)<minNthres | (P~=1)
            Iexclude_compare = DataNArray_compare(:,1)<minNthres | (P~= 0)
            
        case 3 %A vs B2 ineffective GG (terrible airway?)
            Iexclude = DataNArray(:,1)<minNthres | (P~=1)
            Iexclude_compare = DataNArray_compare(:,1)<minNthres | (P~=-1)
            
        case 4 %Compare B1 and B2
            Iexclude = DataNArray_compare(:,1)<minNthres | (P~= 0)
            Iexclude_compare = DataNArray_compare(:,1)<minNthres | (P~=-1)
            
        case 5 %All rising and falling
            Iexclude = DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | isnan(P)
            Iexclude_compare = Iexclude;
            
        case 6 %A rising and falling
            Iexclude = DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | (P~=1)
            Iexclude_compare = Iexclude;
            
        case 7 %B rising and falling
            Iexclude = DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | (P==1) | isnan(P);
            Iexclude_compare = Iexclude;
            
        case 8  %B1 rising and falling
            Iexclude = DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | (P~= 0)
            Iexclude_compare = Iexclude;
            
        case 9  %B2 rising and falling
            Iexclude = DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | (P~= -1)
            Iexclude_compare = Iexclude;
          case 10  %Combine N2 N3
            Iinclude = (DataNArray(:,1)>minNthres & DataNArray_compare(:,1)>minNthres) | (DataNArray2(:,1)>minNthres & DataNArray_compare2(:,1)>minNthres);
            Iexclude = Iinclude==0
          %  Iexclude = (DataNArray(:,1)<minNthres & DataNArray2(:,1)<minNthres) | (DataNArray_compare(:,1)<minNthres;
           %  Iexclude = (DataNArray(:,1)<minNthres) | DataNArray_compare(:,1)<minNthres;
            Iexclude_compare = DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres;
              
        case -11  %Combine N2 N3
            Iinclude = zeros(60,1);
            Iinclude(T4.Subject(T4.GGphenotype==-1))=1 %(red)
            Iexclude = Iinclude==0
            
            Iinclude = zeros(60,1);
            Iinclude(T4.Subject(T4.GGphenotype==1))=1 %(blue)
            Iexclude_compare = Iinclude==0
             
       
        case -1
                 Iexclude = (1:60)~=N %DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | isnan(P)
    Iexclude_compare = Iexclude
        case -2
                 Iexclude = DataNArray(:,1)<minNthres |  isnan(P); %DataNArray(:,1)<minNthres | DataNArray_compare(:,1)<minNthres | isnan(P)
    Iexclude_compare = Iexclude
    
    end 
    
    VEdeciles_compare(:,Iexclude_compare==1)=NaN;
    Vdrivedeciles_compare(:,Iexclude_compare==1)=NaN;
    DataArrayIncl_compare=DataArray_compare;
    DataArrayIncl_compare(Iexclude_compare==1,:)=NaN;
    GGpdeciles_compare(:,Iexclude_compare==1)=NaN;
    GGtdeciles_compare(:,Iexclude_compare==1)=NaN;
    EdecilesMean_compare(:,Iexclude_compare==1)=NaN;
else
    Iexclude = DataNArray(:,6)<2;
end

VEdeciles(:,Iexclude==1)=NaN;
Vdrivedeciles(:,Iexclude==1)=NaN;
DataArrayIncl=DataArray;
DataArrayIncl(Iexclude==1,:)=NaN;
GGpdeciles(:,Iexclude==1)=NaN;
GGtdeciles(:,Iexclude==1)=NaN;
GGdataIncl = GGdata;
GGdataIncl(Iexclude==1,:)=NaN;
EdecilesMean(:,Iexclude==1)=NaN;

GGdataIncl_compare = GGdata_compare;
GGdataIncl_compare(Iexclude_compare==1,:)=NaN;

GGpdecilesPmax = GGpdeciles;
GGtdecilesPmax = GGtdeciles;
GGpdeciles_comparePmax = GGpdeciles_compare;
GGtdeciles_comparePmax = GGtdeciles_compare;

if 0 %for plots, but doesn't break stats in mixed models (they use GGpdecilesPmax)
    GGpdeciles = 100*GGpdeciles./GGdataIncl(:,1)';
    GGtdeciles = 100*GGtdeciles./GGdataIncl(:,1)';
    GGpdeciles_compare = 100*GGpdeciles_compare./GGdataIncl_compare(:,1)';
    GGtdeciles_compare = 100*GGtdeciles_compare./GGdataIncl_compare(:,1)';
end
if 1 %plots only, don't use for stats
    GGpdeciles = GGpdeciles./GGdataIncl(:,1)'*nanmedian(GGdataIncl(:,1)');
    GGtdeciles = GGtdeciles./GGdataIncl(:,1)'*nanmedian(GGdataIncl(:,1)');
    GGpdeciles_compare = GGpdeciles_compare./GGdataIncl_compare(:,1)'*nanmedian(GGdataIncl_compare(:,1)');
    GGtdeciles_compare = GGtdeciles_compare./GGdataIncl_compare(:,1)'*nanmedian(GGdataIncl_compare(:,1)');
end

% normalisation required for between phenotype comparison, none for rising/falling
if 0 %comparisonMode < 5
    GGMeanFactor = nanmedian(GGdataIncl(:,1)')/nanmedian(GGdataIncl_compare(:,1)');
    %GGMeanFactor = nanmean(GGdataIncl(:,1)')/nanmean(GGdataIncl_compare(:,1)');
else
    GGMeanFactor = 1;
end

%%
upper95CI = @(x) nanmean(x) + 1.96*nanstd(x)./sum(~isnan(x)).^0.5;
lower95CI = @(x) nanmean(x) - 1.96*nanstd(x)./sum(~isnan(x)).^0.5;



upperIQR = @(x) prctile2(x,75);
lowerIQR = @(x) prctile2(x,25);


if 0
    upperlimit = upper95CI;
    lowerlimit = lower95CI;
    average = @(x) nanmean(x);
else
    upperlimit = upperIQR;
    lowerlimit = lowerIQR;
    average = @(x) nanmedian(x);
end


EdecilesMean_compare=1-EdecilesMean_compare;
EdecilesMean=1-EdecilesMean;
if settings.runcompare
    EdecilesMeanUpper_compare = upperlimit(EdecilesMean_compare')';
    EdecilesMeanLower_compare = lowerlimit(EdecilesMean_compare')';
    EdecilesMeanLower_compare(EdecilesMeanLower_compare<0) = 0;
    VdrivedecilesMedian_compare = average(Vdrivedeciles_compare')';
    EdecilesMeanMedian_compare = average(EdecilesMean_compare')';
end

 EdecilesMeanUpper = upperlimit(EdecilesMean')';
    EdecilesMeanLower = lowerlimit(EdecilesMean')';
    VdrivedecilesMedian = average(Vdrivedeciles')';
    EdecilesMeanMedian = average(EdecilesMean')';
    



%% GGp
figure(10); 
if holdfigs ~=1
clf(10);
end
%subplot(2,1,1);
if settings.runcompare
    GGpdecilesUpper_compare = upperlimit(GGpdeciles_compare')';
    GGpdecilesLower_compare = lowerlimit(GGpdeciles_compare')';
    VdrivedecilesMedian_compare = average(Vdrivedeciles_compare')';
    GGpdecilesMedian_compare = average(GGpdeciles_compare')';
    
    set(gcf,'color',[1 1 1]);
    if 1
        hold('on')
        box('off');
        set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[GGpdecilesUpper_compare;flipud(GGpdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,GGpdecilesMedian_compare,'k','linewidth',2);
   
    [~,ii] = unique(VdrivedecilesMedian_compare);
    GGppassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),GGpdecilesMedian_compare(ii),100,'linear');
    arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
   
    if ~ignoreGGdotsandlines
        GGpactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),GGpdecilesMedian_compare(ii),arthresMedian_compare,'linear');
        plot(100,GGppassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
        plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,GGpactiveMedian_compare)],'g'); %Arthres line
        %plot(arthresMedian_compare,GGpactiveMedian_compare,'r.','markersize',25); %plots Vactive dot if required
    end
end

GGpdecilesUpper = GGMeanFactor*upperlimit(GGpdeciles')';
GGpdecilesLower = GGMeanFactor*lowerlimit(GGpdeciles')';
VdrivedecilesMedian = average(Vdrivedeciles')';
GGpdecilesMedian = GGMeanFactor*average(GGpdeciles')';

if ~settings.runcompare
    plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
    hold('on')
    %plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
end

fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGpdecilesUpper;flipud(GGpdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,GGpdecilesMedian,'k','linewidth',2);
if displayColorPlot
DispEventLikelihoodOnEndogram(VdrivedecilesMedian,GGpdecilesMedian,EdecilesMeanMedian,plotColorBar)
DispEventLikelihoodOnEndogram(VdrivedecilesMedian_compare,GGpdecilesMedian_compare,EdecilesMeanMedian_compare,plotColorBar)
end

%     VpassiveMedian = nanmedian(DataArrayIncl(:,6));
arthresMedian = nanmedian(DataArrayIncl(:,5));
%     VactiveMedian = nanmedian(DataArrayIncl(:,7));

[~,ii] = unique(VdrivedecilesMedian);
GGppassiveMedian=interp1(VdrivedecilesMedian(ii),GGpdecilesMedian(ii),100,'linear');
GGpactiveMedian = interp1(VdrivedecilesMedian(ii),GGpdecilesMedian(ii),arthresMedian,'linear');


if ~ignoreGGdotsandlines
    plot(100,GGppassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian arthresMedian],[0 GGpactiveMedian],'g'); %Arthres line
   % plot(arthresMedian,GGpactiveMedian,'r.','markersize',25);  %Vactive
end
%saveas(10,[settings.directoryout 'Fig9_' settings.Comparename '_' num2str(i) '_s' num2str(settings.selectstate) '.png'])

if settings.runcompare
    ylim([0 max(max(GGpdecilesUpper),max(GGpdecilesUpper_compare))])
else
    ylim([0 (max(GGpdecilesUpper))])
end

%subplot(2,1,2);

if settings.runcompare
    GGtdecilesUpper_compare = upperlimit(GGtdeciles_compare')';
    GGtdecilesLower_compare = lowerlimit(GGtdeciles_compare')';
    VdrivedecilesMedian_compare = average(Vdrivedeciles_compare')';
    GGtdecilesMedian_compare = average(GGtdeciles_compare')';
    
    set(gcf,'color',[1 1 1]);
    if 0
        plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
        hold('on')
        %plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
        box('off');
        set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[GGtdecilesUpper_compare;flipud(GGtdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,GGtdecilesMedian_compare,'k','linewidth',2);
    
    %VpassiveMedian_compare = nanmedian(DataArrayIncl_compare(:,6));
    
    [~,ii] = unique(VdrivedecilesMedian_compare);
    GGtpassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),GGtdecilesMedian_compare(ii),100,'linear');
    
    arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
    %     VactiveMedian_compare = nanmedian(DataArrayIncl_compare(:,7));
    if ~ignoreGGdotsandlines
        GGtactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),GGtdecilesMedian_compare(ii),arthresMedian_compare,'linear');
        plot(100,GGtpassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
        plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,GGtactiveMedian_compare)],'g'); %Arthres line
       % plot(arthresMedian_compare,GGtactiveMedian_compare,'r.','markersize',25);  %Vactive
    end
    
end


GGtdecilesUpper = GGMeanFactor*upperlimit(GGtdeciles')';
GGtdecilesLower = GGMeanFactor*lowerlimit(GGtdeciles')';

VdrivedecilesMedian = average(Vdrivedeciles')';
GGtdecilesMedian = GGMeanFactor*average(GGtdeciles')';

if ~settings.runcompare
    %plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
    hold('on')
    %plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
end
fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGtdecilesUpper;flipud(GGtdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,GGtdecilesMedian,'k','linewidth',2);
if displayColorPlot
DispEventLikelihoodOnEndogram(VdrivedecilesMedian,GGtdecilesMedian,EdecilesMeanMedian,plotColorBar)
    DispEventLikelihoodOnEndogram(VdrivedecilesMedian_compare,GGtdecilesMedian_compare,EdecilesMeanMedian_compare,plotColorBar)
end
%     VpassiveMedian = nanmedian(DataArrayIncl(:,6));
arthresMedian = nanmedian(DataArrayIncl(:,5));
%     VactiveMedian = nanmedian(DataArrayIncl(:,7));

[~,ii] = unique(VdrivedecilesMedian);
GGtpassiveMedian=interp1(VdrivedecilesMedian(ii),GGtdecilesMedian(ii),100,'linear');
GGtactiveMedian = interp1(VdrivedecilesMedian(ii),GGtdecilesMedian(ii),arthresMedian,'linear');


if ~ignoreGGdotsandlines
    plot(100,GGtpassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian arthresMedian],[0 GGtactiveMedian],'g'); %Arthres line
   % plot(arthresMedian,GGtactiveMedian,'r.','markersize',25);  %Vactive
end

if settings.runcompare
    ylim([0 max(max(GGpdecilesUpper),max(GGpdecilesUpper_compare))])
else
    ylim([0 (max(GGpdecilesUpper))])
end

xlim([0 fixedxupper]);
set(gcf,'position',setpositionGG);

%saveas(10,[settings.directoryout 'Fig10_' settings.Comparename '_s' num2str(settings.selectstate) '.png'])

%% Individual plots for debug
if 0
    
    figure(799); clf(799);
    for i=1:size(Vdrivedeciles_compare,2)
        plot(Vdrivedeciles_compare(:,i),GGpdeciles_compare(:,i));
        hold on
        i
        pause
    end
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[GGpdecilesUpper_compare;flipud(GGpdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    
    figure(789); clf(789);
    for i=1:size(Vdrivedeciles,2)
        plot(Vdrivedeciles(:,i),GGpdeciles(:,i));
        hold on
        i
        pause
    end
    fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGpdecilesUpper;flipud(GGpdecilesLower)],color2,'edgecolor','none','facealpha',0.5);
end



%%    Ventilation
figure(9); 
if holdfigs ~=1
clf(9);
end

if settings.runcompare
    VEdecilesUpper_compare = upperlimit(VEdeciles_compare')';
    VEdecilesLower_compare = lowerlimit(VEdeciles_compare')';
    VdrivedecilesMedian_compare = average(Vdrivedeciles_compare')';
    VEdecilesMedian_compare = average(VEdeciles_compare')';
    
    set(gcf,'color',[1 1 1]);
    hold('on')
    box('off');
    set(gca,'tickdir','out')
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[VEdecilesUpper_compare;flipud(VEdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,VEdecilesMedian_compare,'k','linewidth',2);
    
    %VpassiveMedian_compare = nanmedian(DataArrayIncl_compare(:,6));
    [~,ii] = unique(VdrivedecilesMedian_compare);
    VpassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),100,'linear');
    
    arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
    VactiveMedian_compare = nanmedian(DataArrayIncl_compare(:,7));
    
    VactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),arthresMedian_compare,'linear');
    if ~ignoreVentdotsandlines
    plot(100,VpassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,VactiveMedian_compare)],'g'); %Arthres line
   % plot(arthresMedian_compare,VactiveMedian_compare,'r.','markersize',25);  %Vactive
     if 1
        plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
        plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
        
     end
    end
   
    
end

VEdecilesUpper = upperlimit(VEdeciles')';
VEdecilesLower = lowerlimit(VEdeciles')';
VdrivedecilesMedian = average(Vdrivedeciles')';
VEdecilesMedian = average(VEdeciles')';

if ~settings.runcompare
    plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
    hold('on')
    plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
end
fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[VEdecilesUpper;flipud(VEdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,VEdecilesMedian,'k','linewidth',2);
 
if displayColorPlot
DispEventLikelihoodOnEndogram(VdrivedecilesMedian,VEdecilesMedian,EdecilesMeanMedian,plotColorBar)
DispEventLikelihoodOnEndogram(VdrivedecilesMedian_compare,VEdecilesMedian_compare,EdecilesMeanMedian_compare,plotColorBar)
end

VpassiveMedian = nanmedian(DataArrayIncl(:,6));
arthresMedian = nanmedian(DataArrayIncl(:,5));
VactiveMedian = nanmedian(DataArrayIncl(:,7));

[~,ii] = unique(VdrivedecilesMedian);
VpassiveMedian=interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),100,'linear');
VactiveMedian = interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),arthresMedian,'linear');


   if ~ignoreVentdotsandlines
plot(100,VpassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
plot([arthresMedian arthresMedian],[0 max(arthresMedian,VactiveMedian)],'g'); %Arthres line
%plot(arthresMedian,VactiveMedian,'r.','markersize',25);  %Vactive
   end
   
if settings.runcompare
    ylim([0 max(max(VEdecilesUpper),max(VEdecilesUpper_compare))])
else
    ylim([0 (max(VEdecilesUpper))])
end


xlim([0 fixedxupper]);

set(gcf,'position',setpositionflowdrive);

%saveas(9,[settings.directoryout 'Fig9_' settings.Comparename '_s' num2str(settings.selectstate) '.png'])

%%
if 1
figure(12); clf(12);
%EdecilesMean=TdecilesMean;
%EdecilesMean_compare=TdecilesMean_compare;
if settings.runcompare
    EdecilesMeanUpper_compare = upperlimit(EdecilesMean_compare')';
    EdecilesMeanLower_compare = lowerlimit(EdecilesMean_compare')';
    EdecilesMeanLower_compare(EdecilesMeanLower_compare<0) = 0;
    VdrivedecilesMedian_compare = average(Vdrivedeciles_compare')';
    EdecilesMeanMedian_compare = average(EdecilesMean_compare')';
    
    set(gcf,'color',[1 1 1]);
    if 1
        plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
        hold('on')
        plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
        box('off');
        set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[EdecilesMeanUpper_compare;flipud(EdecilesMeanLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,EdecilesMeanMedian_compare,'k','linewidth',2);
    
%    
%     [~,ii] = unique(VdrivedecilesMedian_compare);
%     VpassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),100,'linear');
%     
%     arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
%     VactiveMedian_compare = nanmedian(DataArrayIncl_compare(:,7));
%     
%     VactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),arthresMedian_compare,'linear');
%     plot(100,VpassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
%     plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,VactiveMedian_compare)],'g'); %Arthres line
%     plot(arthresMedian_compare,VactiveMedian_compare,'r.','markersize',25);  %Vactive
%     
%     
%     
end

 EdecilesMeanUpper = upperlimit(EdecilesMean')';
    EdecilesMeanLower = lowerlimit(EdecilesMean')';
    VdrivedecilesMedian = average(Vdrivedeciles')';
    EdecilesMeanMedian = average(EdecilesMean')';
    

if ~settings.runcompare
    plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
    hold('on')
    plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
end
hold on
fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[EdecilesMeanUpper;flipud(EdecilesMeanLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,EdecilesMeanMedian,'k','linewidth',2);
% 
% [~,ii] = unique(VdrivedecilesMedian);
% VpassiveMedian=interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),100,'linear');
% VactiveMedian = interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),arthresMedian,'linear');
% 
% 
% 
% plot(100,VpassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
% plot([arthresMedian arthresMedian],[0 max(arthresMedian,VactiveMedian)],'g'); %Arthres line
% plot(arthresMedian,VactiveMedian,'r.','markersize',25);  %Vactive

if settings.runcompare
    ylim([0 max(max(EdecilesMeanUpper),max(EdecilesMeanUpper_compare))])
else
    ylim([0 (max(EdecilesMeanUpper))])
end


xlim([0 fixedxupper]);

set(gcf,'position',setpositionflowdrive);

%saveas(9,[settings.directoryout 'Fig9_' settings.Comparename '_s' num2str(settings.selectstate) '.png'])
end
%%
if 1
%     rowLow=1;
%     rowMid=6;
%     rowHigh=10;
%     
%     temp=Vdrivedeciles(rowLow,:)
%     % LowestDriveA = [Vdrivedeciles(rowLow,P==1)]';
%     % LowestDriveB1 = [Vdrivedeciles(rowLow,P==0)]';
%     % LowestDriveB2 = [Vdrivedeciles(rowLow,P==-1)]';
%     [p2]=ranksum(LowestDriveA,[LowestDriveB1;LowestDriveB2])
%     for i=1:-1:-2
%         out(2-i,:)= [nanmean(temp(:,P==i)) nanstd(temp(:,P==i)) nanmedian(temp(:,P==i)) tsnaniqr(temp(:,P==i))]
%     end
%     [p2]=ranksum(temp(:,P==1),temp(:,(P==0|P==-1|P==-2)))
%     
%     
%     for i=1:-1:-2
%         out(2-i,:)= [nanmean(temp(P==i)) nanstd(temp(P==i)) nanmedian(temp(P==i)) tsnaniqr(temp(P==i))]
%     end
%     [p2]=ranksum(temp(P==1),temp((P==0|P==-1|P==-2)))
%     [nanmean(temp(P==1)) nanstd(temp(P==1))]
%     [nanmean(temp((P==0|P==-1|P==-2)))   nanstd(temp(P==0|P==-1|P==-2))]
%     
    
    %% ORIGINAL stats
     rowLow=1;
     rowMid=6;
     rowHigh=10;
    
    LowestDrive = [Vdrivedeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    LowestDrive_Compare = [Vdrivedeciles_compare(rowLow,:)]';
    LowestDrives = [LowestDrive LowestDrive_Compare];
    [~,p]=ttest2(LowestDrive,LowestDrive_Compare);
    [p2]=ranksum(LowestDrive,LowestDrive_Compare)
    
    [nanmean(LowestDrives) nanstd(LowestDrives) nanmedian(LowestDrives) tsnaniqr(LowestDrives)]
    
    MedianDrive = [Vdrivedeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    MedianDrive_Compare = [Vdrivedeciles_compare(rowMid,:)]';
    MedianDrives = [MedianDrive MedianDrive_Compare];
    [~,p]=ttest2(MedianDrive,MedianDrive_Compare);
    [p2]=ranksum(MedianDrive,MedianDrive_Compare)
    [nanmean(MedianDrives) nanstd(MedianDrives) nanmedian(MedianDrives) tsnaniqr(MedianDrives)]
    
    HighestDrive = [Vdrivedeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    HighestDrive_Compare = [Vdrivedeciles_compare(rowHigh,:)]';
    HighestDrives = [HighestDrive HighestDrive_Compare];
    [~,p]=ttest2(HighestDrive,HighestDrive_Compare);
    [p2]=ranksum(HighestDrive,HighestDrive_Compare)
    
    [nanmean(HighestDrives) nanstd(HighestDrives) nanmedian(HighestDrives) tsnaniqr(HighestDrives)]
    
    
    %data are in percent of eupneic reference value, from Vpassive in compare (e.g. NREM)
    LowestGGt = [GGtdeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    LowestGGt_Compare = [GGtdeciles_compare(rowLow,:)]';
    LowestGGts = [LowestGGt LowestGGt_Compare];
    [~,p]=ttest2(LowestGGt,LowestGGt_Compare);
    [p2]=ranksum(LowestGGt,LowestGGt_Compare)
    
    [nanmean(LowestGGts) nanstd(LowestGGts) nanmedian(LowestGGts) tsnaniqr(LowestGGts)]
    
    
    MedianGGt = [GGtdeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    MedianGGt_Compare = [GGtdeciles_compare(rowMid,:)]';
    MedianGGts = [MedianGGt MedianGGt_Compare];
    [~,p]=ttest2(MedianGGt,MedianGGt_Compare);
    [p2]=ranksum(MedianGGt,MedianGGt_Compare)
    
    [nanmean(MedianGGts) nanstd(MedianGGts) nanmedian(MedianGGts)]
    
    HighestGGt = [GGtdeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    HighestGGt_Compare = [GGtdeciles_compare(rowHigh,:)]';
    HighestGGts = [HighestGGt HighestGGt_Compare];
    [~,p]=ttest2(HighestGGt,HighestGGt_Compare);
    [p2]=ranksum(HighestGGt,HighestGGt_Compare)
    
    [nanmean(HighestGGts) nanstd(HighestGGts) nanmedian(HighestGGts)]
    
    
    LowestGGp = [GGpdeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    LowestGGp_Compare = [GGpdeciles_compare(rowLow,:)]';
    LowestGGps = [LowestGGp LowestGGp_Compare];
    [~,p]=ttest2(LowestGGp,LowestGGp_Compare);
    [p2]=ranksum(LowestGGp,LowestGGp_Compare)
    
    [nanmean(LowestGGps) nanstd(LowestGGps) nanmedian(LowestGGps) tsnaniqr(LowestGGps)]
    
    MedianGGp = [GGpdeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    MedianGGp_Compare = [GGpdeciles_compare(rowMid,:)]';
    MedianGGps = [MedianGGp MedianGGp_Compare];
    [~,p]=ttest2(MedianGGp,MedianGGp_Compare);
    [p2]=ranksum(MedianGGp,MedianGGp_Compare)
    
    [nanmean(MedianGGps) nanstd(MedianGGps) nanmedian(MedianGGps)]
    
    HighestGGp = [GGpdeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    HighestGGp_Compare = [GGpdeciles_compare(rowHigh,:)]';
    HighestGGps = [HighestGGp HighestGGp_Compare];
    [~,p]=ttest2(HighestGGp,HighestGGp_Compare);
    [p2]=ranksum(HighestGGp,HighestGGp_Compare)
    
    [nanmean(HighestGGps) nanstd(HighestGGps) nanmedian(HighestGGps)]
    
    
    
    GGtResp = (HighestGGt-LowestGGt)./(HighestDrive-LowestDrive);
    GGtResp_Compare = (HighestGGt_Compare-LowestGGt_Compare)./(HighestDrive_Compare-LowestDrive_Compare);
    GGtResps = [GGtResp GGtResp_Compare];
    [~,p]=ttest2(GGtResp,GGtResp_Compare);
    [p2]=ranksum(GGtResp,GGtResp_Compare)
    [nanmean(GGtResps) nanstd(GGtResps) nanmedian(GGtResps)]
    
    
    
    GGpResp = (HighestGGp-LowestGGp)./(HighestDrive-LowestDrive);
    GGpResp_Compare = (HighestGGp_Compare-LowestGGp_Compare)./(HighestDrive_Compare-LowestDrive_Compare);
    GGpResps = [GGpResp GGpResp_Compare];
    [~,p]=ttest2(GGpResp,GGpResp_Compare);
    [p2]=ranksum(GGpResp,GGpResp_Compare)
    [nanmean(GGpResps) nanstd(GGpResps) nanmedian(GGpResps)]
    
    
    
    LowestVE = [VEdeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    LowestVE_Compare = [VEdeciles_compare(rowLow,:)]';
    LowestVEs = [LowestVE LowestVE_Compare];
    [~,p]=ttest2(LowestVE,LowestVE_Compare);
    [p2]=ranksum(LowestVE,LowestVE_Compare);
    
    [nanmean(LowestVEs) nanstd(LowestVEs) nanmedian(LowestVEs) tsnaniqr(LowestVEs)]
    
    MedianVE = [VEdeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    MedianVE_Compare = [VEdeciles_compare(rowMid,:)]';
    MedianVEs = [MedianVE MedianVE_Compare];
    [~,p]=ttest2(MedianVE,MedianVE_Compare);
    [p2]=ranksum(MedianVE,MedianVE_Compare);
    
    [nanmean(MedianVEs) nanstd(MedianVEs) nanmedian(MedianVEs)]
    
    HighestVE = [VEdeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    HighestVE_Compare = [VEdeciles_compare(rowHigh,:)]';
    HighestVEs = [HighestVE HighestVE_Compare];
    [~,p]=ttest2(HighestVE,HighestVE_Compare);
    [p2]=ranksum(HighestVE,HighestVE_Compare);
    
    [nanmean(HighestVEs) nanstd(HighestVEs) nanmedian(HighestVEs) tsnaniqr(HighestVEs)]
    
    
    VEResp = (HighestVE-LowestVE)./(HighestDrive-LowestDrive);
    VEResp_Compare = (HighestVE_Compare-LowestVE_Compare)./(HighestDrive_Compare-LowestDrive_Compare);
    VEResps = [VEResp VEResp_Compare];
    [~,p]=ttest2(VEResp,VEResp_Compare);
    [p2]=ranksum(VEResp,VEResp_Compare);
    
    [nanmean(VEResps) nanstd(VEResps) nanmedian(VEResps)];
    
    
    
    
    LowestVE = [VEdeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    LowestVE_Compare = [VEdeciles_compare(rowLow,:)]';
    LowestVEs = [LowestVE LowestVE_Compare];
    [~,p]=ttest2(LowestVE,LowestVE_Compare);
    [p2]=ranksum(LowestVE,LowestVE_Compare);
    
    %
    % ArThresX = DataArray(:,5);
    % ArThresY = DataArray_compare(:,5);
    % Vdrivedeciles
    % Colstorun = find(~isnan(Vdrivedeciles(1,:)));
    %
    % rows = ones(1,size(Vdrivedeciles,2));
    % for i=1:length(Colstorun)
    %     rows(Colstorun(i)) = find(Vdrivedeciles(:,Colstorun(i))>=ArThresX(Colstorun(i)),1);
    % end
    
    %GGdataInclNorm = GGdataIncl./GGdataIncl(:,3)*100;
    %GGdataIncl_compareNorm = GGdataIncl_compare./GGdataIncl_compare(:,3)*100;
    
    GGpassivep = GGdataIncl(:,1);
    GGpassivep_compare = GGdataIncl_compare(:,1);
    
    %GGpassivep = [GGdataInclNorm(:,1)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    %GGpassivep_compare = [GGdataIncl_compareNorm(:,1)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGpassiveps = [GGpassivep GGpassivep_compare];
    [~,p]=ttest2(GGpassivep,GGpassivep_compare);
    [p2]=ranksum(GGpassivep,GGpassivep_compare);
    
    [nanmean(GGpassiveps) nanstd(GGpassiveps) nanmedian(GGpassiveps) tsnaniqr(GGpassiveps)]
    
    
    GGpassivet = [GGdataIncl(:,3)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGpassivet_compare = [GGdataIncl_compare(:,3)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGpassivets = [GGpassivet GGpassivet_compare];
    [~,p]=ttest2(GGpassivet,GGpassivet_compare);
    [p2]=ranksum(GGpassivet,GGpassivet_compare);
    
    [nanmean(GGpassivets) nanstd(GGpassivets) nanmedian(GGpassivets) tsnaniqr(GGpassivets)]
    
    
    GGactivep = [GGdataIncl(:,2)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGactivep_compare = [GGdataIncl_compare(:,2)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGactiveps = [GGactivep GGactivep_compare];
    [~,p]=ttest2(GGactivep,GGactivep_compare);
    [p2]=ranksum(GGactivep,GGactivep_compare)
    
    [nanmean(GGactiveps) nanstd(GGactiveps) nanmedian(GGactiveps) tsnaniqr(GGactiveps)]
    
    
    GGactivet = [GGdataIncl(:,4)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGactivet_compare = [GGdataIncl_compare(:,4)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
    GGactivets = [GGactivet GGactivet_compare];
    [~,p]=ttest2(GGactivet,GGactivet_compare);
    [p2]=ranksum(GGactivet,GGactivet_compare);
    
    [nanmean(GGactivets) nanstd(GGactivets) nanmedian(GGactivets) tsnaniqr(GGactivets)]
    
    Comp = GGactivet - LowestGGt;
    Comp_compare = GGactivet_compare - LowestGGt_Compare;
    [p2]=ranksum(Comp,Comp_compare)
    Comps = [Comp Comp_compare]
    [nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]
    
    Comp = GGactivep - LowestGGp;
    Comp_compare = GGactivep_compare - LowestGGp_Compare;
    [p2]=ranksum(Comp,Comp_compare)
    Comps = [Comp Comp_compare]
    [nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]
    
    arthres = DataArrayIncl(:,5);
    arthres_compare = DataArrayIncl_compare(:,5);
    [p2]=ranksum(arthres,arthres_compare)
    arthress = [arthres arthres_compare];
    [nanmean(arthress) nanstd(arthress) nanmedian(arthress) tsnaniqr(arthress)]
    
    %Vactive
    arthres = DataArrayIncl(:,7);
    arthres_compare = DataArrayIncl_compare(:,7);
    [p2]=ranksum(arthres,arthres_compare)
    arthress = [arthres arthres_compare];
    [nanmean(arthress) nanstd(arthress) nanmedian(arthress) tsnaniqr(arthress)]
    
    % %Vactive
    % VRA = VRAmean;
    % VRA_compare = VRAmean;
    % VRA(Iexclude==1)=NaN;
    % VRA_compare(Iexclude_compare==1)=NaN;
    % [p2]=ranksum(VRA,VRA_compare)
    % VRAs = [VRA VRA_compare];
    % [nanmean(VRAs) nanstd(VRAs) nanmedian(VRAs) tsnaniqr(VRAs)]
    
    
    
    Comp = (GGactivet - LowestGGt)./(DataArray(:,5) - LowestDrive);
    Comp_compare = (GGactivet_compare - LowestGGt_Compare)./(DataArray_compare(:,5) - LowestDrive_Compare);
    
    Comp = (GGactivep - LowestGGp)./(DataArray(:,5) - LowestDrive)
    Comp_compare = (GGactivep_compare - LowestGGp_Compare)./(DataArray_compare(:,5) - LowestDrive_Compare)
    
    
    nanmedian([GGpassivep;GGpassivep_compare])
    
    [p2]=ranksum(Comp,Comp_compare)
    [p2]=ranksum(Comp,Comp_compare)
    [~,p]=ttest2(Comp,Comp_compare)
    [~,p]=ttest2(log10(Comp),log10(Comp_compare))
    Comps = [Comp Comp_compare];
    [nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]    
    figure(99); hist(Comps(:),60)
    
    Comp = (DataArray(:,7) - LowestVE)./(DataArray(:,5) - LowestDrive);
    Comp_compare = (DataArray_compare(:,7) - LowestVE_Compare)./(DataArray_compare(:,5) - LowestDrive_Compare);
    
    
    Comp = (GGactivep - GGpassivep)./(DataArray(:,5) - 100)
    Comp_compare = (GGactivep_compare - GGpassivep_compare)./(DataArray_compare(:,5) - 100)
    
    [p2]=ranksum(Comp,Comp_compare)
    Comps = [Comp Comp_compare];
    [nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]
    
    
    MuscleEffectiveness = (GGactivep - GGpassivep)./(arthresMedian - 100)
    MuscleEffectiveness = (VactiveMedian - VpassiveMedian)./(arthresMedian  - 100)
    MuscleEffectiveness_compare = (VactiveMedian_compare  - VpassiveMedian_compare )./(arthresMedian_compare - 100)
    
    
    Comp_compare = (GGactivep_compare - GGpassivep_compare)./(DataArray_compare(:,5) - 100)
    
end

%%
    %MuscleEffectiveness = (GGactivep - GGpassivep)./(arthresMedian - 100)
    
    MuscleEffectiveness = (DataArray(:,7) - DataArray(:,6))./(DataArray(:,5)  - 100) *100;
    if 0
    ME=nan(size(DataArray,1),2);
    for i=1:size(DataArray,1)
        figure(101); clf(101);
        Ydata = VEdeciles(:,i);
        Xdata = Vdrivedeciles(:,i);
        if sum(isnan(Ydata)|isnan(Xdata))>0
            continue
        end
        plot(Xdata,Ydata)
        [p,S] = polyfit(Xdata-100,Ydata,1);
        Xline = [min(Xdata-100) max(Xdata-100)];
        Yline = polyval(p,Xline);
        hold on
        plot(Xline+100,Yline);
        ME(i,:)=[p(1)*100 p(2)];
        %pause(0.1)
    end
    
    
    GGEffectiveness = (GGdata(:,2) - GGdata(:,1))./GGdata(:,1)*100./(DataArray(:,5)  - 100) *100;
    
    GGpdecilesAll = GGpdeciles;
    X = ~isnan(GGpdeciles_compare(1,:));
    GGpdecilesAll(:,X) = GGpdeciles_compare(:,X);
    
    VdrivedecilesAll = Vdrivedeciles;
    X = ~isnan(Vdrivedeciles_compare(1,:));
    VdrivedecilesAll(:,X) = Vdrivedeciles_compare(:,X);
      
    
    
    GGE=nan(size(GGdata,1),2);
    for i=1:size(GGdata,1)
        figure(101); clf(101);
        %Ydata = GGpdecilesAll(:,i)/GGdata(i,1)*100;
        Ydata = GGpdecilesAll(:,i);
        Xdata = VdrivedecilesAll(:,i);
        if sum(isnan(Ydata)|isnan(Xdata))>0
            continue
        end
        plot(Xdata,Ydata)
        [p,S] = polyfit(Xdata-100,Ydata,1);
        Xline = [min(Xdata-100) max(Xdata-100)];
        Yline = polyval(p,Xline);
        hold on
        plot(Xline+100,Yline);
        GGE(i,:)=[p(1) p(2)];
        pause(0.1)
    end
    end
    GGEffectiveness2 = (GGdata(:,2) - GGdata(:,1))./(DataArray(:,5)  - 100) *100;
    
    %MuscleEffectiveness_compare = (VactiveMedian_compare  - VpassiveMedian_compare )./(arthresMedian_compare - 100)
    
    %%
%     figure(103)
%     scatter(ME(:,1),MuscleEffectiveness,15,[1 0 0])
%     
%      figure(105)
%     scatter(GGE(:,1),GGEffectiveness,15,[0 1 0])

% %%
% displaytext='Saving Summary Analysis';
% disp(displaytext); set(handletext,'String',displaytext); drawnow;
% savefilename = [settings.directoryout 'SummaryAnalysis_' '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '_' datestr(now,'yyyymmdd HHMMSS')];
% savefilename(end-6)='T';
% save(savefilename)
% savefilename = [settings.directoryout 'EndoplotStats_AvB'];
% save(savefilename)

GGpassiveAll = GGdata(:,1);

%%
addpath(genpath('E:\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta\'))
try
    h=figure(101)
    %Data = MuscleEffectiveness;
    %Data = fArThresT(DataArray(:,5));
    %Data = (DataArray(:,6));
    %Data = (DataArray(:,1));
    %Data = GGEffectiveness2;
    %swtest(Data)
    
    %Data = ME(:,1)
    Data = GGE(:,1).^.5
    %Data = GGpassiveAll(:,1).^.5
    crit2 = P==1;
    crit = P==0|P==-1|P==-2;
    [~,pval,ci]=ttest2(Data(crit2),Data(crit));
    [pval]=ranksum(Data(crit2),Data(crit))
    meanandsemdiff = [nanmean(Data(crit2))-nanmean(Data(crit)),ci',pval]
    PlotThreeGroupIndividualData(h,Data,crit,crit2)
    pranksum = ranksum(Data(crit),Data(crit2));
    %temp = fArThresBT(Data)
end

%%
addpath(genpath('E:\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta\'))
try
    h=figure(101)
    Data = MuscleEffectiveness;
    Data = fArThresT(DataArray(:,5));
   % Data = (DataArray(:,8));
    Data = (DataArray(:,8));
    %Data = (DataArray(:,1));
    %Data = GGEffectiveness;
    %Data = GGEffectiveness2;
    %swtest(Data)
    
    %Data = ME(:,1)
    %Data = GGE(:,1)
    %Data = GGE
    crit2 = P==1;
    crit = P==-1 |P==0 |P==-2; %-1
    crit3=~isnan(P)
    [~,pval,ci]=ttest2(Data(crit2),Data(crit))
    
    [nanmean(Data(crit3)) nanstd(Data(crit3)) nanmedian(Data(crit3)) tsnaniqr(Data(crit3))]
    [nanmean(Data(crit)) nanstd(Data(crit)) nanmedian(Data(crit)) tsnaniqr(Data(crit))]
    [nanmean(Data(crit2)) nanstd(Data(crit2)) nanmedian(Data(crit2)) tsnaniqr(Data(crit2))]
    
    meanandsemdiff = [nanmean(Data(crit2))-nanmean(Data(crit)),ci',pval]
    PlotThreeGroupIndividualData(h,Data,crit,crit2)
    pranksum = ranksum(Data(crit),Data(crit2))
    %temp = fArThresBT(Data)
end


% P(raw=="A")=1;
% P(raw=="B1")=0; %nonresponsive [ME down, Comp down]
% P(raw=="B2")=-1; %inefficient NM [collapsibility down]
% P(raw=="B")=-2;

%%

%RisingFalling stats
% 
% X(:,1)= VdrivedecilesMedian_compare;
% X(:,2)= VEdecilesMedian_compare;
% X(:,3)= VdrivedecilesMedian;
% X(:,4)= VEdecilesMedian;




