%% Load Summary Analysis: s5p2
clear all
%%
%choose path
path1 = 'C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\'
path1 = 'G:\Partners Healthcare Dropbox\SATP Group\'
path1 = 'C:\Users\ludov\Dropbox\'

%load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_PhasicREMvNREM_20201110T003801.mat']) % PhasicREM
%load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_20201110T001911.mat']) % REM
%load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_TonicREMvNREM_20201111T164124.mat']) % Tonic REM

%retest R1
load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_20201110T001911.mat']) 
%load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_R1retest_20210429T094526.mat'])% REM retest
%load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_R1supine_20210429T101004.mat']) % REM supine
%load([path1 'PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_R1lateral_20210505T235028.mat']) % REM lateral
%edit these as needed:

%load('G:\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\REMvNREM\SummaryAnalysis_TonicREMvNREM_20201111T164124.mat')
%load('G:\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_Pes_20201111T171656.mat')
%load('G:\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\REMvNREM\SummaryAnalysis_REMvNREM_Pmus_20201111T173356.mat')
%
load([path1 'PhenotypeDrive2018\REMvNREM\getData_20210119T105333.mat'])

addpath(genpath('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PUPbeta_git'))
addpath(genpath('C:\Users\ludov\Dropbox\PUPbeta_git'))
%%

Isubject = [2 10 11 16 19 20 25 26 27 30 34 37 40 43 44 46 48 49 50 51 53 55 58 59 61];
temp=T;
temp{:,:}=NaN;
temp(Isubject,:)=T(Isubject,:);
T = temp;

%%
rowLow = 1; %contains values 0-20th centile, centered on 10th centile
rowMid = 5;
rowHigh = 10; %16,66,116
%note 1st value contains 0-5th centile, centered on 2.5th centile.

usemedianCI=0;

temp = nanmean(Vdrivedeciles_compare([5 6],:))';
temp = nanmean(Vdrivedeciles([5 6],:))';
I2 = find(~isnan(temp))' %have this data
%41
Nsubjects = length(I2);
I=[]; %additional test exclusions (off)
%19 26 55 missing phasic?

%% Get the AHI data from tables


AHInrem = AHIdata(:,88); 
AHIrem = AHIdata(:,96);
AHItotal = AHIdata(:,80);

ArInrem = AHIdata(:,82); 
ArIrem = AHIdata(:,90);
ArItotal = AHIdata(:,74);


TSTnrem = AHIdata(:,81); 
TSTrem = AHIdata(:,89);
TSTtotal = AHIdata(:,73);

[nanmean(TSTtotal)  nanstd(TSTtotal)]
[nanmean(TSTnrem)  nanstd(TSTnrem)]
[nanmean(TSTrem)  nanstd(TSTrem)]

TSTnremS = AHIdata(:,9); 
TSTremS = AHIdata(:,17);
TSTtotalS = AHIdata(:,1);

[nanmean(TSTtotalS./TSTtotal)  nanstd(TSTtotalS./TSTtotal)]*100
[nanmean(TSTnremS./TSTnrem)  nanstd(TSTnremS./TSTnrem)]*100
[nanmean(TSTremS./TSTrem)  nanstd(TSTremS./TSTrem)]*100

OAInrem = AHIdata(:,83); 
OAIrem = AHIdata(:,91);
OAItotal = AHIdata(:,75);
CAInrem = AHIdata(:,84); 
CAIrem = AHIdata(:,92);
CAItotal = AHIdata(:,76);
OHInrem = AHIdata(:,85); 
OHIrem = AHIdata(:,93);
OHItotal = AHIdata(:,77);
MAInrem = AHIdata(:,86); 
MAIrem = AHIdata(:,94);
MAItotal = AHIdata(:,78);
CHInrem = AHIdata(:,87); 
CHIrem = AHIdata(:,95);
CHItotal = AHIdata(:,79);

HBnrem = T.HBnrem; 
HBrem = T.HBrem;
HBtotal = T.HBtotal;

EDnrem = T.EventDurationMeannrem; 
EDrem = T.EventDurationMeanrem;
EDtotal = T.EventDurationMean;



AHInremT = AHInrem.^0.5;
AHIremT = AHIrem.^0.5;
AHItotalT = AHItotal.^0.5;
ArInremT = ArInrem.^0.5;
ArIremT = ArIrem.^0.5;
ArItotalT = ArItotal.^0.5;
OAInremT = OAInrem.^0.5;
OAIremT = OAIrem.^0.5;
OAItotalT = OAItotal.^0.5;
CAInremT = CAInrem.^0.5;
CAIremT = CAIrem.^0.5;
CAItotalT = CAItotal.^0.5;
OHInremT = OHInrem.^0.5;
OHIremT = OHIrem.^0.5;
OHItotalT = OHItotal.^0.5;
MAInremT = MAInrem.^0.5;
MAIremT = MAIrem.^0.5;
MAItotalT = MAItotal.^0.5;
CHInremT = CHInrem.^0.5;
CHIremT = CHIrem.^0.5;
CHItotalT = CHItotal.^0.5;


HBnremT = HBnrem.^0.5;
HBremT = HBrem.^0.5;
HBtotalT = HBtotal.^0.5;

% %minAHIforREMnREMbalance=5;
% I = find(round(AHItotal)<=5 & round(AHIrem)<=5) %total AHI >5, or AHIrem >5
% %I = find((AHInrem+AHIrem)/2 < minAHIforREMnREMbalance)
% AHInrem(I)=NaN;
% AHIrem(I)=NaN;
% %AHIratioT = 100*(AHIrem-AHInrem)./(AHIrem+AHInrem);
% %nanmean(AHIratioT)
% %[~,p]=ttest(AHIratioT,zeros(61,1))

Y = HBnremT;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT.^2;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

Y = HBremT;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT.^2;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

Y = HBtotalT;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT.^2;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

[~,p]=ttest(HBnremT,HBremT)
swtest(HBrem.^0.5)


Y = ArItotalT;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT.^2;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

Y = ArInremT;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT.^2;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

Y = ArIremT;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT.^2;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

[~,p]=ttest(ArInremT,ArIremT)

swtest(AHIrem.^0.5)


Y = EDtotal;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

Y = EDnrem;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

Y = EDrem;
yT = nanmean(Y) + nanstd(Y).*[-1 0 +1];
yUT = yT;
ymeanandsdUT = [yUT(2) nanmean(diff(yUT))]

[~,p]=ttest(EDnrem,EDrem)

swtest(EDtotal)


% find(round((AHInrem+AHIrem)/2)<5) %average of AHIrem and AHInrem >5
% find(((AHInrem+AHIrem)/2)<5) %average of AHIrem and AHInrem >5
% find(((AHInrem+AHIrem)/2)<=5)
% find(round(AHItotal)<=5 & round(AHIrem)<=5) %total AHI >5, or AHIrem >5





%% FlowVsDrive PhenoPlot Average
if settings.runcompare
    Iexclude_compare = DataNArray_compare(:,6)<2 | DataNArray(:,6)<2;
    Iexclude = DataNArray_compare(:,6)<2 | DataNArray(:,6)<2;
    Iexclude_compare(I)=1;
    Iexclude(I)=1;
    VEdeciles_compare(:,Iexclude_compare==1)=NaN;
    Vdrivedeciles_compare(:,Iexclude_compare==1)=NaN;
    DataArrayIncl_compare=DataArray_compare;
    DataArrayIncl_compare(Iexclude_compare==1,:)=NaN;
    GGpdeciles_compare(:,Iexclude_compare==1)=NaN;
    GGtdeciles_compare(:,Iexclude_compare==1)=NaN;
    
    if nanmedian(DataArray_compare(:,9))<10
        DataArray_compare(:,9) = DataArray_compare(:,9)*100;
    end
    
    if nanmedian(DataArray(:,9))<10
        DataArray(:,9) = DataArray(:,9)*100;
    end
else
    Iexclude = DataNArray(:,6)<2;
    Iexclude(I)=1;
    
    if nanmedian(DataArray(:,9))<10
        DataArray(:,9) = DataArray(:,9)*100;
    end
end
VEdeciles(:,Iexclude==1)=NaN;
Vdrivedeciles(:,Iexclude==1)=NaN;
DataArrayIncl=DataArray;
DataArrayIncl(Iexclude==1,:)=NaN;
GGpdeciles(:,Iexclude==1)=NaN;
GGtdeciles(:,Iexclude==1)=NaN;

GGdataIncl = GGdata;
GGdataIncl(Iexclude==1,:)=NaN;

if settings.runcompare
    GGdataIncl_compare = GGdata_compare;
    GGdataIncl_compare(Iexclude==1,:)=NaN;
end

if settings.runcompare %normalize by VpassiveCompare
    GGpdeciles = 100*GGpdeciles./GGdataIncl_compare(:,1)';
    GGtdeciles = 100*GGtdeciles./GGdataIncl_compare(:,1)';
    GGpdeciles_compare = 100*GGpdeciles_compare./GGdataIncl_compare(:,1)';
    GGtdeciles_compare = 100*GGtdeciles_compare./GGdataIncl_compare(:,1)';
end
% if positioncompare
%     Iexclude_compare1 = DataNArray_compare1(:,6)<2;
%     VEdeciles_compare1(:,Iexclude_compare1==1)=NaN;
%     Vdrivedeciles_compare1(:,Iexclude_compare1==1)=NaN;
%     DataArrayIncl_compare1=DataArray_compare1;
%     DataArrayIncl_compare1(Iexclude_compare1==1,:)=NaN;
% end

%% GGp
figure(10); clf(10);

%prep data
GGpdecilesUpper_compare = prctile(GGpdeciles_compare',75)';
GGpdecilesLower_compare = prctile(GGpdeciles_compare',25)';
VdrivedecilesMedian_compare = prctile(Vdrivedeciles_compare',50)';
GGpdecilesMedian_compare = prctile(GGpdeciles_compare',50)';
GGpdecilesUpper = prctile(GGpdeciles',75)';
GGpdecilesLower = prctile(GGpdeciles',25)';
VdrivedecilesMedian = prctile(Vdrivedeciles',50)';
GGpdecilesMedian = prctile(GGpdeciles',50)';
if 1
    GGpdecilesUpper_compare = filter121(GGpdecilesUpper_compare);
    GGpdecilesLower_compare = filter121(GGpdecilesLower_compare);
    GGpdecilesUpper = filter121(GGpdecilesUpper);
    GGpdecilesLower = filter121(GGpdecilesLower);
end
if usemedianCI
    [~,temp] = medianCI(GGpdeciles'); GGpdecilesUpper=filter121(temp(2,:)'); GGpdecilesLower=filter121(temp(1,:)');
    [~,temp] = medianCI(GGpdeciles_compare'); GGpdecilesUpper_compare=filter121(temp(2,:)'); GGpdecilesLower_compare=filter121(temp(1,:)');
end

GGtdecilesUpper_compare = prctile(GGtdeciles_compare',75)';
GGtdecilesLower_compare = prctile(GGtdeciles_compare',25)';
GGtdecilesMedian_compare = prctile(GGtdeciles_compare',50)';
GGtdecilesUpper = prctile(GGtdeciles',75)';
GGtdecilesLower = prctile(GGtdeciles',25)';
GGtdecilesMedian = prctile(GGtdeciles',50)';
if 1
    GGtdecilesUpper_compare = filter121(GGtdecilesUpper_compare);
    GGtdecilesLower_compare = filter121(GGtdecilesLower_compare);
    GGtdecilesUpper = filter121(GGtdecilesUpper);
    GGtdecilesLower = filter121(GGtdecilesLower);
end
if usemedianCI
    [~,temp] = medianCI(GGtdeciles_compare'); GGtdecilesUpper_compare=temp(2,:)'; GGtdecilesLower_compare=temp(1,:)';
    [~,temp] = medianCI(GGtdeciles'); GGtdecilesUpper=temp(2,:)'; GGtdecilesLower=temp(1,:)';
end

arthresMedian = nanmedian(DataArrayIncl(:,5));
arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));

[~,ii] = unique(VdrivedecilesMedian_compare);
GGppassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),GGpdecilesMedian_compare(ii),100,'linear');
GGpactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),GGpdecilesMedian_compare(ii),arthresMedian_compare,'linear');

[~,ii] = unique(VdrivedecilesMedian);
GGppassiveMedian = interp1(VdrivedecilesMedian(ii),GGpdecilesMedian(ii),100,'linear');
GGpactiveMedian = interp1(VdrivedecilesMedian(ii),GGpdecilesMedian(ii),arthresMedian,'linear');

set(gcf,'color',[1 1 1]);
if 1
    plot([100 100],[0 110],'--','color',0.5*[1 1 1]);
    hold('on')
    plot([0 300],[100 100],'--','color',0.5*[1 1 1]);
    box('off');
    set(gca,'tickdir','out')
end

plot([arthresMedian_compare arthresMedian_compare],[0 GGpactiveMedian_compare],'g'); %Arthres line
plot([arthresMedian arthresMedian],[0 GGpactiveMedian],'g'); %Arthres line

fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[GGpdecilesUpper_compare;flipud(GGpdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian_compare,GGpdecilesMedian_compare,'k','linewidth',2);

fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGpdecilesUpper;flipud(GGpdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,GGpdecilesMedian,'k','linewidth',2);

DispEventLikelihoodOnEndogram(VdrivedecilesMedian_compare,GGpdecilesMedian_compare,nanmean(1-EdecilesMean_compare')',0);
DispEventLikelihoodOnEndogram(VdrivedecilesMedian,GGpdecilesMedian,nanmean(1-EdecilesMean')',0)

if settings.runcompare
    ylim([0 max(max(GGpdecilesUpper),max(GGpdecilesUpper_compare))])
else
    ylim([0 (max(GGpdecilesUpper))])
end

plot(100,GGppassiveMedian,'.','markersize',10,'color',[0 0 0]); %Vpassive
plot(100,GGppassiveMedian_compare,'.','markersize',10,'color',[0 0 0]); %Vpassive

%   GGt

[~,ii] = unique(VdrivedecilesMedian_compare);
GGtpassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),GGtdecilesMedian_compare(ii),100,'linear');
GGtactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),GGtdecilesMedian_compare(ii),arthresMedian_compare,'linear');

[~,ii] = unique(VdrivedecilesMedian);
GGtpassiveMedian=interp1(VdrivedecilesMedian(ii),GGtdecilesMedian(ii),100,'linear');
GGtactiveMedian = interp1(VdrivedecilesMedian(ii),GGtdecilesMedian(ii),arthresMedian,'linear');

fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[GGtdecilesUpper_compare;flipud(GGtdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian_compare,GGtdecilesMedian_compare,'k','linewidth',2);

fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGtdecilesUpper;flipud(GGtdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,GGtdecilesMedian,'k','linewidth',2);



DispEventLikelihoodOnEndogram(VdrivedecilesMedian,GGtdecilesMedian,nanmean(1-EdecilesMean')',0)
DispEventLikelihoodOnEndogram(VdrivedecilesMedian_compare,GGtdecilesMedian_compare,nanmean(1-EdecilesMean_compare')',0)

plot(100,GGtpassiveMedian,'.','markersize',10,'color',[0 0 0]); %Vpassive
plot(100,GGtpassiveMedian_compare,'.','markersize',10,'color',[0 0 0]); %Vpassive

xlim([0 275])
ylim([0 max(max(GGpdecilesUpper),max(GGpdecilesUpper_compare))])

yticks([0:50:250]);

%%    Ventilation
figure(9); clf(9);

set(gcf,'color',[1 1 1]);
if 1
    plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
    hold('on')
    plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
end

    VEdecilesUpper_compare = prctile(VEdeciles_compare',75)';
    VEdecilesLower_compare = prctile(VEdeciles_compare',25)';
    VdrivedecilesMedian_compare = prctile(Vdrivedeciles_compare',50)';
    VEdecilesMedian_compare = prctile(VEdeciles_compare',50)';
    
    VEdecilesUpper = prctile(VEdeciles',75)';
    VEdecilesLower = prctile(VEdeciles',25)';
    VdrivedecilesMedian = prctile(Vdrivedeciles',50)';
    VEdecilesMedian = prctile(VEdeciles',50)';

    if 1
        VEdecilesUpper_compare = filter121(VEdecilesUpper_compare);
        VEdecilesLower_compare = filter121(VEdecilesLower_compare);
        VEdecilesUpper = filter121(VEdecilesUpper);
        VEdecilesLower = filter121(VEdecilesLower);
    end
    if usemedianCI
        [~,temp] = medianCI(VEdeciles_compare'); VEdecilesUpper_compare=temp(2,:)'; VEdecilesLower_compare=temp(1,:)';
        [~,temp] = medianCI(VEdeciles'); VEdecilesUpper=temp(2,:)'; VEdecilesLower=temp(1,:)';
    end
    
[~,ii] = unique(VdrivedecilesMedian_compare);
arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
VpassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),100,'linear');
VactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),arthresMedian_compare,'linear');

[~,ii] = unique(VdrivedecilesMedian);
arthresMedian = nanmedian(DataArrayIncl(:,5));
VpassiveMedian=interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),100,'linear');
VactiveMedian = interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),arthresMedian,'linear');

plot([arthresMedian_compare arthresMedian_compare],[0 VactiveMedian_compare],'g'); %Arthres line
plot([arthresMedian arthresMedian],[0 VactiveMedian],'g'); %Arthres line


   
    
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[VEdecilesUpper_compare;flipud(VEdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,VEdecilesMedian_compare,'k','linewidth',2);

fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[VEdecilesUpper;flipud(VEdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,VEdecilesMedian,'k','linewidth',2);

DispEventLikelihoodOnEndogram(VdrivedecilesMedian_compare,VEdecilesMedian_compare,nanmean(1-EdecilesMean_compare')',1)
DispEventLikelihoodOnEndogram(VdrivedecilesMedian,VEdecilesMedian,nanmean(1-EdecilesMean')',1)

plot(100,VpassiveMedian_compare,'.','markersize',10,'color',[0 0 0]); %Vpassive
plot(100,VpassiveMedian,'.','markersize',10,'color',[0 0 0]); %Vpassive

ylim([0 125])
xlim([0 275])
%%
% figure(10); clf(10);
% if positioncompare
%     VEdecilesUpper_compare1 = prctile(VEdeciles_compare1',75)';
%     VEdecilesLower_compare1 = prctile(VEdeciles_compare1',25)';
%     VdrivedecilesMedian_compare1 = prctile(Vdrivedeciles_compare1',50)';
%     VEdecilesMedian_compare1 = prctile(VEdeciles_compare1',50)';
%
%     set(gcf,'color',[1 1 1]);
%     if 1
%         plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
%         hold('on')
%         plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
%         box('off');
%         set(gca,'tickdir','out')
%     end
%     fill([VdrivedecilesMedian_compare1;flipud(VdrivedecilesMedian_compare1)],[VEdecilesUpper_compare1;flipud(VEdecilesLower_compare1)],[0.9500 0.2500 0.1],'edgecolor','none','facealpha',0.5);
%     plot(VdrivedecilesMedian_compare1,VEdecilesMedian_compare1,'k','linewidth',2);
%
%     %VpassiveMedian_compare1 = nanmedian(DataArrayIncl_compare1(:,6));
%
%     VpassiveMedian_compare1=interp1(VdrivedecilesMedian_compare1,VEdecilesMedian_compare1,100,'linear');
%
%     arthresMedian_compare1 = nanmedian(DataArrayIncl_compare1(:,5));
%     VactiveMedian_compare1 = nanmedian(DataArrayIncl_compare1(:,7));
%
%     VactiveMedian_compare1 = interp1(VdrivedecilesMedian_compare1,VEdecilesMedian_compare1,arthresMedian_compare1,'linear');
%     plot(100,VpassiveMedian_compare1,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
%     plot([arthresMedian_compare1 arthresMedian_compare1],[0 max(arthresMedian_compare1,VactiveMedian_compare1)],'g'); %Arthres line
%     plot(arthresMedian_compare1,VactiveMedian_compare1,'r.','markersize',25);  %Vactive
%
% end
%
% VEdecilesUpper = prctile(VEdeciles',75)';
% VEdecilesLower = prctile(VEdeciles',25)';
% VdrivedecilesMedian = prctile(Vdrivedeciles',50)';
% VEdecilesMedian = prctile(VEdeciles',50)';
%
%     if ~positioncompare
%     plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
%     hold('on')
%     plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
%     box('off');
%     set(gca,'tickdir','out')
%     end
%     fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[VEdecilesUpper;flipud(VEdecilesLower)],[0.100 0.8500 0.6],'edgecolor','none','facealpha',0.5);
% plot(VdrivedecilesMedian,VEdecilesMedian,'k','linewidth',2);
%
%
%     VpassiveMedian = nanmedian(DataArrayIncl(:,6));
%     arthresMedian = nanmedian(DataArrayIncl(:,5));
%     VactiveMedian = nanmedian(DataArrayIncl(:,7));
%
%     VpassiveMedian=interp1(VdrivedecilesMedian,VEdecilesMedian,100,'linear');
%     VactiveMedian = interp1(VdrivedecilesMedian,VEdecilesMedian,arthresMedian,'linear');
%
%
%
%     plot(100,VpassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
%     plot([arthresMedian arthresMedian],[0 max(arthresMedian,VactiveMedian)],'g'); %Arthres line
%     plot(arthresMedian,VactiveMedian,'r.','markersize',25);  %Vactive
%
%
%     saveas(10,[settings.directoryout 'Fig10_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png'])
%
%% Analysis
%Vdrivedeciles_comparex = Vdrivedeciles_compare; temp

LowestDrive = [Vdrivedeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
LowestDrive_Compare = [Vdrivedeciles_compare(rowLow,:)]';
LowestDrives = [LowestDrive LowestDrive_Compare];
[~,p]=ttest2(LowestDrive,LowestDrive_Compare);
[p2]=ranksum(LowestDrive,LowestDrive_Compare);

[nanmean(LowestDrives) nanstd(LowestDrives) nanmedian(LowestDrives) tsnaniqr(LowestDrives)]

MedianDrive = [Vdrivedeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
MedianDrive_Compare = [Vdrivedeciles_compare(rowMid,:)]';
MedianDrives = [MedianDrive MedianDrive_Compare];
[~,p]=ttest2(MedianDrive,MedianDrive_Compare);
[p2]=ranksum(MedianDrive,MedianDrive_Compare);
[nanmean(MedianDrives) nanstd(MedianDrives) nanmedian(MedianDrives) tsnaniqr(MedianDrives)]

HighestDrive = [Vdrivedeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
HighestDrive_Compare = [Vdrivedeciles_compare(rowHigh,:)]';
HighestDrives = [HighestDrive HighestDrive_Compare];
[~,p]=ttest2(HighestDrive,HighestDrive_Compare);
[p2]=ranksum(HighestDrive,HighestDrive_Compare);

[nanmean(HighestDrives) nanstd(HighestDrives) nanmedian(HighestDrives) tsnaniqr(HighestDrives)]


%data are in percent of eupneic reference value, from Vpassive in compare (e.g. NREM)
LowestGGt = [GGtdeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
LowestGGt_Compare = [GGtdeciles_compare(rowLow,:)]';
LowestGGts = [LowestGGt LowestGGt_Compare];
[~,p]=ttest2(LowestGGt,LowestGGt_Compare);
[p2]=ranksum(LowestGGt,LowestGGt_Compare);

[nanmean(LowestGGts) nanstd(LowestGGts) nanmedian(LowestGGts) tsnaniqr(LowestGGts)]


MedianGGt = [GGtdeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
MedianGGt_Compare = [GGtdeciles_compare(rowMid,:)]';
MedianGGts = [MedianGGt MedianGGt_Compare];
[~,p]=ttest2(MedianGGt,MedianGGt_Compare);
[p2]=ranksum(MedianGGt,MedianGGt_Compare);

[nanmean(MedianGGts) nanstd(MedianGGts) nanmedian(MedianGGts)]

HighestGGt = [GGtdeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
HighestGGt_Compare = [GGtdeciles_compare(rowHigh,:)]';
HighestGGts = [HighestGGt HighestGGt_Compare];
[~,p]=ttest2(HighestGGt,HighestGGt_Compare);
[p2]=ranksum(HighestGGt,HighestGGt_Compare);

[nanmean(HighestGGts) nanstd(HighestGGts) nanmedian(HighestGGts)]


LowestGGp = [GGpdeciles(rowLow,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
LowestGGp_Compare = [GGpdeciles_compare(rowLow,:)]';
LowestGGps = [LowestGGp LowestGGp_Compare];
[~,p]=ttest2(LowestGGp,LowestGGp_Compare);
[p2]=ranksum(LowestGGp,LowestGGp_Compare);

[nanmean(LowestGGps) nanstd(LowestGGps) nanmedian(LowestGGps) tsnaniqr(LowestGGps)]

MedianGGp = [GGpdeciles(rowMid,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
MedianGGp_Compare = [GGpdeciles_compare(rowMid,:)]';
MedianGGps = [MedianGGp MedianGGp_Compare];
[~,p]=ttest2(MedianGGp,MedianGGp_Compare);
[p2]=ranksum(MedianGGp,MedianGGp_Compare);

[nanmean(MedianGGps) nanstd(MedianGGps) nanmedian(MedianGGps)]

HighestGGp = [GGpdeciles(rowHigh,:)]'; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
HighestGGp_Compare = [GGpdeciles_compare(rowHigh,:)]';
HighestGGps = [HighestGGp HighestGGp_Compare];
[~,p]=ttest2(HighestGGp,HighestGGp_Compare);
[p2]=ranksum(HighestGGp,HighestGGp_Compare);

[nanmean(HighestGGps) nanstd(HighestGGps) nanmedian(HighestGGps)]



GGtResp = (HighestGGt-LowestGGt)./(HighestDrive-LowestDrive);
GGtResp_Compare = (HighestGGt_Compare-LowestGGt_Compare)./(HighestDrive_Compare-LowestDrive_Compare);
GGtResps = [GGtResp GGtResp_Compare];
[~,p]=ttest2(GGtResp,GGtResp_Compare);
[p2]=ranksum(GGtResp,GGtResp_Compare);

[nanmean(GGtResps) nanstd(GGtResps) nanmedian(GGtResps)]



GGpResp = (HighestGGp-LowestGGp)./(HighestDrive-LowestDrive);
GGpResp_Compare = (HighestGGp_Compare-LowestGGp_Compare)./(HighestDrive_Compare-LowestDrive_Compare);
GGpResps = [GGpResp GGpResp_Compare];
[~,p]=ttest2(GGpResp,GGpResp_Compare);
[p2]=ranksum(GGpResp,GGpResp_Compare);

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
GGdataInclNorm = GGdataIncl./GGdataIncl_compare(:,1)*100;

GGdataIncl_compareNorm = GGdataIncl_compare./GGdataIncl_compare(:,1)*100;


GGpassivep = [GGdataInclNorm(:,1)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGpassivep_compare = [GGdataIncl_compareNorm(:,1)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGpassiveps = [GGpassivep GGpassivep_compare];
[~,p]=ttest2(GGpassivep,GGpassivep_compare);
[p2]=ranksum(GGpassivep,GGpassivep_compare);

[nanmean(GGpassiveps) nanstd(GGpassiveps) nanmedian(GGpassiveps) tsnaniqr(GGpassiveps)]


GGpassivet = [GGdataInclNorm(:,3)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGpassivet_compare = [GGdataIncl_compareNorm(:,3)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGpassivets = [GGpassivet GGpassivet_compare];
[~,p]=ttest2(GGpassivet,GGpassivet_compare);
[p2]=ranksum(GGpassivet,GGpassivet_compare);

[nanmean(GGpassivets) nanstd(GGpassivets) nanmedian(GGpassivets) tsnaniqr(GGpassivets)]


GGactivep = [GGdataInclNorm(:,2)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGactivep_compare = [GGdataIncl_compareNorm(:,2)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGactiveps = [GGactivep GGactivep_compare];
[~,p]=ttest2(GGactivep,GGactivep_compare);
[p2]=ranksum(GGactivep,GGactivep_compare);

[nanmean(GGactiveps) nanstd(GGactiveps) nanmedian(GGactiveps) tsnaniqr(GGactiveps)]


GGactivet = [GGdataInclNorm(:,4)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGactivet_compare = [GGdataIncl_compareNorm(:,4)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGactivets = [GGactivet GGactivet_compare];
[~,p]=ttest2(GGactivet,GGactivet_compare);
[p2]=ranksum(GGactivet,GGactivet_compare);

[nanmean(GGactivets) nanstd(GGactivets) nanmedian(GGactivets) tsnaniqr(GGactivets)]

Comp = GGactivet - LowestGGt;
Comp_compare = GGactivet_compare - LowestGGt_Compare;
[p2]=ranksum(Comp,Comp_compare)
Comps = [Comp Comp_compare];
[nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]

Comp = GGactivep - LowestGGp;
Comp_compare = GGactivep_compare - LowestGGp_Compare;
[p2]=ranksum(Comp,Comp_compare)
Comps = [Comp Comp_compare];
[nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]

arthres = DataArray(:,5);
arthres_compare = DataArray_compare(:,5);
[p2]=ranksum(arthres,arthres_compare)
arthress = [arthres arthres_compare];
[nanmean(arthress) nanstd(arthress) nanmedian(arthress) tsnaniqr(arthress)]

%Vactive
arthres = DataArray(:,6);
arthres_compare = DataArray_compare(:,6);
[p2]=ranksum(arthres,arthres_compare)
arthress = [arthres arthres_compare];
[nanmean(arthress) nanstd(arthress) nanmedian(arthress) tsnaniqr(arthress)]

if 1
    tempArthres = DataArray(:,5);
    tempArthres_compare = DataArray_compare(:,5);
else
    tempArthres = DataArray(:,9);
    tempArthres_compare = DataArray_compare(:,9);
end
CompGGt = 100*(GGactivet - LowestGGt)./(tempArthres - LowestDrive);
Comp_compareGGt = 100*(GGactivet_compare - LowestGGt_Compare)./(tempArthres_compare - LowestDrive_Compare);

CompGGp = 100*(GGactivep - LowestGGp)./(tempArthres - LowestDrive);
Comp_compareGGp = 100*(GGactivep_compare - LowestGGp_Compare)./(tempArthres_compare - LowestDrive_Compare);

CompVE = 100*(DataArray(:,7) - LowestVE)./(tempArthres - LowestDrive);
Comp_compareVE = 100*(DataArray_compare(:,7) - LowestVE_Compare)./(tempArthres_compare - LowestDrive_Compare);

[p2]=ranksum(CompVE,Comp_compareVE)
Comps = [CompVE Comp_compareVE];
[nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]

[p2]=ranksum(CompGGp,Comp_compareGGp)
Comps = [CompGGp Comp_compareGGp];
[nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]

%% 
NN=61;
%% Mixed Model Stats

X = size(VEdeciles,1)*size(VEdeciles,2)
Flow = [reshape(VEdeciles,X,1);reshape(VEdeciles_compare,X,1)];
Subj = [reshape(repmat([1:NN],10,1),X,1);reshape(repmat([1:NN],10,1),X,1)];
Drive = 0.01*([reshape(Vdrivedeciles,X,1);reshape(Vdrivedeciles_compare,X,1)] - 100);
Drive2 = ([reshape(Vdrivedeciles,X,1);reshape(Vdrivedeciles_compare,X,1)]);
GGt = [reshape(GGtdeciles,X,1);reshape(GGtdeciles_compare,X,1)];
GGp = [reshape(GGpdeciles,X,1);reshape(GGpdeciles_compare,X,1)];
Events = 1-[reshape(EdecilesMean,X,1);reshape(EdecilesMean_compare,X,1)];
REM = [ones(X,1);zeros(X,1)];
Flateral = 1-[reshape(repmat(Fsupine(:,6)',10,1),X,1);reshape(repmat(Fsupine_compare(:,6)',10,1),X,1)];
Decile = [repmat([1:10]',size(VEdeciles,2),1);repmat([1:10]',size(VEdeciles,2),1)] - 1;
DecileM = [repmat([1:10]',size(VEdeciles,2),1);repmat([1:10]',size(VEdeciles,2),1)] - 5.5;
Decile10 = [repmat([1:10]',size(VEdeciles,2),1);repmat([1:10]',size(VEdeciles,2),1)] - 10;
DecileN = nominal(Decile);

%AHIratio = [reshape(repmat(AHIratioT',10,1),X,1);reshape(repmat((AHIratioT*0)',10,1),X,1)];

linkpower=1/2; F.Link = @(mu) mu.^linkpower; F.Derivative= @(mu) linkpower*mu.^(linkpower-1); F.Inverse = @(mu) mu.^(1/linkpower); %link
AHI = F.Link([reshape(repmat(AHIrem',10,1),X,1);reshape(repmat((AHInrem)',10,1),X,1)]);

T1 = table(Subj,Flow,Drive,GGp,GGt,REM,Flateral,Decile,Events,DecileN,DecileM,Decile10,Drive2);
T1 = table(Subj,Flow,Drive,GGp,GGt,REM,Flateral,Decile,Events,DecileN,AHI,DecileM,Decile10,Drive2);

T2 = T1(T1.Decile==0,:);
T2.Drive1 = T1{T1.Decile==0,'Drive'}; %0.01*(T2.Drive-nanmean(T2.Drive(T2.REM==0)));
T2.Drive2 = T1{T1.Decile==1,'Drive'};
T2.Drive3 = T1{T1.Decile==2,'Drive'};
T2.DriveQ1 = (T2.Drive2 + T2.Drive3)/2;

T2.GGp1 = 0.01*(T2.GGp-nanmean(T2.GGp(T2.REM==0)));
T2.GGp2 = T1{T1.Decile==1,'GGp'};
T2.GGp3 = T1{T1.Decile==2,'GGp'};
T2.GGpQ1 = (T2.GGp2 + T2.GGp3)/2;

T2.Drive5 = T1{T1.Decile==4,'Drive'};
T2.Drive6 = T1{T1.Decile==5,'Drive'};
T2.Drive10 = T1{T1.Decile==9,'Drive'};
T2.DriveQ1 = (T2.Drive2 + T2.Drive3)/2;
T2.DriveM = (T2.Drive5 + T2.Drive6)/2;
T2.DriveMless1 = T2.DriveM - T2.Drive1;
T2.Drive10less1 = T2.Drive10 - T2.Drive1;
%standardize?
T2.Drive1c = (T2.Drive1-nanmean(T2.Drive1(T2.REM==0)));
T2.DriveMc = (T2.DriveM-nanmean(T2.DriveM(T2.REM==0)));
T2.Drive10c = (T2.Drive10-nanmean(T2.Drive10(T2.REM==0)));

T2.Drive1c = (T2.Drive1-nanmean(T2.Drive1(T2.REM==0)));
T2.DriveMc = (T2.DriveM-nanmean(T2.DriveM(T2.REM==0)));
T2.Drive10c = (T2.Drive10-nanmean(T2.Drive10(T2.REM==0)));
norm2SD = @(x) (x-nanmean(x))./(2*nanstd(x));
norm1SD = @(x) (x-nanmean(x))./(nanstd(x));

T2.AHI1SD = norm1SD(T2.AHI);
T2.Drive12SD = norm2SD(T2.Drive1);
T2.DriveM2SD = norm2SD(T2.DriveM);
T2.Drive102SD = norm2SD(T2.Drive10);
if 0
    T2.Drive2 = 0.01*(T2.Drive2-nanmean(T2.Drive2(T2.REM==0)));
    T2.Drive3 = 0.01*(T2.Drive3-nanmean(T2.Drive3(T2.REM==0)));
    T2.DriveQ1 = 0.01*(T2.DriveQ1-nanmean(T2.DriveQ1(T2.REM==0)));
    T2.DriveM = 0.01*(T2.DriveM-nanmean(T2.DriveM(T2.REM==0)));
    T2.Drive10 = 0.01*(T2.Drive10-nanmean(T2.Drive10(T2.REM==0)));
end

T2.ArThres = [arthres;arthres_compare];
T2.LG1 = [DataArray(:,1);DataArray_compare(:,1)];

%% Alt Models
%comment: adding Flateral to models has no impact

%VE-Vdrive curves are the same
mdl = fitglme(T1,'Flow ~ REM  + (1 | Subj)') %%Mediation: Flow is reduced in REM (10, 21%), but not after adjusting for Drive
mdl = fitglme(T1,'Flow ~ Drive + Drive^2 + (1 + Drive | Subj)') %curve
mdl = fitglme(T1,'Flow ~ REM  + Drive^2 + (1 + Drive | Subj)') %REM does not change baseline of the curve
mdl = fitglme(T1,'Flow ~ Drive*REM  + Drive^2 + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve

X = randomEffects(mdl)


%Simple %%%%%%%%%%%%%%%
mdl = fitglme(T1,'GGp ~ REM + (1  | Subj)')
mdl = fitglme(T1,'GGp ~ Drive + Drive^2 + (1 + Drive | Subj)')
mdl = fitglme(T1,'GGp ~ Drive + REM + Drive^2 + (1 + Drive | Subj)')
mdl = fitglme(T1,'GGp ~ Drive*REM + Drive^2 + (1 + Drive | Subj)')

%Simple %%%%%%%%%%%%%%%
mdl = fitglme(T1,'GGt ~ REM + (1 | Subj)')
mdl = fitglme(T1,'GGt ~ Drive + Drive^2 + (1 + Drive | Subj)')
mdl = fitglme(T1,'GGt ~ Drive + REM + Drive^2 + (1 + Drive | Subj)')
mdl = fitglme(T1,'GGt ~ Drive*REM + Drive^2 + (1 + Drive | Subj)')

%Is flow dependent on GG? P and T independently?
mdl = fitglme(T1,'Flow ~ GGp + GGt + (1 | Subj)')
mdl = fitglme(T1,'Flow ~ GGp + GGt + Drive + Drive^2 + (1 | Subj)')

%% Models
%comment: adding Flateral to models has no impact

%VE-Vdrive curves are the same
mdl = fitglme(T1,'Flow ~ REM  + (1 | Subj)') %%Mediation: Flow is reduced in REM (10, 21%), but not after adjusting for Drive
mdl = fitglme(T1,'Flow ~ Drive + Drive^2 + (1 | Subj)') %curve
mdl = fitglme(T1,'Flow ~ REM  + Drive^2 + (1 | Subj)') %REM does not change baseline of the curve
mdl = fitglme(T1,'Flow ~ Drive*REM  + Drive^2 + (1 | Subj)')  %REM does not change baseline or reduce the slope of the curve

X = randomEffects(mdl)


%Simple %%%%%%%%%%%%%%%
mdl = fitglme(T1,'GGp ~ REM + (1 | Subj)')
mdl = fitglme(T1,'GGp ~ Drive + Drive^2 + (1 | Subj)')
mdl = fitglme(T1,'GGp ~ Drive + REM + Drive^2 + (1 | Subj)')
mdl = fitglme(T1,'GGp ~ Drive*REM + Drive^2 + (1 | Subj)')

%Simple %%%%%%%%%%%%%%%
mdl = fitglme(T1,'GGt ~ REM + (1 | Subj)')
mdl = fitglme(T1,'GGt ~ Drive + Drive^2 + (1 | Subj)')
mdl = fitglme(T1,'GGt ~ Drive + REM + Drive^2 + (1 | Subj)')
mdl = fitglme(T1,'GGt ~ Drive*REM + Drive^2 + (1 | Subj)')

%Is flow dependent on GG? P and T independently?
mdl = fitglme(T1,'Flow ~ GGp + GGt + (1 | Subj)')
mdl = fitglme(T1,'Flow ~ GGp + GGt + Drive + Drive^2 + (1 | Subj)')

%% Not fully compelling, but works
%REM is associated with greater event likelihood, but a little less so
%after accounting for drive.
mdl = fitglme(T1,'Events ~ REM + (1 | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ REM + GGp + (1 | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ REM + GGt + (1 | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ REM + Drive^2 + (1 | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ REM + Drive^2 + GGt + (1 | Subj)','Link','logit')

%alt model:
mdl = fitglme(T1,'Events ~ REM + (1 | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ Drive^2 + (1 + Drive | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ REM + Drive^2 + (1 + Drive | Subj)','Link','logit')
mdl = fitglme(T1,'Events ~ REM + Drive^2 + REM*Drive + (1 + Drive | Subj)','Link','logit')


%% All first decile values are lowered in phasic REM
%Check these, missing some squared terms, should try to be consistent

mdlD = fitglme(T1,'Drive2 ~ Decile*REM + Decile^2 + (1 | Subj)') %[-22, -24]
mdlF = fitglme(T1,'Flow ~ Decile*REM + (1 | Subj)') %[-26, -32]
mdlP = fitglme(T1,'GGp ~ Decile*REM + (1 | Subj)') %[-5 ns, -20]
mdlT = fitglme(T1,'GGt ~ Decile*REM + (1 | Subj)') %[-8, -13]

mdlD = fitglme(T1,'Drive2 ~ DecileM*REM + DecileM^2 + (1 | Subj)') %[-45, -60]
mdlF = fitglme(T1,'Flow ~ DecileM*REM + (1 | Subj)') %[-10, -22]
mdlP = fitglme(T1,'GGp ~ DecileM*REM + (1 | Subj)') %[+10]
mdlT = fitglme(T1,'GGt ~ DecileM*REM + (1 | Subj)') %[-10]

mdlD = fitglme(T1,'Drive2 ~ Decile10*REM + (1 | Subj)') %[-70]
mdlF = fitglme(T1,'Flow ~ Decile10*REM + (1 | Subj)') %[+5]
mdlP = fitglme(T1,'GGp ~ Decile10*REM + (1 | Subj)') %[+25]
mdlT = fitglme(T1,'GGt ~ Decile10*REM + (1 | Subj)') %[-12]

%% Ignore for now
mdlD = fitglme(T2,'Drive1 ~ REM + (1 | Subj)') %[]
mdlP = fitglme(T2,'GGp ~ REM + (1 | Subj)') %[]

mdlD = fitglme(T2,'DriveQ1 ~ REM + (1 | Subj)') %[]
mdlP = fitglme(T2,'GGpQ1 ~ REM + Flateral + (1 | Subj)') %[]


%% AHI, nadir drive, and REM
mdlA = fitglme(T2,'AHI ~ REM + (1 | Subj)'); BackTransformModel(mdlA,F) %%%
mdlB = fitglme(T2,'AHI ~ Drive1c + (1 | Subj)'); BackTransformModel(mdlB,F) %%%
mdlC = fitglme(T2,'AHI ~ REM + Drive1c + (1 | Subj)'); BackTransformModel(mdlC,F) %%%


mdlC = fitglme(T2,'AHI ~ REM + DriveMc + (1 | Subj)'); BackTransformModel(mdlC,F) %%%
mdlC = fitglme(T2,'AHI ~ REM + Drive10c + (1 | Subj)'); BackTransformModel(mdlC,F) %%%


%% ArThres, add to T2:
mdlC = fitglme(T2,'AHI ~ REM + ArThres + (1 | Subj)'); BackTransformModel(mdlC,F) 
mdlC = fitglme(T2,'AHI ~ REM + LG1 + (1 | Subj)'); BackTransformModel(mdlC,F) 

mdlC = fitglme(T2,'AHI ~ REM + Drive1c + ArThres + (1 | Subj)'); BackTransformModel(mdlC,F) %%% Interesting. Low ArThres is bad, and slightly attenuates
mdlC = fitglme(T2,'AHI ~ REM + Drive1c + LG1 + (1 | Subj)'); BackTransformModel(mdlC,F) 

mdlC = fitglme(T2,'AHI ~ REM + Drive1c + ArThres + LG1 + (1 | Subj)'); BackTransformModel(mdlC,F) 
mdlC = fitglme(T2,'AHI ~ REM + Drive1c + ArThres + LG1 + (1 | Subj)'); BackTransformModel(mdlC,F) 


%% AHI, median drive, and REM
mdl = fitglme(T2,'AHI1SD ~ REM + (1 | Subj)')
mdl.Rsquared

mdlB = fitglme(T2,'AHI1SD ~ Drive12SD + (1 | Subj)')
mdlB.Rsquared
mdlC = fitglme(T2,'AHI1SD ~ REM + Drive12SD + (1 | Subj)')
mdlC.Rsquared
mdlB = fitglme(T2,'AHI1SD ~ DriveM2SD + (1 | Subj)')
mdlC = fitglme(T2,'AHI1SD ~ REM + DriveM2SD + (1 | Subj)')

mdlB = fitglme(T2,'AHI1SD ~ Drive102SD + (1 | Subj)')
mdlC = fitglme(T2,'AHI1SD ~ REM + Drive102SD + (1 | Subj)')

%%
figure(34)
plotregressionwithSEM(T2.AHI,100*(T2.Drive1+1))



%%
if 0
    mdl = fitglme(T2,'AHI ~ DriveQ1 + (1 | Subj)')
    mdl = fitglme(T2,'AHI ~ Drive2 + (1 | Subj)')
    mdl = fitglme(T2,'AHI ~ Drive3 + (1 | Subj)')
    mdl = fitglme(T2,'AHI ~ DriveM + (1 | Subj)')
    mdl = fitglme(T2,'AHI ~ DriveMless1 + (1 | Subj)')
    mdl = fitglme(T2,'AHI ~ Drive10less1 + (1 | Subj)')
    
    mdl = fitglme(T2,'AHI ~ REM + (1 | Subj)'); BackTransformModel(mdl,F) %%%
    mdl = fitglme(T2,'AHI ~ DriveQ1 + (1 | Subj)')
    
    mdl = fitglme(T2,'AHI ~ REM + DriveQ1 + (1 | Subj)'); BackTransformModel(mdl,F)
    mdl = fitglme(T2,'AHI ~ REM + DriveM + (1 | Subj)'); BackTransformModel(mdl,F)
    mdl = fitglme(T2,'AHI ~ REM + DriveMless1 + (1 | Subj)'); BackTransformModel(mdl,F)
    mdl = fitglme(T2,'AHI ~ REM + Drive10 + (1 | Subj)'); BackTransformModel(mdl,F)
    
    mdl = fitglme(T2,'AHI ~ REM + Drive1 + (1 | Subj)'); BackTransformModel(mdl,F) %%%
end

%%
I=I2;
figure(100); clf(100);
Nplots = 6;

subplot(1,Nplots,1);
X = [Vdrivedeciles_compare(1,:)',Vdrivedeciles(1,:)'] - Vdrivedeciles_compare(1,:)';
X = X(I,:);
DotAcrossPlot(X)
ylabel('\fontsize{13} \fontname{Arial Narrow} \Delta Nadir Drive [Withdrawal] (%)');

[~,p]=ttest2(Vdrivedeciles_compare(1,:)',Vdrivedeciles(1,:)');
[p2]=ranksum(Vdrivedeciles_compare(1,:)',Vdrivedeciles(1,:)');

subplot(1,Nplots,2);
X = [DataArrayIncl_compare(:,6),DataArrayIncl(:,6)] - DataArrayIncl_compare(:,6);
X = X(I,:);
DotAcrossPlot(X)
ylabel('\fontsize{13} \fontname{Arial Narrow} \Delta Baseline Ventilation [Collapsibility] (%)');

[~,p]=ttest2(DataArrayIncl_compare(:,6),DataArrayIncl(:,6));
[p2]=ranksum(DataArrayIncl_compare(:,6),DataArrayIncl(:,6));

subplot(1,Nplots,4);
X = ([Comp_compareVE,CompVE] - Comp_compareVE);
X = X(I,:);
DotAcrossPlot(X)
ylabel('\fontsize{13} \fontname{Arial Narrow} \Delta Flow-Drive Slope [Functional Responsiveness] (%)');

[~,p]=ttest2(Comp_compareVE,CompVE);
[p2]=ranksum(Comp_compareVE,CompVE);

subplot(1,Nplots,3);
X = [GGpassivep_compare,GGpassivep] - GGpassivep_compare;
X = X(I,:);
DotAcrossPlot(X)
ylabel('\fontsize{13} \fontname{Arial Narrow} \Delta Baseline EMGgg [Muscle Activity] (%)');


[~,p]=ttest2(GGpassivep_compare,GGpassivep);
[p2]=ranksum(GGpassivep_compare,GGpassivep);

subplot(1,Nplots,5);
X = ([Comp_compareGGp,CompGGp] - Comp_compareGGp);
X = X(I,:);
DotAcrossPlot(X)
ylabel('\fontsize{13} \fontname{Arial Narrow} \Delta EMGgg-Drive Slope [Responsiveness] (%) ');

[~,p]=ttest2(Comp_compareGGp,CompGGp);
[p2]=ranksum(Comp_compareGGp,CompGGp);

% 
% subplot(1,Nplots,6);
% X = [DataArrayIncl_compare(:,5),DataArrayIncl(:,5)] - DataArrayIncl_compare(:,5);
% X = X(I,:);
% DotAcrossPlot(X)
% ylabel('? Arousal Threshold (%)');

subplot(1,Nplots,6);
X = [AHInrem,AHIrem] - AHInrem;
X = X(I,:);
DotAcrossPlot(X)
ylabel('\fontsize{13} \fontname{Arial Narrow} \Delta Apnea-Hypopnea Index [OSA Severity] (events/h)');

[~,p]=ttest2(AHInrem,AHIrem);
[p2]=ranksum(AHInrem,AHIrem);

% subplot(1,Nplots,7);
% X = [DataArrayIncl_compare(:,1),DataArrayIncl(:,1)] - DataArrayIncl_compare(:,1);
% X = X(I,:);
% DotAcrossPlot(X)

GGpassivep = [GGdataInclNorm(:,1)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.
GGpassivep_compare = [GGdataIncl_compareNorm(:,1)]; %data in 0-10th ctile i.e. first decile, i.e. ~5th ctile.



% 
% CompGGt = (GGactivet - LowestGGt)./(tempArthres - LowestDrive);
% Comp_compareGGt = (GGactivet_compare - LowestGGt_Compare)./(tempArthres_compare - LowestDrive_Compare);
% 
% CompGGp = (GGactivep - LowestGGp)./(tempArthres - LowestDrive);
% Comp_compareGGp = (GGactivep_compare - LowestGGp_Compare)./(tempArthres_compare - LowestDrive_Compare);
% 
% CompVE = (DataArray(:,7) - LowestVE)./(tempArthres - LowestDrive);
% Comp_compareVE = (DataArray_compare(:,7) - LowestVE_Compare)./(tempArthres_compare - LowestDrive_Compare);

