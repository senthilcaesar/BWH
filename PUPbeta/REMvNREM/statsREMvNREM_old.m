
rowLow = 16; %contains values 0-20th centile, centered on 10th centile
rowMid = 66;
rowHigh = 116;
%note 1st value contains 0-5th centile, centered on 2.5th centile.

%% FlowVsDrive PhenoPlot Average
if settings.runcompare
    Iexclude_compare = DataNArray_compare(:,6)<2 | DataNArray(:,6)<2;
    Iexclude = DataNArray_compare(:,6)<2 | DataNArray(:,6)<2;
    VEdeciles_compare(:,Iexclude_compare==1)=NaN;
    Vdrivedeciles_compare(:,Iexclude_compare==1)=NaN;
    DataArrayIncl_compare=DataArray_compare;
    DataArrayIncl_compare(Iexclude_compare==1,:)=NaN;
    GGpdeciles_compare(:,Iexclude_compare==1)=NaN;
    GGtdeciles_compare(:,Iexclude_compare==1)=NaN;
else
    Iexclude = DataNArray(:,6)<2;
end
VEdeciles(:,Iexclude==1)=NaN;
Vdrivedeciles(:,Iexclude==1)=NaN;
DataArrayIncl=DataArray;
    DataArrayIncl(Iexclude==1,:)=NaN;
    GGpdeciles(:,Iexclude_compare==1)=NaN;
    GGtdeciles(:,Iexclude_compare==1)=NaN;
    
    GGdataIncl = GGdata;
        GGdataIncl(Iexclude==1,:)=NaN;
    GGdataIncl_compare = GGdata_compare;
        GGdataIncl_compare(Iexclude==1,:)=NaN;

    if 1 %normalize by VpassiveCompare
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
subplot(2,1,1);
if settings.runcompare
    GGpdecilesUpper_compare = prctile(GGpdeciles_compare',75)';
    GGpdecilesLower_compare = prctile(GGpdeciles_compare',25)';
    VdrivedecilesMedian_compare = prctile(Vdrivedeciles_compare',50)';
    GGpdecilesMedian_compare = prctile(GGpdeciles_compare',50)';
   
    set(gcf,'color',[1 1 1]);
    if 1
        plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
        hold('on')
        %plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
        box('off');
        set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[GGpdecilesUpper_compare;flipud(GGpdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,GGpdecilesMedian_compare,'k','linewidth',2);
    
    %VpassiveMedian_compare = nanmedian(DataArrayIncl_compare(:,6));
    
    [~,ii] = unique(VdrivedecilesMedian_compare);
    GGppassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),GGpdecilesMedian_compare(ii),100,'linear');
    
    arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
%     VactiveMedian_compare = nanmedian(DataArrayIncl_compare(:,7));
    
    GGpactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),GGpdecilesMedian_compare(ii),arthresMedian_compare,'linear');
    plot(100,GGppassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,GGpactiveMedian_compare)],'g'); %Arthres line
    plot(arthresMedian_compare,GGpactiveMedian_compare,'r.','markersize',25);  %Vactive
    
end

GGpdecilesUpper = prctile(GGpdeciles',75)';
GGpdecilesLower = prctile(GGpdeciles',25)';
VdrivedecilesMedian = prctile(Vdrivedeciles',50)';
GGpdecilesMedian = prctile(GGpdeciles',50)';

    if ~settings.runcompare
    plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
    hold('on')
    %plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGpdecilesUpper;flipud(GGpdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,GGpdecilesMedian,'k','linewidth',2);

    
%     VpassiveMedian = nanmedian(DataArrayIncl(:,6));
    arthresMedian = nanmedian(DataArrayIncl(:,5));
%     VactiveMedian = nanmedian(DataArrayIncl(:,7));
    
    [~,ii] = unique(VdrivedecilesMedian);
    GGppassiveMedian=interp1(VdrivedecilesMedian(ii),GGpdecilesMedian(ii),100,'linear');
    GGpactiveMedian = interp1(VdrivedecilesMedian(ii),GGpdecilesMedian(ii),arthresMedian,'linear');

    

    plot(100,GGppassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian arthresMedian],[0 GGpactiveMedian],'g'); %Arthres line
    plot(arthresMedian,GGpactiveMedian,'r.','markersize',25);  %Vactive

    saveas(10,[settings.directoryout 'Fig9_' settings.Comparename '_' num2str(i) '_s' num2str(settings.selectstate) '.png'])
    
    ylim([0 max(max(GGpdecilesUpper),max(GGpdecilesUpper_compare))])
    
subplot(2,1,2);
if settings.runcompare
    GGtdecilesUpper_compare = prctile(GGtdeciles_compare',75)';
    GGtdecilesLower_compare = prctile(GGtdeciles_compare',25)';
    VdrivedecilesMedian_compare = prctile(Vdrivedeciles_compare',50)';
    GGtdecilesMedian_compare = prctile(GGtdeciles_compare',50)';
   
    set(gcf,'color',[1 1 1]);
    if 1
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
    
    GGtactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),GGtdecilesMedian_compare(ii),arthresMedian_compare,'linear');
    plot(100,GGtpassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,GGtactiveMedian_compare)],'g'); %Arthres line
    plot(arthresMedian_compare,GGtactiveMedian_compare,'r.','markersize',25);  %Vactive
    
end

GGtdecilesUpper = prctile(GGtdeciles',75)';
GGtdecilesLower = prctile(GGtdeciles',25)';
VdrivedecilesMedian = prctile(Vdrivedeciles',50)';
GGtdecilesMedian = prctile(GGtdeciles',50)';

    if ~settings.runcompare
    plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
    hold('on')
    %plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[GGtdecilesUpper;flipud(GGtdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,GGtdecilesMedian,'k','linewidth',2);

    
%     VpassiveMedian = nanmedian(DataArrayIncl(:,6));
     arthresMedian = nanmedian(DataArrayIncl(:,5));
%     VactiveMedian = nanmedian(DataArrayIncl(:,7));
    
[~,ii] = unique(VdrivedecilesMedian);
    GGtpassiveMedian=interp1(VdrivedecilesMedian(ii),GGtdecilesMedian(ii),100,'linear');
    GGtactiveMedian = interp1(VdrivedecilesMedian(ii),GGtdecilesMedian(ii),arthresMedian,'linear');

    

    plot(100,GGtpassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian arthresMedian],[0 GGtactiveMedian],'g'); %Arthres line
    plot(arthresMedian,GGtactiveMedian,'r.','markersize',25);  %Vactive

    saveas(10,[settings.directoryout 'Fig10_' settings.Comparename '_s' num2str(settings.selectstate) '.png'])
    
    ylim([0 max(max(GGpdecilesUpper),max(GGpdecilesUpper_compare))])
        
%%    Ventilation
figure(9); clf(9);
if settings.runcompare
    VEdecilesUpper_compare = prctile(VEdeciles_compare',75)';
    VEdecilesLower_compare = prctile(VEdeciles_compare',25)';
    VdrivedecilesMedian_compare = prctile(Vdrivedeciles_compare',50)';
    VEdecilesMedian_compare = prctile(VEdeciles_compare',50)';
   
    set(gcf,'color',[1 1 1]);
    if 1
        plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
        hold('on')
        plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
        box('off');
        set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian_compare;flipud(VdrivedecilesMedian_compare)],[VEdecilesUpper_compare;flipud(VEdecilesLower_compare)],color2,'edgecolor','none','facealpha',0.5);
    plot(VdrivedecilesMedian_compare,VEdecilesMedian_compare,'k','linewidth',2);
    
    %VpassiveMedian_compare = nanmedian(DataArrayIncl_compare(:,6));
    [~,ii] = unique(VdrivedecilesMedian_compare);
    VpassiveMedian_compare=interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),100,'linear');
    
    arthresMedian_compare = nanmedian(DataArrayIncl_compare(:,5));
    VactiveMedian_compare = nanmedian(DataArrayIncl_compare(:,7));
    
    VactiveMedian_compare = interp1(VdrivedecilesMedian_compare(ii),VEdecilesMedian_compare(ii),arthresMedian_compare,'linear');
    plot(100,VpassiveMedian_compare,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian_compare arthresMedian_compare],[0 max(arthresMedian_compare,VactiveMedian_compare)],'g'); %Arthres line
    plot(arthresMedian_compare,VactiveMedian_compare,'r.','markersize',25);  %Vactive
    
end

VEdecilesUpper = prctile(VEdeciles',75)';
VEdecilesLower = prctile(VEdeciles',25)';
VdrivedecilesMedian = prctile(Vdrivedeciles',50)';
VEdecilesMedian = prctile(VEdeciles',50)';

    if ~settings.runcompare
    plot([100 100],[0 180],'--','color',[0.7 0.7 0.7]);
    hold('on')
    plot([0 300],[100 100],'--','color',[0.7 0.7 0.7]);
    box('off');
    set(gca,'tickdir','out')
    end
    fill([VdrivedecilesMedian;flipud(VdrivedecilesMedian)],[VEdecilesUpper;flipud(VEdecilesLower)],color1,'edgecolor','none','facealpha',0.5);
plot(VdrivedecilesMedian,VEdecilesMedian,'k','linewidth',2);

    
    VpassiveMedian = nanmedian(DataArrayIncl(:,6));
    arthresMedian = nanmedian(DataArrayIncl(:,5));
    VactiveMedian = nanmedian(DataArrayIncl(:,7));
    
    [~,ii] = unique(VdrivedecilesMedian);
    VpassiveMedian=interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),100,'linear');
    VactiveMedian = interp1(VdrivedecilesMedian(ii),VEdecilesMedian(ii),arthresMedian,'linear');

    

    plot(100,VpassiveMedian,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthresMedian arthresMedian],[0 max(arthresMedian,VactiveMedian)],'g'); %Arthres line
    plot(arthresMedian,VactiveMedian,'r.','markersize',25);  %Vactive

    saveas(9,[settings.directoryout 'Fig9_' settings.Comparename '_s' num2str(settings.selectstate) '.png'])
        
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
%% 
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

[nanmean(LowestGGps) nanstd(LowestGGps) nanmedian(LowestGGps)]

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
arthres = DataArray(:,7);
arthres_compare = DataArray_compare(:,7);
[p2]=ranksum(arthres,arthres_compare)
arthress = [arthres arthres_compare];
[nanmean(arthress) nanstd(arthress) nanmedian(arthress) tsnaniqr(arthress)]


Comp = (GGactivet - LowestGGt)./(DataArray(:,5) - LowestDrive);
Comp_compare = (GGactivet_compare - LowestGGt_Compare)./(DataArray_compare(:,5) - LowestDrive_Compare);

Comp = (GGactivep - LowestGGp)./(DataArray(:,5) - LowestDrive);
Comp_compare = (GGactivep_compare - LowestGGp_Compare)./(DataArray_compare(:,5) - LowestDrive_Compare);

Comp = (DataArray(:,7) - LowestVE)./(DataArray(:,5) - LowestDrive);
Comp_compare = (DataArray_compare(:,7) - LowestVE_Compare)./(DataArray_compare(:,5) - LowestDrive_Compare);

[p2]=ranksum(Comp,Comp_compare)
Comps = [Comp Comp_compare];
[nanmean(Comps) nanstd(Comps) nanmedian(Comps) tsnaniqr(Comps)]

