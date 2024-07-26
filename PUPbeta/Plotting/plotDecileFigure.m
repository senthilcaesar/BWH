
%% decile function
function [ydataAverage,ydataUpper,ydataLower,xdataAverage]= plotDecileFigure(xdata,ydata,colorX, averagingMethod,plotVpassive,plotVactive,plotArThline, ArTh, Nfigure, holdonfig)
%Function to plot nDecile endotype plots;
%Requires data from summary analysis in form xdata=Vdrivedeciles,
%ydata=GGpdeciles or VEdeciles etc
%colorX is color for fill region on endoplot figure
%ArTh is group median arousal threshold (from DataArray(:,5) in SummaryAnanlysis output)
%Averaging method 1=median,IQR; 2=mean,95%CI
%can choose whether to plot dots for Vpassive, Vactive and ArTh line
%Nfigure is to assign figure number
%holdonfig=1 allows to plot additional data for comparisons

%for example: plotDecileFigure(Vdrivedeciles,GGpdeciles,color3, 1,1,1, ArTh, 100, 1, 1)

%To add the colorbar for event liklihood to endoplot figure, run
%DisEventLikelihoodOnEndogram.m after
%DispEventLikelihoodOnEndogram(ydataAverage,xdataAverage,EdecilesMeanMedian,plotColorBar)

if ~exist('plotVpassive')|isempty('plotVpassive')
    plotVpassive=0
end

if ~exist('plotVactive')|isempty('plotVactive')
    plotVactive=0
end
if ~exist('plotArThline')|isempty('plotArThline')
    plotArThline=0
end
 if ~exist('Nfigure')|isempty('Nfigure')
     Nfigure=101 %assign new figure number
 end
 
 if ~exist('holdonfig')| isempty('holdonfig')
 holdonfig=0
 end
 
 if ~exist('averagingMethod')|isempty('averagingMethod')
     averagingMethod=1 %set to medians, IQR
 end

upper95CI = @(x) nanmean(x) + 1.96*nanstd(x)./sum(~isnan(x)).^0.5;
lower95CI = @(x) nanmean(x) - 1.96*nanstd(x)./sum(~isnan(x)).^0.5;
upperIQR = @(x) prctile2(x,75);
lowerIQR = @(x) prctile2(x,25);

if averagingMethod==1   
    upperlimit = upperIQR;
    lowerlimit = lowerIQR;
    average = @(x) nanmedian(x);
elseif averagingMethod==2
    upperlimit = upper95CI;
    lowerlimit = lower95CI;
    average = @(x) nanmean(x);
end

ydataUpper=upperlimit(ydata');
ydataAverage=average(ydata');
ydataLower=lowerlimit(ydata');
xdataAverage=average(xdata');



%% GGp
figure(Nfigure)
if holdonfig==0
clf(Nfigure);
end
    set(gcf,'color',[1 1 1]);
    if 1
        hold('on')
        box('off');
        set(gca,'tickdir','out')
    end
    
    fill([xdataAverage';flipud(xdataAverage')],[ydataUpper';flipud(ydataLower')], colorX,'edgecolor','none','facealpha',0.5);
    plot(xdataAverage,ydataAverage,'k','linewidth',2);
    
      if plotVpassive==1
    [~,ii] = unique(ydataAverage);
    yVpassive=interp1(xdataAverage(ii),ydataAverage(ii),100,'linear'); %find Vpassive
    plot(100,yVpassive,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
      end
      
      if plotArThline==1
          yVactive=interp1(xdataAverage(ii),ydataAverage(ii),ArTh,'linear');
          plot([ArTh ArTh],[0 yVactive],'g')
      end
      
      if plotVactive==1
         plot(ArTh,yVactive,'r.','markersize',25); %Vpassive
      end
          
if 0
    plot([100 100],[0 15],'--','color',[0.7 0.7 0.7]);
    hold('on')
    box('off');
    set(gca,'tickdir','out')
end



end

%saveas(10,[settings.directoryout 'Fig10_' settings.Comparename '_s' num2str(settings.selectstate) '.png'])
