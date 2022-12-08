function PlotThreeGroupIndividualData(h,Data,crit,crit2,crit3,FShrinkMarker)
N=length(Data);
% color3 = [0.1 0.7 0.2];
% color1 = [0.2 0.1 1];
% color2 = [1 0.3 0.2];
color3 = [0.5 0.5 0.5];
color2 = [0.2 0.1 1];
color1 = [1 0.3 0.2];

if ~exist('FShrinkMarker')
    FShrinkMarker=1;
end
    
figure(77777)
bw=0.05/(2^log10(N/100));
trans = 0.7 - 0.2* log10(N/100);
    trans(trans>1)=1;
    trans(trans<0.1)=0.1;
mn = 1-(1./(N.^0.5));
data = Data(crit);
catIdx = ones(length(data),1);
[~,dataOut,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
if exist('crit2')
    data = Data(crit2);
    catIdx = ones(length(data),1);
    [~,dataOut2,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
end
if exist('crit3')
    data = Data(crit3);
    catIdx = ones(length(data),1);
    [~,dataOut3,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
end
close(77777)


if h.Type=="figure"
    figure(h); 
else
    subplot(h); 
end
if 1 %clear axes
    hh=get(gca,'children'); delete(hh);
end

set(gcf,'color',[1 1 1]);

scatter(dataOut(:,1),dataOut(:,2),50*FShrinkMarker,color1,'filled','markerfacealpha',trans);
xlim([0.5 1.5]);
if exist('crit2')
hold on
scatter(dataOut2(:,1)+1,dataOut2(:,2),50*FShrinkMarker,color2,'filled','markerfacealpha',trans);
xlim([0.5 2.5]);
end
if exist('crit3')
hold on
scatter(dataOut3(:,1)+2,dataOut3(:,2),50*FShrinkMarker,color3,'filled','markerfacealpha',trans);
xlim([0.5 3.5]);
end
hold off

%set(gcf,'position',[ 991   545   249   433])
set(gca,'tickdir','out')

% ylabel('AHI')
%  set(gca,'xtick',[1 2 3],'xticklabels',{'1','2'},'fontsize',12);
 set(gca,'xtick',[1 2 3],'xticklabels',{'','',''},'fontsize',12);
% xtickangle(45);