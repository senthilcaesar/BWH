function [rho, pval]=DoQuickScatter(FL, F1, F2, n, str1, str2, str3)
OutlierThreshold = inf;%1000; % should use some percentage of median signal
ScatterData1 = table2array(FL(:,F1)); 
Scat1Median=abs(nanmedian(ScatterData1));
Scat1Threshold = Scat1Median+(Scat1Median*OutlierThreshold/100);
%ScatterData2 = data_all;
ScatterData2 = table2array(FL(:,F2)); 
Scat2Median=abs(nanmedian(ScatterData2));
Scat2Threshold = Scat2Median+(Scat2Median*OutlierThreshold/100);
ScatterData1(abs(ScatterData1)>Scat1Threshold)=NaN;
ScatterData2(abs(ScatterData2)>Scat2Threshold)=NaN;
% mods
% if n == [15]
%     %invert
%     % 15 is Peak and Mid insp flow ratios
%     ScatterData2 = 1./ScatterData2;
% end
ax = figure(200+n); clf(figure(200+n)); %figure()
scatter(ScatterData1, ScatterData2, 'r.'); hold on;
lsline; xlabel(str1); ylabel(str2); title(str3);
% plot line of identity
refline(1,0);
axis('square');
[rho, pval] = corr(ScatterData1, ScatterData2);
end


