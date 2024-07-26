function plotFeatStats(FeatStatsTable, dAHIper, LabelsSub, stat, FeatName, SubsIn)

figure('pos', [-934 48 819 910]);
subplot(2,1,1)
boxplot(FeatStatsTable.(stat), LabelsSub, 'Symbol', '')
hold on
scatter(LabelsSub+1, FeatStatsTable.(stat))
title(stat)
xticklabels({'NR', 'R'})
ylabel(FeatName)
set(gca,'Box','On','FontSize',12)

subplot(2,1,2)
scatter(FeatStatsTable.(stat)(LabelsSub == 0), dAHIper(LabelsSub == 0), 'rx')
hold on
scatter(FeatStatsTable.(stat)(LabelsSub == 1), dAHIper(LabelsSub == 1),'bo')

for ii = 1:length(SubsIn)
    text(FeatStatsTable.(stat)(ii), dAHIper(ii), num2str(SubsIn(ii)))
end

title(stat)
xlabel(FeatName)
ylabel('Change in AHI')
set(gca,'Box','On','FontSize',12)