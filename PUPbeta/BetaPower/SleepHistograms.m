%% Sleep Histograms
ArTemp = WStbl.EventsAr;
logit = @(p) log(p./(1-p));

Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
Edges(end)=Inf;

[h1,edges] = histcounts((XvalH(WStbl.Epochs==4&ArTemp>0.95&~WStbl.ExcludeAR)),Edges); 
[h2,edges] = histcounts((XvalH(WStbl.Epochs==2&ArTemp<0.05&~WStbl.ExcludeAR)),Edges);
[h3,edges] = histcounts((XvalH(WStbl.Epochs==1&ArTemp<0.05&~WStbl.ExcludeAR)),Edges);
[h4,edges] = histcounts((XvalH(WStbl.Epochs==0&ArTemp<0.05&~WStbl.ExcludeAR)),Edges);
[h5,edges] = histcounts((XvalH(WStbl.Epochs==3&ArTemp<0.05&~WStbl.ExcludeAR)),Edges);
[hAr,edges] = histcounts((XvalH(WStbl.Epochs<4&ArTemp>0.95&~WStbl.ExcludeW)),Edges); %AR

area = sum([h1,h2,h3,h4,h5])*dStep;
h1=h1/area;
h2=h2/area;
h3=h3/area;
h4=h4/area;
h5=h5/area;
hAr=hAr/area;

bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1); hold('on');
bar(Centers,h2,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h3,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h4,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h5,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
stairs(Centers-dStep/2,hAr,'g');

h=legend('W','N1','N2','N3','REM','AR');
set(h,'box','off');
box('off')
set(gca,'tickdir','out')

ax1=gca
ylim([0 0.1]);
xticks = get(gca,'xtick');
xlims = get(gca,'xlim');
set(gca,'xtick',xticks,'xlim',xlims);

ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
set(ax2,'xtick',xticks,'xlim',xlims,'xticklabels',round(logitinverse(xticks),4),'xcolor',[0 0 0],'tickdir','out','ycolor',[1 1 1]);

ylim([0 0.1]);