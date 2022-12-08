
chArPr = find(strcmp(ChannelsList,'ArPr'));
chWPr = find(strcmp(ChannelsList,'WPr'));
chArIntensity = find(strcmp(ChannelsList,'ArIntensity'));
chArIntensityOrig = find(strcmp(ChannelsList,'ArIntensityOrig'));
chEpochs = find(strcmp(ChannelsList,'Epochs'));
chEventsAr = find(strcmp(ChannelsList,'EventsAr'));
chHR = find(strcmp(ChannelsList,'HR'));

I = find(DataEventHypnog_Mat(:,chEpochs)<4 & DataEventHypnog_Mat(:,chEventsAr)==1);
%I = find(DataEventHypnog_Mat(:,chEpochs)<4 );

I = I(1:100:end);
PlotOption = 2; 
switch PlotOption
    case 1
        X = logit(DataEventHypnog_Mat(I,chArPr));
        Y = DataEventHypnog_Mat(I,chArIntensity);
    case 1.5
        X = logit(DataEventHypnog_Mat(I,chWPr));
        Y = DataEventHypnog_Mat(I,chArIntensity);
    case 1.8
        X = logit(DataEventHypnog_Mat(I,chArPr));
        Y = DataEventHypnog_Mat(I,chArIntensityOrig);
    case 2
        X = (DataEventHypnog_Mat(I,chHR));
        X(X<45)=NaN;
        Y = DataEventHypnog_Mat(I,chArIntensity);
    case 3
        X = (DataEventHypnog_Mat(I,chHR));
        X(X<45)=NaN;
        Y = logit(DataEventHypnog_Mat(I,chArPr));
    case 3.5
        X = (DataEventHypnog_Mat(I,chHR));
        Y = logit(DataEventHypnog_Mat(I,chWPr));
end

if 1 %swap X and Y
    temp=X;
    X=Y;
    Y=temp;
end

figure(67); clf(67);
scatter(X,Y,5,[0.8 0.2 0.1],'filled','markerfacealpha',0.2)

if 1 %addQuantiles
    Nciles=25;
    clear ybindata
    for j=1:Nciles
        lower = (j-1)*(100/Nciles); upper = j*(100/Nciles); if upper>100, upper=100; end
        I = (X>=prctile(X,lower)&X<=prctile(X,upper));
        ybindata(j,:)=[nanmedian(Y(I)),prctile(Y(I),25),prctile(Y(I),75),sum(I),nanmedian(X(I))];
    end
    hold('on')
    fill([ybindata(:,5);flipud(ybindata(:,5))],[ybindata(:,2);flipud(ybindata(:,3))],[0.9500 0.2000  0.0500],'edgecolor','none','facealpha',0.5);
    plot(ybindata(:,5),ybindata(:,1),'k-','linewidth',2);
end

I2 = ~isnan(X)&~isnan(Y);
corrXY = corr(X(I2),Y(I2),'Type','Spearman')

