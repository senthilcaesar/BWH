plotchanneltext = { ...
    {'Epochs','EventsResp'},{'Flow'},{'SpO2'},{'EventsAr'},{'EventsArWS'},{'WPr'},{'WSBalance'},{'Pbetalogfilt','Pdeltalogfilt'} };

run PlotXHzData

WSBalance=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'WSBalance')==1));
Epochs=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Epochs')==1));

meanWSBalance = nanmean(WSBalance(Epochs<4))

