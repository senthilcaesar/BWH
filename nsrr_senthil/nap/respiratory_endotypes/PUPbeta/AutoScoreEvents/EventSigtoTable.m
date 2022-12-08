clear I I1 I2 I1N I2N
I = diff([NaN;1*(ApneaHypopnea_>0)]);
I1 = find(I==1);
I2 = find(I==-1);
[I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
evntdur=(I2N-I1N)*(1/ChannelsFs);    


EventStart=EventStart(:); 
    EventDuration=EventDuration(:);
    
    clear Evts
    EventCodesList=EventCodesList(:);
    Evts.EventCodesList = EventCodesList;
    Evts.EventCategoriesLabels = EventCategoriesLabels;
    Evts.Table1 = table(EventCodesList,EventStart,EventDuration);
    
    if ~isempty(EventCodesList)
    clear match
    for m=1:length(EventCategoriesLabels)% each column is an event type
        if EventExact1orStartsWith2(m)==1 %exact
            match(:,m)=string(EventCategoriesLabels(m))==string(EventCodesList');
        elseif EventExact1orStartsWith2(m)==2 %startswith
            match(:,m)=startsWith(string(EventCodesList'),string(EventCategoriesLabels(m)));
        end
    end
    temp = match.*EventCategoriesCodes;
    temp(temp==0)=NaN;
    EventCodes=min(temp'); %just in case two criteria are detected, will leave the lowest code value in list
    
    else
        EventCodes=[];
        disp(['Events Empty']);
    end
    Evts.Table1.EventCodes=EventCodes(:);
    
    %check for events that should be in the list
    if ~isempty(EventCodesList(isnan(EventCodes)))
        disp(['Events found but unused: ' strjoin(unique(EventCodesList(isnan(EventCodes))),', ')]);
    else
        disp(['Events OK']);
    end

    %     for m=1:length(EventStart) %
    %         temp = strcmp(EventCategoriesLabels,EventCodesList(m));
    %         if sum(temp)>0
    %             EventCodes(m)=EventCategoriesCodes(temp);
    %         end
    %     end
    I=isnan(EventCodes);
    EventCodes(I)=[]; 
    EventStart(I)=[]; 
    EventDuration(I)=[];
    
    %remove very short resp events
    I=EventDuration(:)<2&EventCodes(:)>1;
    EventCodes(I)=[]; 
    EventStart(I)=[]; 
    EventDuration(I)=[];  
    
    if sum(EventCodes==1)==0
        displaytext=['Warning: No scored arousals'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    if std(EventDuration(EventCodes==1))<1
        displaytext=['Warning: Arousals are all similar durations'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    %initialize
    EventsResp=0*Time;
    EventsAr=0*Time;
    for m=1:length(EventStart) %Arousals
        lefti=round((EventStart(m)-StartTime)*F_samp)+1;
        righti=lefti+round((EventDuration(m))*F_samp);
        if lefti<1, lefti=1; end
        if righti>length(Time), righti=length(Time); end
        if EventCodes(m)>1
            EventsResp(lefti:righti)=EventCodes(m);
        elseif EventCodes(m)==1
            EventsAr(lefti:righti)=1;
        end
    end
    
    if 1 %set NaN, others, to 8
        Hypnogram(isnan(Hypnogram))=8; %no stage code
        Hypnogram(Hypnogram<0)=8; %no stage code
    end
    Hypnogram_t = StartTime+starttimediff+(0:30:(30*(length(Hypnogram)-1)))';
    Epochs = interp1(Hypnogram_t+15,Hypnogram,Time,'nearest','extrap');
    
    Evts.Hypnogram = Hypnogram;
    Evts.Hypnogram_t = Hypnogram_t;
%     if 1
%         Epochs(isnan(Epochs))=8;
%     end
    %figure(); plot(Time,Epochs);
    
    ChannelsList = [ChannelsList,{'Epochs','EventsAr','EventsResp'}];
    Channels_Fs = [Channels_Fs; F_samp;F_samp;F_samp];