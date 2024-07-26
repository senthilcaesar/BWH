function ReplaceHumanArousalTimingWithAuto(Evts,SigT)

clear ArTA ArTH

Time = SigT.Time;
EventsAr=0*Time;
Fs = 1./(Time(2)-Time(1));

figure (3333); clf(3333); plot(SigT.Time,SigT.EventsAr); hold on
plot(SigT.Time,SigT.EventsArWS,'k')

% auto scored
ArTA = Evts.EvtsArAuto.ArT;
X = ArTA.EventStart;
Y = ArTA.EventEnd;

% human scored
ArTH = Evts.ArT;
ArTH.TimingEdited=zeros(height(ArTH),1);
ArTH.EventStartNew=ArTH.EventStart;
ArTH.EventEndNew=ArTH.EventEnd;

for i=1:height(ArTH)
    A = ArTH.EventStart(i);
    B = ArTH.EventEnd(i);
    %find qualifying auto arousals (distance,d<5s)
    try
        clear DistancetoAr I ind val
        DistancetoAr=[abs((A-X)./Fs),abs((A-Y)./Fs),abs((B-X)./Fs),abs((B-Y)./Fs)];
        I=any(DistancetoAr<5,2);
        if nansum(I)>0
            if nansum(I)>1
                [val,ind]=min(min(DistancetoAr(I),[],2));
                ind=find(DistancetoAr==val,1,'first');
            else
                ind=I;
            end
        else
            ind =[];
        end
    end
    
    % old code
    %    I = ~(A>Y | B<X);
    %    ind = find(I==1);
    %ArTA.OverlapTemp=I;
    
    if ~isempty(ind)
        ArTH.TimingEdited(i)=1;
        ArTH.EventStartNew(i) = ArTA.EventStart(ind(1));
        ArTH.EventEndNew(i) = ArTA.EventEnd(ind(end));
    end
end

ArTH.EventDurationNew = ArTH.EventEndNew-ArTH.EventStartNew;

%% continuous Arousal Signal
for m=1:length(ArTH.EventStart) 
    lefti=round((ArTH.EventStart(m)-Time(1))*Fs)+1;
    righti=lefti+round((ArTH.EventDuration(m))*Fs);
    if lefti<1, lefti=1; end
    if righti>length(Time), righti=length(Time); end
    EventsAr(lefti:righti)=1;
end

Evts.ArT2 = EventRespSigtoRespT(EventsAr,Time);
[Evts.ArT,Evts.ArTinfo]=getArT(DataEventHypnog_Mat,ChannelsList,"EventsAr",Evts.Hypnogram,Evts.Hypnogram_t);

