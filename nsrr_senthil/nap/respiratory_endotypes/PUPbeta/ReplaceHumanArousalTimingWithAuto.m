function [Evts,SigT]=ReplaceHumanArousalTimingWithAuto(Evts,SigT)

clear ArTA ArTH

Time = SigT.Time;
EventsAr=0*Time;
Fs = 1./(Time(2)-Time(1));



% auto scored
ArTA = Evts.EvtsArAuto.ArT;
X = ArTA.EventStart;
Y = ArTA.EventEnd;

% human scored
ArTH = Evts.ArT;
ArTH.TimingEdited=zeros(height(ArTH),1);
ArTH.EventDurationOrig = ArTH.EventDuration;
ArTH.EventEndOrig = ArTH.EventEnd;
ArTH.EventStartOrig = ArTH.EventStart;

for i=1:height(ArTH)
    A = ArTH.EventStartOrig(i); %original manual
    B = ArTH.EventEndOrig(i);
    %find qualifying auto arousals (distance,d<5s)
    try
        clear DistancetoAr I ind val
        DistancetoAr1=[abs((A-X)),abs((A-Y)),abs((B-X)),abs((B-Y))];
        DistancetoAr = min(DistancetoAr1,[],2);
        Overlap = (X>A & X<B) | (Y>A & Y<B) | (A>X & A<Y) | (B>X & B<Y);
        DistancetoAr(Overlap==1)=0;
        
        I=find(any(DistancetoAr<5,2));
        if length(I)>0
            if length(I)>1
                [val,~]=min(min(DistancetoAr(I,:),[],2)); % min(DistancetoAr(I,:),[],2) will provide min in each row 
                ind=find(DistancetoAr==val,1,'first'); %id row of selected arousal
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
        ArTH.EventStart(i) = ArTA.EventStart(ind(1));
        ArTH.EventEnd(i) = ArTA.EventEnd(ind(end));
    end
end

ArTH.EventDuration = ArTH.EventEnd-ArTH.EventStart;

%%


%%
Evts.ArTOrig = Evts.ArT;
Evts.ArTinfoOrig=Evts.ArTinfo;
SigT.EventsArOrig = SigT.EventsAr;

%% continuous Arousal Signal
for m=1:length(ArTH.EventStart) 
    lefti=round((ArTH.EventStart(m)-Time(1))*Fs)+1;
    righti=lefti+round((ArTH.EventDuration(m))*Fs);
    if lefti<1, lefti=1; end
    if righti>length(Time), righti=length(Time); end
    EventsAr(lefti:righti)=1;
end

SigT.EventsAr = EventsAr;
Evts.ArT = EventRespSigtoRespT(EventsAr,Time); %oldone
[Evts.ArT,Evts.ArTinfo]=getArT(SigT,SigT.Properties.VariableNames,"EventsAr",Evts.Hypnogram,Evts.Hypnogram_t);

figure (3333); clf(3333); 
plot(SigT.Time,[SigT.EventsArOrig SigT.EventsAr-1.1 SigT.EventsArWS-2.2]); 
