function RespT = EventRespSigtoRespT(EventsResp,Time)
%Will also work for 
%Evts.ArT = EventRespSigtoRespT(EventsAr,Time)
%Evts.EvtsAuto.ArT = EventRespSigtoRespT(EventsArWS,Time)
%Evts.EvtsAuto.ArTB = EventRespSigtoRespT(EventsArWSB,Time)

if istable(Time)
Time=Time{:,:};
end

if istable(EventsResp)
  EventsResp=EventsResp{:,:};
end

dt = Time(2) - Time(1);

starttimesi=1+find([diff(EventsResp)>0]); %start is first sample inside event
    starttimesi(EventsResp(starttimesi-1)>0)=[]; %removes events if the pre-event sample is not zero [fix earlier]
endtimesi=1+find([diff(EventsResp)<0]); %end is first sample after event
    endtimesi(EventsResp(endtimesi)>0)=[];
[starttimesi,endtimesi] = TidyStartEndEventList(starttimesi,endtimesi,length(EventsResp));
% 
% figure(1); clf(1);
% plot(EventsResp)

EventStart = Time(starttimesi);
EventEnd = Time(endtimesi);

EventDuration = EventEnd - EventStart;
EventCodes = EventsResp(starttimesi);

RespT = table(EventCodes,EventStart,EventDuration,EventEnd,starttimesi,endtimesi);
RespT = sortrows(RespT,'EventStart');
