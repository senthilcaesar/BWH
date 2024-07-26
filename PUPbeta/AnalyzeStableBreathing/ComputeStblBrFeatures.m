function EvtVEFtrs2 = ComputeStblBrFeatures(EvtVE)

%% Find start and end of average event based on VEmean
EnTime=EvtVE.EnTime;
VEmean = nanmean(EvtVE.VI);
% plotchanneltext = { ...
%     {'VI','VdriveEdiNorm'}, ...
%     {'GGpeak','GGtonic'}, ...
%     {'Arousal'} ...
%     }; %,'WakeSleepVe' ,{'DeltaPes'}     {'Dsat'} ...
% 
% C = { ...
%     {[0.2 0.2 1],[0.2 0.2 0.3]}, ...
%     {[0.95 0.5 0.1],[0.95 0.5 0.1]}, ...
%     {[0.4 0.95 0.3]} ...
%     };
%    {[0.8 0.1 0.1]} ...

thres=1;
temp = VEmean>thres;
I=find(diff(temp)==1);
[~,i]=min(abs(EnTime(I)));
I2=I(i);

deltaT = interp1(VEmean([I2 I2+1]),EnTime([I2 I2+1]),thres);

I3=find(diff(temp)==-1);%added 1 to get the next second, within the event
I3(I3>I2)=[];
I3=I3(end);

deltaT2 = interp1(VEmean([I3 I3+1]),EnTime([I3 I3+1]),thres);

if 1
    EnTime = EnTime - deltaT;
    deltaT2=deltaT2-deltaT;
    deltaT=0;
end

%% Find mean event depth
EventDepth = 1-mean(VEmean(EnTime>=deltaT2 & EnTime <=0));
EventDepth2Min = 1-min(VEmean(EnTime>=deltaT2 & EnTime <=0));

% Store in struct
EvtVEFtrs2.EventDepth = EventDepth;
EvtVEFtrs2.EventDepth2Min = EventDepth2Min;