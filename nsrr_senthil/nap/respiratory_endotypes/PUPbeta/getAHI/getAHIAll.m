function [Evts]=getAHIAll(SigT,ROImask,Evts,ChannelsList,settings,EventsChannelName,ArChannelName)

if ~exist('EventsChannelName')
    EventsChannelName = 'EventsResp';
end
if ~exist('ArChannelName')
    ArChannelName = 'EventsAr';
end

if isempty(ROImask)
    ROImask=ones(size(SigT(:,1)));
end

Time=SigT.Time;


SpO2ch=find(strcmp(ChannelsList,'SpO2')==1);
if ~isempty(SpO2ch)
    SpO2=SigT.SpO2;
else
    SpO2=NaN*Time;
end
SpO2Available = ~isnan(SpO2)*1;

            
[AHIdata{1},~]=getAHI(SigT,ROImask,ChannelsList,EventsChannelName);
%displaytext=['Total AHI: ' num2str(AHIdata{1}(58),4) ' per hr'];
%disp(displaytext);
%             set(handletext,'String',displaytext); drawnow;
%displaytext=['NREM supine AHI: ' num2str(AHIdata{1}(16),4) ' per hr'];
%disp(displaytext);
%             set(handletext,'String',displaytext); drawnow;

try %only works when Evts.RespT already exists
    I = Evts.RespT.starttimesi > size(SigT,1);
    if sum(I)>0
        disp('warning: events scored after end of recording');
    end
    Evts.RespT(I,:)=[];
catch me
    if settings.verbose
        disp('Notice: Evts.RespT not found');
        % disp(me.messge); % uncomment this if you want a little more
        % disp(me.getReport); % uncomment this if you want a lot more
    end
end

try %only works when Evts.ArT already exists
    I = Evts.ArT.starttimesi > size(SigT,1);
    if sum(I)>0
        disp('warning: arousals scored after end of recording');
    end
    Evts.ArT(I,:)=[];
catch me
    if settings.verbose
        disp('Notice: Evts.ArT not found');
        % disp(me.messge); % uncomment this if you want a little more
        % disp(me.getReport); % uncomment this if you want a lot more
    end
end

%             if 0 %% New process plan

Evts = EventRespTable(SigT,Evts,ChannelsList,EventsChannelName,ArChannelName);
% addding position codes for each events depending on protocol
[Position,PositionRaw] = getPos(SigT,ChannelsList,settings);
Evts.RespT.Position=Position(Evts.RespT.starttimesi);
Evts.RespT.PositionRaw=PositionRaw(Evts.RespT.starttimesi);

%% all events
AHI=struct();
EpochsXHz=SigT.Epochs;
Time=SigT.Time;
[AHI.Total,Evts,rowNames]=getAHIEvtSubset(Evts,EpochsXHz,Time,Position,ROImask);

AHITable1=array2table(AHI.Total,'VariableNames',rowNames,'RowNames',{'AHITotal'});
%                 displaytext=['Total AHI: ' num2str(AHI.Total(58),4) ' per hr'];
%                 disp(displaytext); set(handletext,'String',displaytext); drawnow;
%                 displaytext=['NREM supine AHI: ' num2str(AHI.Total(16),4) ' per hr'];
%                 disp(displaytext); set(handletext,'String',displaytext); drawnow;
%

%% Event subset, e.g. 4pc AHI
desatlist=[3.1 4 3]; %note 4percent will not include arousal

RO1mask2 = (ROImask==1 & SpO2Available==1)*1;


[AHI.ThreePA,Evts,~] = getDesatArSubset(Evts,EpochsXHz,Time,Position,RO1mask2,desatlist(1));

[AHI.FourP,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,RO1mask2,desatlist(2));

[AHI.ThreeP,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,RO1mask2,desatlist(3));

if 1
[AHI.ThreePAb,Evts,~] = getDesatArSubset(Evts,EpochsXHz,Time,Position,ROImask,desatlist(1));
[AHI.FourPb,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,ROImask,desatlist(2));
[AHI.ThreePb,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,ROImask,desatlist(3));
end

AHITable2=array2table(AHI.ThreePA,'VariableNames',rowNames,'RowNames',{'AHI3PA'});
AHITable3=array2table(AHI.FourP,'VariableNames',rowNames,'RowNames',{'AHI4P'});
AHITable4=array2table(AHI.ThreeP,'VariableNames',rowNames,'RowNames',{'AHI3P'});
AHIdata2=[AHITable1;AHITable2;AHITable3;AHITable4];

if 1
AHITable2b=array2table(AHI.ThreePAb,'VariableNames',rowNames,'RowNames',{'AHI3PAb'});
AHITable3b=array2table(AHI.FourPb,'VariableNames',rowNames,'RowNames',{'AHI4Pb'});
AHITable4b=array2table(AHI.ThreePb,'VariableNames',rowNames,'RowNames',{'AHI3Pb'});
AHIdata2b=[AHITable2b;AHITable3b;AHITable4b];
AHIdata2 = [AHIdata2;AHIdata2b];
end

AHItotal = AHIdata2.AllSleepAllPahi(1);
AHI3pa = AHIdata2.AllSleepAllPahi(2);
AHI4 = AHIdata2.AllSleepAllPahi(3);
AHI3 = AHIdata2.AllSleepAllPahi(4);
displaytext=['AHItotal : ' num2str(round(AHI3pa,1)) ' events/hr'];
disp(displaytext);
displaytext=['AHI3pa Per PUP: ' num2str(round(AHI3pa,1)) ' events/hr'];
disp(displaytext);
displaytext=['AHI4: ' num2str(round(AHI4,1)) ' events/hr'];
disp(displaytext);
displaytext=['AHI3: ' num2str(round(AHI3,1)) ' events/hr'];
disp(displaytext);

%% Edits
%EvtsData{1}=Evts;
Evts.AHIdata = AHIdata;
Evts.AHIdata2 = AHIdata2;


%% hypoxic burden
HBtotal = nanmean(Evts.RespT.HBarea)/60*AHIdata2.AllSleepAllPahi(1);
displaytext=['HBtotal: ' num2str(HBtotal,4) ' %min per hr'];
disp(displaytext);

Evts.HBtotal = HBtotal;