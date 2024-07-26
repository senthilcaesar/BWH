function [AHItemp,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,CPAPoff,desat)
try
HypCodes=[4 6 10 11]; % hypopnea codes
HypInd=any(Evts.RespT.EventCodes==HypCodes,2); % hypopneas

ApnInd=(Evts.RespT.EventCodes==2)|(Evts.RespT.EventCodes==3)|(Evts.RespT.EventCodes==5); %All Apneas

clear tempEvt ArDesatInd ApHypInd;

tempEvt=Evts;
if desat==3.1
    Criteria=((Evts.RespT.SpO2DeltaE>=3)| (Evts.RespT.ArE==1)) & (Evts.RespT.Epochs<=3); % arousals or 3 % desat and sleep
    newname = {'InclAHI3a'};
    % CPAPoff
elseif desat==3
    Criteria=(Evts.RespT.SpO2DeltaE>=3)  & (Evts.RespT.Epochs<=3); % 3% desat and sleep
    newname = {'InclAHI3'};
elseif desat==4
    Criteria=(Evts.RespT.SpO2DeltaE>=4)  & (Evts.RespT.Epochs<=3); % 4% desat and sleep
    newname = {'InclAHI4'};
end
ApHypInd=(ApnInd)|(HypInd & Criteria); % real apneas and hypopneas
tempEvt.RespT(~ ApHypInd,:)=[]; % delete those events without real apneas/hypopneas
[AHItemp,~,rowNames]=getAHIEvtSubset(tempEvt,EpochsXHz,Time,Position,CPAPoff);

%newname=cellstr([num2str('ApHypInd') '_' num2str(desat) 'pAr']);
Tnew=table(ApHypInd,'VariableNames',newname);
if any(Evts.RespT.Properties.VariableNames==string(newname)) %check if newname is already in Evts.RespT table
    Evts.RespT(:,newname)=Tnew; %replace
else
    Evts.RespT=[Evts.RespT Tnew]; %add column
end
catch
    disp('Error within getAHIAll-getDesatArSubset');
end

