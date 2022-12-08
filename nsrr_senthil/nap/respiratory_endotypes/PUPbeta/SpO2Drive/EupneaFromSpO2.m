%% EupneaFromSpO2

%Wake SpO2 from "BreathDataTable"
W=length(BreathDataTable{1});
SpO2WakeSyncData=[];
clear me
for w=1:W
    try
%         if criteria(w)==0 %pos and state etc are irrelevant here
%             continue
%         end
        TableTemp=BreathDataTable{1}{w};
        if sum(size(TableTemp)==[1 1])==2
            continue
        end
        spo2_w = TableTemp.spo2;
        hyp_w = TableTemp.hypnog_B;
        SpO2WakeSyncData = [SpO2WakeSyncData;[spo2_w hyp_w];[NaN NaN]];
    catch me
    end
end
spo2_w=SpO2WakeSyncData(:,1);
hyp_w=SpO2WakeSyncData(:,2);
criteria_ = hyp_w==4;

% figure(89)
% histogram(spo2_w)

temp = spo2_w(criteria_);
I = temp>prctile(temp,90)|temp<prctile(temp,10);
temp(I)=[];
spo2wakemean10to90 = nanmean(temp);
spo2wake75th = prctile(temp,75);

%% 
%then run FVeupnea when needed

