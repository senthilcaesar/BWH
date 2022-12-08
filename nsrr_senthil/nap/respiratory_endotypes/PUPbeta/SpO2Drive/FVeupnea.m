function Alpha = FVeupnea(SpO2,SpO2wake)
%uses S_O2toP_O2

fudgeS = SpO2wake-97.6; %use for patients who show SpO2eupnea > 98: e.g. 1723: Seupnea=97.6 (find stable wake meanSpO2 and subtract difference from 97.6)

%% FVeupnea_est
maxsat=98.828;
SaO2_=SpO2-fudgeS;
SaO2_(SaO2_>maxsat)=maxsat;
SaO2_(isnan(SaO2_))=[];
temp=S_O2toP_O2(1,SaO2_,23400,150);
PpO2_mean = nanmean(temp);
FVeupnea_est=1-(150-PpO2_mean)/50; %Veupnea = mean(VI_rs)*(1-FVeupneaF); FVeupnea=1-Veupnea/meanVE
Alpha = 1/(1-FVeupnea_est);

Alpha(Alpha>1.1)=1.1;

