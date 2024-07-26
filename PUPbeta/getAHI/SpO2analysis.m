function [SpO2data,SaO2XHz]=SpO2analysis(SigT,ChannelsList,CPAPoff)

% Added ChannelsList as input to function---needed for calling it in spo2
% analysis in ConvertToMat --9/28/2020 RMA

global settings
%%
showfigures=settings.plotfigure;

Time = SigT.Time;
dt = 1/settings.Fs;
EpochsXHz = SigT.Epochs;
if ~isnan(find(strcmp(ChannelsList,'Position')==1))
    Position = round(SigT.Position);
else
    Position = NaN * SigT{:,1};
end
SaO2XHz =SigT.SpO2;


Supine = sum((Position==settings.supinepositioncode),2)==1;

%NREM: (sum((EpochsXHz==NREMcodes),2)==1);
NREMcodes=[0 1 2];
REMcode=3;
Sleepcodes=[0 1 2 3];
Wakecode=4;

if showfigures
    figure(101);
    ax(2)=subplot(4,1,2);
    hold('off');
    plot(Time,SaO2XHz,'k'); hold('on');
    set(gca,'FontName','Arial Narrow','FontSize',10,'box','off');
    ylabel('SpO2 %');
end

SaO2XHz3=SpO2ArtifactReject(SaO2XHz,dt);

if showfigures
    plot(Time,SaO2XHz3,'r');
end

% This part is quick:
%uses sleep as a criteria (considers delay) [i.e. current or last epoch must be sleep]
lastepochsXHz=[NaN*zeros(round(30/dt),1);EpochsXHz((1:(length(EpochsXHz)-round(30/dt))))];
[nadir_desat_sleep,tempi] = min(SaO2XHz3(((sum((EpochsXHz==Sleepcodes),2)==1)|(sum((lastepochsXHz==Sleepcodes),2)==1))&~isnan(SaO2XHz3)&CPAPoff));
temptime = Time(((sum((EpochsXHz==Sleepcodes),2)==1)|(sum((lastepochsXHz==Sleepcodes),2)==1))&~isnan(SaO2XHz3)&CPAPoff);
nadir_desat_sleep_time=temptime(tempi);

clear tempi temptime
% Tidy up desat
SaO2XHz=SaO2XHz3;
clear SaO2XHz1 SpO2neighborhoodD_dt SpO2neighborhoodD_i SaO2XHz2 SaO2info1 SaO2info2 proximity_score
clear I I1 I2 I3;
%Plot

if showfigures
    plot(nadir_desat_sleep_time,nadir_desat_sleep,'k.');  %,Time,SaO2XHz3,
    plot(Time,100*(1-CPAPoff),'g');
    plot(Time,100+5*EpochsXHz,'b');
end

%% ODI
Include(:,1) = ~isnan(SaO2XHz)&EpochsXHz>=0&EpochsXHz<4&CPAPoff; %All Sleep (0,1,2,3)
Include(:,2) = ~isnan(SaO2XHz)&EpochsXHz>=0&EpochsXHz<3&CPAPoff; %NREM (0,1,2)
Include(:,3) = ~isnan(SaO2XHz)&EpochsXHz==3&CPAPoff; %REM (3)

[ODI3Array,ODI4Array] = CalcODI(SaO2XHz,dt,Include);

Include(:,4) = ~isnan(SaO2XHz)&EpochsXHz==4; %Wake

%% Desat info
SpO2mean_wake=nanmean(SaO2XHz(EpochsXHz==Wakecode));
nadir_desat_sleep_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1),[0 1 5 25 50 75 95 99 100]);
SpO2mean_sleep=nanmean(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
SpO2mean_sup=nanmean(SaO2XHz(Supine&CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
SpO2below90p_sleep_prct = 100*length(SaO2XHz(SaO2XHz<90&CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1))/length(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
SpO2below90p_sleep_prct_NREMsup = 100*length(SaO2XHz(SaO2XHz<90&Supine&CPAPoff&sum((EpochsXHz==NREMcodes),2)==1))/length(SaO2XHz(Supine&CPAPoff&sum((EpochsXHz==NREMcodes),2)==1));
SpO2below95p_sleep_prct = 100*length(SaO2XHz(SaO2XHz<95&CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1))/length(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
nadir_desat_sup_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&EpochsXHz>0&(Supine)),[0 1 5 25 50 75 95 99 100]);
nadir_desat_NREMsup_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&sum((EpochsXHz==NREMcodes),2)==1&(Supine)),[0 1 5 25 50 75 95 99 100]);
nadir_desat_REMsup_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&(EpochsXHz==REMcode)&(Supine)),[0 1 5 25 50 75 95 99 100]);
nadir_desat_NREMlat_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&(sum((EpochsXHz==NREMcodes),2)==1)&~(Supine)),[0 1 5 25 50 75 95 99 100]);
SpO2mean_NREMsup=mean(SaO2XHz(CPAPoff&(sum((EpochsXHz==NREMcodes),2)==1)&(Supine)&~isnan(SaO2XHz)));
SpO2mean_REMsup=mean(SaO2XHz(CPAPoff&(EpochsXHz==REMcode)&(Supine)&~isnan(SaO2XHz)));
SpO2mean_NREMlat=mean(SaO2XHz(CPAPoff&(sum((EpochsXHz==NREMcodes),2)==1)&~(Supine)&~isnan(SaO2XHz)));

% clear SpO2data
% SpO2data = [...
%     SpO2mean_wake ...
%     nadir_desat_sleep ...
%     SpO2below90p_sleep_prct ...
%     SpO2mean_sleep ...
%     SpO2mean_sup ...
%     SpO2mean_NREMsup ...
%     SpO2mean_REMsup ...
%     SpO2mean_NREMlat ...
%     SpO2below90p_sleep_prct_NREMsup ...
%     nadir_desat_NREMsup_prct_0_1_5_25_50_75_95_99 ...
%     ]';
%%

% Baseline SpO2
% "Min swing" <2
% Before first sleep onset


ODI3=ODI3Array(1);
ODI3nrem=ODI3Array(2);
ODI3rem=ODI3Array(3);
ODI4=ODI4Array(1);
ODI4nrem=ODI4Array(2);
ODI4rem=ODI4Array(3);
SpO2mean=nanmean(SaO2XHz(Include(:,1)));
SpO2meannrem=nanmean(SaO2XHz(Include(:,2)));
SpO2meanrem=nanmean(SaO2XHz(Include(:,3)));
SpO2meanwake=nanmean(SaO2XHz(Include(:,4)));
%SpO2medianwake=nanmedian(SaO2XHz(Include(:,4)));
SpO2nadir=prctile(SaO2XHz(Include(:,1)),0);
SpO2nadirnrem=prctile(SaO2XHz(Include(:,2)),0);
SpO2nadirrem=prctile(SaO2XHz(Include(:,3)),0);

SpO2data = table(ODI3,ODI3nrem,ODI3rem,ODI4,ODI4nrem,ODI4rem,SpO2mean,SpO2meannrem,SpO2meanrem,SpO2meanwake,SpO2nadir,SpO2nadirnrem,SpO2nadirrem,SpO2below90p_sleep_prct,SpO2below95p_sleep_prct);
