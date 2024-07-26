function SpO2data=SpO2analysis_DV(DataEventHypnog_Mat,CPAPoff)
global settings ChannelsList
%%
showfigures=0;

Time = DataEventHypnog_Mat(:,1);
dt = 1/settings.Fs;
EpochsXHz = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Epochs')==1));
if ~isnan(find(strcmp(ChannelsList,'Position')==1))
    Position = round(DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Position')==1)));
else 
    Position = NaN * DataEventHypnog_Mat(:,1);
end
SaO2XHz = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1));


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


%% Desat info
SpO2mean_wake=nanmean(SaO2XHz(EpochsXHz==Wakecode));
nadir_desat_sleep_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1),[0 1 5 25 50 75 95 99 100]);
SpO2mean_sleep=nanmean(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
SpO2mean_sup=nanmean(SaO2XHz(Supine&CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
SpO2below90p_sleep_prct = 100*length(SaO2XHz(SaO2XHz<90&CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1))/length(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
SpO2below90p_sleep_prct_NREMsup = 100*length(SaO2XHz(SaO2XHz<90&Supine&CPAPoff&sum((EpochsXHz==NREMcodes),2)==1))/length(SaO2XHz(Supine&CPAPoff&sum((EpochsXHz==NREMcodes),2)==1));
SpO2below90p_sleep_prct = 100*length(SaO2XHz(SaO2XHz<95&CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1))/length(SaO2XHz(CPAPoff&sum((EpochsXHz==Sleepcodes),2)==1));
nadir_desat_sup_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&EpochsXHz>0&(Supine)),[0 1 5 25 50 75 95 99 100]);
nadir_desat_NREMsup_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&sum((EpochsXHz==NREMcodes),2)==1&(Supine)),[0 1 5 25 50 75 95 99 100]);
nadir_desat_REMsup_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&(EpochsXHz==REMcode)&(Supine)),[0 1 5 25 50 75 95 99 100]);
nadir_desat_NREMlat_prct_0_1_5_25_50_75_95_99 = prctile(SaO2XHz(CPAPoff&(sum((EpochsXHz==NREMcodes),2)==1)&~(Supine)),[0 1 5 25 50 75 95 99 100]);
SpO2mean_NREMsup=mean(SaO2XHz(CPAPoff&(sum((EpochsXHz==NREMcodes),2)==1)&(Supine)&~isnan(SaO2XHz)));
SpO2mean_REMsup=mean(SaO2XHz(CPAPoff&(EpochsXHz==REMcode)&(Supine)&~isnan(SaO2XHz)));
SpO2mean_NREMlat=mean(SaO2XHz(CPAPoff&(sum((EpochsXHz==NREMcodes),2)==1)&~(Supine)&~isnan(SaO2XHz)));

clear SpO2data
SpO2data = [...
    SpO2mean_wake ...
    nadir_desat_sleep ...
    SpO2below90p_sleep_prct ...
    SpO2mean_sleep ...
    SpO2mean_sup ...
    SpO2mean_NREMsup ...
    SpO2mean_REMsup ...
    SpO2mean_NREMlat ...
    SpO2below90p_sleep_prct_NREMsup ...
    nadir_desat_NREMsup_prct_0_1_5_25_50_75_95_99 ...
    ]';

% SpO2table = table(...
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
%     )';



