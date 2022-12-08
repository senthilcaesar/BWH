%% PhenoDrive Phenotypes

%
global settings
settings.PlotEventData2Option=1;

cd 'G:\Dropbox (Personal)\PhenotypeDrive2018\Analyzed\EventAnalyzed'
addpath(genpath('G:\Dropbox (Personal)\PUPbeta_git\PUPbeta20190629\'));

%Import manual assessment
[~,~,raw]=xlsread('G:\Dropbox (Personal)\PhenotypeDrive2018\PUPstart\AMasterSpreadsheet','BG4:BG53');
raw = string(raw);
P = NaN*ones(50,1);
    P(raw=="A")=1;
    P(raw=="B")=0;

%Import AHI
[AHI,~,~]=xlsread('G:\Dropbox (Personal)\PhenotypeDrive2018\PUPstart\AMasterSpreadsheet','IU4:IU53');


    
Time = [-90:1:90]';
%%
R1 = NaN*ones(50,1);
FlowvsDriveNadirTimeDiff = NaN*ones(50,1);
FlowvsDriveNadirTimeDiffF = NaN*ones(50,1);
Vdrive_VENadir= NaN*ones(50,1);
    VE_VENadir= NaN*ones(50,1);
    Vdrive_start= NaN*ones(50,1);
    VE_start= NaN*ones(50,1);
    
%%
for i=29%1:50
    i 
    try
    load(['PhenoDrive2019_EvtData_' num2str(i)])
    
    PlotEventData2(EvtData);
    Ievent = EvtData.EnTime>=EvtData.deltaT2*1 & round(EvtData.EnTime)<=EvtData.deltaT2*0.05;
    %Ievent = EvtData.EnTime>=round(EvtData.deltaT2) & round(EvtData.EnTime)<=-2;
    
    h=figure(15);
    allaxes = findall(h, 'type', 'axes');
    set(h,'CurrentAxes',allaxes(3))

    I = Ievent;
    hold on
    Y = 100*nanmean(EvtData.VI(:,I)); 
    X = 100*nanmean(EvtData.VdriveEdiNorm(:,I));
    plot(EvtData.EnTime(I),Y,'r'); 
    plot(EvtData.EnTime(I),X,'r');
    
    [~,Xmini]=min(X);
    [~,Ymini]=min(Y);
    
    Vdrive_VENadir(i)=X(Xmini);
    VE_VENadir(i)=Y(Xmini);
    Vdrive_start(i)=X(1);
    VE_start(i)=Y(1);
    
    FlowvsDriveNadirTimeDiff(i) = (Ymini-Xmini);
    FlowvsDriveNadirTimeDiffF(i) = FlowvsDriveNadirTimeDiff(i)/abs(EvtData.deltaT2);
    
    R1backup
    R1(i)=corr(Y(:),X(:))
    
    figure(2);
    plot(X,Y,'k-.')
    
    catch me
        'fail'
    end
    %pause
end


    
sum(R1>0.5)/sum(R1>-Inf)

I = ~isnan(R1.*FlowvsDriveNadirTimeDiff);
corr(FlowvsDriveNadirTimeDiff(I),R1(I))
%%
figure(6); set(gcf,'color',[1 1 1]);
plot(FlowvsDriveNadirTimeDiffF(P==1),R1(P==1),'.','markersize',18);
box off
hold on
plot(FlowvsDriveNadirTimeDiffF(P==0),R1(P==0),'.','markersize',18);
%%
figure(7); set(gcf,'color',[1 1 1]);
plot(FlowvsDriveNadirTimeDiff(P==1),R1(P==1),'.','markersize',18);
box off
hold on
plot(FlowvsDriveNadirTimeDiff(P==0),R1(P==0),'.','markersize',18);

figure(8); set(gcf,'color',[1 1 1]);
plot(FlowvsDriveNadirTimeDiff(P==1),R1(P==1),'.','markersize',18);
box off
hold on
plot(FlowvsDriveNadirTimeDiff(P==0),R1(P==0),'.','markersize',18);


%%
mdltbl = table(R1,FlowvsDriveNadirTimeDiff,P);
mdl = fitglm(mdltbl,'P ~ R1 + FlowvsDriveNadirTimeDiff','distribution','binomial');
predP = predict(mdl,mdltbl);


%%
figure(8); set(gcf,'color',[1 1 1]);
plot(AHI(P==1),R1(P==1),'.','markersize',18);
box off
hold on
plot(AHI(P==0),R1(P==0),'.','markersize',18);

%% N analyzed
sum(P==1)
sum(P==0)
Iincl = P==1|P==0;
max(AHI(Iincl))
min(AHI(Iincl))
nanmedian(AHI(P==1))
nanmedian(AHI(P==0))
[~,p]=ttest2(AHI(P==1),AHI(P==0))

sum(P==1 & AHI>5)
sum(P==0 & AHI>5)

sum(P==1 & AHI>15)
sum(P==0 & AHI>15)

nanmedian(R1(Iincl))
prctile(R1(Iincl),[50 25 75])
prctile(FlowvsDriveNadirTimeDiff(Iincl),[50 25 75])
prctile(Vdrive_VENadir(Iincl),[50 25 75])
prctile(VE_VENadir(Iincl),[50 25 75])
prctile(Vdrive_start(Iincl),[50 25 75])
prctile(VE_start(Iincl),[50 25 75])
prctile(VdriveFall(Iincl),[50 25 75])
prctile(VEFall(Iincl),[50 25 75])
prctile(VdriveVEnadirdifference(Iincl),[50 25 75])

prctile(R1(P==1),[50 25 75])
prctile(FlowvsDriveNadirTimeDiff(P==1),[50 25 75])
prctile(Vdrive_VENadir(P==1),[50 25 75])
prctile(VE_VENadir(P==1),[50 25 75])
prctile(Vdrive_start(P==1),[50 25 75])
prctile(VE_start(P==1),[50 25 75])
prctile(VdriveFall(P==1),[50 25 75])
prctile(VEFall(P==1),[50 25 75])
prctile(VdriveVEnadirdifference(P==1),[50 25 75])

prctile(R1(P==0),[50 25 75])
prctile(FlowvsDriveNadirTimeDiff(P==0),[50 25 75])
prctile(Vdrive_VENadir(P==0),[50 25 75])
prctile(VE_VENadir(P==0),[50 25 75])
prctile(Vdrive_start(P==0),[50 25 75])
prctile(VE_start(P==0),[50 25 75])
prctile(VdriveFall(P==0),[50 25 75])
prctile(VEFall(P==0),[50 25 75])
prctile(VdriveVEnadirdifference(P==0),[50 25 75])


ranksum(R1(P==0),R1(P==1))
ranksum(FlowvsDriveNadirTimeDiff(P==0),FlowvsDriveNadirTimeDiff(P==1))
ranksum(Vdrive_start(P==0),Vdrive_start(P==1))
ranksum(VdriveFall(P==0),VdriveFall(P==1))
ranksum(Vdrive_VENadir(P==0),Vdrive_VENadir(P==1))
ranksum(VdriveFall(P==0),VdriveFall(P==1))

VdriveFall = Vdrive_start - Vdrive_VENadir;
VEFall = VE_start - VE_VENadir;
VdriveVEnadirdifference = Vdrive_VENadir - VE_VENadir;

Vdrive_VENadir(i)=X(Xmini);
VE_VENadir(i)=Y(Xmini);
    Vdrive_start(i)=X(1);
    VE_start(i)=Y(1);
    
    
    
