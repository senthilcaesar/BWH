try
Sup = load('C:\Users\ludov\Dropbox\PhenotypeDrive2018\REMvNREM\SupineT1');
Lat = load('C:\Users\ludov\Dropbox\PhenotypeDrive2018\REMvNREM\LateralT1');
All = load('C:\Users\ludov\Dropbox\PhenotypeDrive2018\REMvNREM\AllT1');

catch
All =     load('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\REMvNREM\AllT1.mat')
Lat = load('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\REMvNREM\LateralT1.mat')
Sup = load('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\REMvNREM\SupineT1.mat')

end

Lat.T1.Lateral = ones(height(Lat.T1),1);
Lat.T1.Supine = zeros(height(Lat.T1),1);

Sup.T1.Lateral = zeros(height(Lat.T1),1);
Sup.T1.Supine = ones(height(Lat.T1),1);

T1both = [Sup.T1; Lat.T1];

%how many samples
NdatapointsSupine = sum(~isnan(Sup.T1.Flow))
NdatapointsLateral = sum(~isnan(Lat.T1.Flow))

mdl = fitglme(Sup.T1,'Flow ~ REM  + Drive^2 + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
% mdl = fitglme(Sup.T1,'Flow ~ Drive*REM  + Drive^2 + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve

mdl = fitglme(Lat.T1,'Flow ~ REM  + Drive^2 + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve

mdl = fitglme(T1both,'Flow ~ REM  + Drive^2 + Lateral + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
mdl = fitglme(T1both,'Flow ~ REM  + Drive^2 + Supine + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
T1both.nREM=1-T1both.REM;
%show this result
% the influence of REM on Flow (i.e. collapsibility / baseline activity) is not different supine v. lateral
% i.e. no effect whether lateral or supine
mdl = fitglme(T1both,'Flow ~ Drive^2 + Supine*REM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
%lateral -1.3317 [-9.0623     6.3989], p=0.73522
% mdl = fitglme(T1both,'Flow ~ Drive^2 + Lateral*REM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
%supine: -0.72571 [-3.9441     2.4927], p=0.65802
%difference (supine v lateral per Supine*REM interaction):
%0.60601 [-7.5819     8.7939], p=0.88447;
mdl = fitglme(T1both,'Flow ~ Drive^2 + Supine*nREM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve


%exaplin that adding the REM*Drive interaction did not change findings
mdlX = fitglme(T1both,'Flow ~ Drive*REM + Drive^2 + Supine*REM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve

%just to show that supine lowers ventilation in REM and nREM

mdlX = fitglme(T1both,'Flow ~ Drive*nREM + Drive^2 + Supine*nREM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve


%explain that Supine*Drive*REM was also not observed 
% the influence of REM on slope (responsivness) is not different supine v. lateral
mdl = fitglme(T1both,'Flow ~ Supine*Drive*REM + Drive^2 + Supine*REM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
mdl = fitglme(T1both,'Flow ~ Lateral*Drive*REM + Drive^2 + Lateral*REM + (1 + Drive | Subj)')  %REM does not change baseline or reduce the slope of the curve
%model analysis estimates that the slope(responsiveness) is larger in REM
%vs nREM both in supine (+19%) and lateral (+19%).

%%
InputT = All.T1;
%modelstring = 'Flow ~ REM*Drive  + Drive^2 + (1 + Drive | Subj)' %D
modelstring = 'Flow ~ REM + Drive^2 + (1 + Drive | Subj)' %C
mdlMain = fitglme(InputT,modelstring)  %REM does not change baseline or reduce the slope of the curve

mdlD = fitglme(InputT,'Drive ~ REM*Decile + Decile^2 + (1 | Subj)')  %REM does not change baseline or reduce the slope of the curve

tempT = table;
tempT.Decile = [0:9]';
tempT.REM = 0*tempT.Decile;
tempT.Subj = 0*tempT.Decile;
tempT.Drive = predict(mdlD,tempT);
[tempT.Flow,tempCI] = predict(mdlMain,tempT);
tempT.FlowCIl = tempCI(:,1);
tempT.FlowCIu = tempCI(:,2);
TnREM=tempT;

tempT = table;
tempT.Decile = [0:9]';
tempT.REM = 1+0*tempT.Decile;
tempT.Subj = 0*tempT.Decile;
tempT.Drive = predict(mdlD,tempT);
[tempT.Flow,tempCI] = predict(mdlMain,tempT);
tempT.FlowCIl = tempCI(:,1);
tempT.FlowCIu = tempCI(:,2);
TREM=tempT;

TbothX = [TnREM;TREM]

figure(2); clf(2);

plot(100*(TREM.Drive+1),TREM.Flow); hold on
plot(100*(TREM.Drive+1),TREM.FlowCIl)
plot(100*(TREM.Drive+1),TREM.FlowCIu)

fill([(100*(TREM.Drive+1));flipud(100*(TREM.Drive+1))],[TREM.FlowCIu;flipud(TREM.FlowCIl)],color1,'edgecolor','none','facealpha',0.5);
plot((100*(TREM.Drive+1)),TREM.Flow,'k','linewidth',1);


fill([(100*(TnREM.Drive+1));flipud(100*(TnREM.Drive+1))],[TnREM.FlowCIu;flipud(TnREM.FlowCIl)],color2,'edgecolor','none','facealpha',0.5);
plot((100*(TnREM.Drive+1)),TnREM.Flow,'k','linewidth',1);

%% Plot difference REM v nREM flow, in ALL Pos

tempT1 = InputT;
DriveRange=[(min(TnREM.Drive):0.05:max(TREM.Drive)) max(TREM.Drive)]+1;
clear deltaFlow nREMFlow
for i=1:length(DriveRange)
tempT1.Drive = InputT.Drive+1-DriveRange(i);
mdl1 = fitglme(tempT1,modelstring);  %REM does not change baseline or reduce the slope of the curve
Irow = find(mdl1.Coefficients.Name=="REM");
deltaFlow(i,:) = [mdl1.Coefficients.Estimate(Irow) mdl1.Coefficients.Lower(Irow) mdl1.Coefficients.Upper(Irow)];
nREMFlow(i)=mdl1.Coefficients.Estimate(1);
end


figure(3); clf(3);
plot((100*(TnREM.Drive+1)),TnREM.Flow, 'linewidth',1); hold on
%plot(100*DriveRange,nREMFlow(:)); hold on
%plot(100*(TREM.Drive+1),TREM.Flow); 


plot(100*DriveRange,nREMFlow(:)+deltaFlow(:,1));
plot(100*DriveRange,nREMFlow(:)+deltaFlow(:,2));
plot(100*DriveRange,nREMFlow(:)+deltaFlow(:,3));

deltaflowUpper = nREMFlow(:)+deltaFlow(:,3)
deltaflowLower = nREMFlow(:)+deltaFlow(:,2)
DR = 100*DriveRange'

fill([DR;flipud(DR)],[deltaflowUpper;flipud(deltaflowLower)],color1,'edgecolor','none','facealpha',0.5);
plot([(100*(TREM.Drive+1)),TREM.Flow], color1,'linewidth',1);