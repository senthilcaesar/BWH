%% setting up modelled V chem
LG0 = 3;
tau1 = 20;
delay = 15;
tau2 = tau1;

%LG0 = LG1*(1+(2*pi*tau1)^2)^0.5;
LG1 = LG0/(1+(2*pi*tau1/60)^2)^0.5 / (1+(2*pi*tau2/60)^2)^0.5;

%LG1 = abs(LG0 * 1/(1 + j*2*pi/60*tau1) * 1/(1 + j*2*pi/60*tau2))


dt = 1;
Time = (-60:dt:210)';
Time=(1:dt:201)';

alpha1=(tau1./dt)./(1+tau1./dt); % Forward Euler approximation of alph
alpha2=1-dt/tau1; % Backward Euler approximation of alpha
beta1=-LG0./(1+tau1./dt); % Forward Euler approximation of beta
beta2=-LG0./(tau1./dt); % Backward Euler approximation of beta
alpha=(alpha1+alpha2)/2; % alpha
beta=(beta1+beta2)/2; % beta

alpha1X=(tau2./dt)./(1+tau2./dt); % Forward Euler approximation of alph
alpha2X=1-dt/tau2; % Backward Euler approximation of alpha
beta1X=1./(1+tau2./dt); % Forward Euler approximation of beta
beta2X=1./(tau2./dt); % Backward Euler approximation of beta
alphaX=(alpha1X+alpha2X)/2; % alpha
betaX=(beta1X+beta2X)/2; % beta

% alpha=1-dt/tau1; % Backward Euler approximation of alpha
% alphaX=1-dt/tau2; % Backward Euler approximation of alpha
% beta=-k./(tau1./dt); % Backward Euler approximation of beta
% betaX=1./(tau2./dt); % Backward Euler approximation of beta

VE = 0*Time;
E1 = 0*Time + 1;
E1(Time>30 & Time<60) = 0;
E1(Time>120 & Time<150) = 0;
VE(E1==0)=-1;

delayedVE = 0*Time;
delayi = round(delay/dt);
delayedVE(delayi+1:end) = VE(1:end-delayi);

Vchem = 0*Time;
Vchem2 = 0*Time;
for k=(delayi+1):length(Time)
    delayedVE(k) = VE(k-delayi);
    Vchem(k) = alpha*Vchem(k-1) + delayedVE(k)*beta;
    Vchem2(k) = alphaX*Vchem2(k-1) + Vchem(k)*betaX;
    %VE(k) = Vchem2(k);
    VE(k) = Vchem2(k);
    if E1(k)==0
        VE(k) = -1;
    end
end
Vchem = Vchem2;

%Vout = Vchem2;
%Vout(E1==0)=VE(E1==0);

figure(1)
plot(Time,[VE Vchem]);


%%
%             - Column 1: Breath-by-Breath Minute Ventilation -- Must already be "mean subtracted" i.e. eupnea=0.
%             - Column 2: Obstructive respiratory event series. Elements=1 for an un-obstructed breath; Elements=0 for an obstructed breath.
%             - Column 3: Arousal Scoring series: defines whether there is an arousal present in each breath (arousal present, AR[n]=1; No arousal present, AR[n]=0);
%             - Column 4: The time that each breath in the series occurs.
%             - Column 5: Scored Central Apnea array series: defines whether each breath is scored as a central apnea. elements=1 for breaths in a scored central apnea; Elements=0 for breaths not classified as a central apnea.
%             - Column 6: T_tot series (currently unused)
%             - Column 7: T_tot series for the previous breath.
% global settings
% settings.plotfigureLGparameters = 1

%OLD - plots event by event within subject
M=10
load(['EventAnalyzed\PhenoDrive2020_' num2str(M)])

Veupnea=1;
VdriveExponent=1; %0.67 perhaps for those with Vdrive>VE at arousal?
settings.tau2ontau1=1;
settings.VchemNonlinearTransform=1;
minareffect=1;


for i=3:size(Boxes.VI,2)
VE=Boxes.VI(i,:)'-mean(Boxes.VI(i,:));
Vchem=Boxes.VdriveEdiNorm(i,:)'-mean(Boxes.VdriveEdiNorm(i,:));
E1=VE*0;
Time=(1:dt:201)';

%ArSig = Boxes.ArPr(i,:)*1; %EventsAr
ArSig = 0*VE;
ArSig = ArSignalEffects(1*(Boxes.ArPr(i,:)'*1 > 0.75),2,20);
ArSig = (Boxes.ArPr(i,:)'-minareffect)/(1-minareffect); ArSig(ArSig<0)=0;

%Data = [VE E1 E1*0 Time E1*0 E1*0+dt E1*0+dt Vchem];

Data = [VE E1 ArSig Time E1*0 E1*0+dt E1*0+dt Vchem];
VraOn=0;
Veupnea=1;
polyfitorder=0;
% 
% I = [-10+101+round(0*EAinfo.deltaT2):round(101+45)]';
% I(I<1)=[];
% Data = Data(I,:);
I=[1:201];

[Error,VchemOut,Varousal,LoopGainParameters,BestParameters,BestSSres,FitQuality,i_1,i_end,lowerSEM,upperSEM,CI_parameters] = ...
    FindBestModelParametersCItau2(Data,VraOn,Veupnea,polyfitorder)

%

% ParametersArray(i,:)=BestParameters;
% figure(2)
% plot(Time,[Vchem VchemOut]);


figure(2); clf(2);
plot(Time(I),1+[VE(I) Vchem(I)]);
hold on
plot(Time(I),1+[VchemOut]);
plot(Time(I),1+[Varousal+VchemOut],'k');
plot(Time(I),2+[Ensembles.ArPr(I) ArSig(I)]);
%legend on
%hold on

pause
end

%% MAIN FUNCTION FOR LOOP GAIN FITS TO EDI
%Run this section
%idcs   = strfind(dir,'\');
 newdir = settings.workdir
 if ~(exist([newdir '\CPAPTraits'], 'dir') == 7)
        mkdir([newdir '\CPAPTraits']);
 end
    

for M=1:60; %15 41 17* 4 46 19 14 23
    
dt = 1;
    try
load(['EventAnalyzed\PhenoDrive2020_' num2str(M)])

Veupnea=1;
VdriveExponent=1; %0.67 perhaps for those with Vdrive>VE at arousal?
settings.tau2ontau1=1;
settings.VchemNonlinearTransform=1;
minareffect=1;

VE=Ensembles.VI - Veupnea; %affects fit
Vchem=Ensembles.VdriveEdiNorm.^VdriveExponent - Veupnea; %doesn't affect fit (shift/scale)
%Vchem=Ensembles.VdrivePesNorm.^VdriveExponent - Veupnea; %doesn't affect fit (shift/scale)
Time=(1:dt:201)';

if 0
I = find(~isnan(sum(Boxes.VdriveEdiNorm,2)));
VE=Boxes.VI(I,:)' - 1;
Vchem=Boxes.VdriveEdiNorm(I,:)' - 1;
VE=VE(:);
Vchem=Vchem(:);
end

ArSig = Ensembles.ArPr*1; %EventsAr
ArSig = 0*VE;
ArSig = ArSignalEffects(1*(Ensembles.ArPr*1 > 0.75),2,20);
try
GGSig = Ensembles.GGpeak
catch
end
ArSig = (Ensembles.ArPr-minareffect)/(1-minareffect); ArSig(ArSig<0)=0;
   %ArSig=ArSig.^0.5 
%ArSig = ArSignalEffects(1*(Ensembles.ArPr*1 > 0.5),0,[]);
%ArSig = ArSignalEffects(1*(Ensembles.ArPr*1 > 0.5),1,20)

% ArSig = 1*(Ensembles.ArPr*1 > 0.5);
% AR1 = ArSig;

%ArSig = 1*(Ensembles.ArPr*1 > 0.5);
%ArSig = Ensembles.ArPrLogit*1; %EventsAr
%Time(:,1)=(1:dt:length(VE));

%ArSig = Ensembles.EventsAr*1; %EventsAr
E1=VE*0 +1;
%E1(Ensembles.ArPr*1 >0.5)=0;
%E1(Ensembles.EventsResp<0.5)=0;

Data = [VE E1 ArSig Time E1*0 E1*0+dt E1*0+dt Vchem];

I = [-10+101+round(0*EAinfo.deltaT2):round(101+45)]';

I(I<1)=[];
Data = Data(I,:);



VraOn=1;
Veupnea=1;
polyfitorder=0;

[Error,VchemOut,Varousal,LoopGainParameters,BestParameters,BestSSres,FitQuality,i_1,i_end,lowerSEM,upperSEM,CI_parameters] = ...
    FindBestModelParametersCItau2(Data,VraOn,Veupnea,polyfitorder)
%LoopGainParameters=[LG0 tau delay LGn Tn LG1 LG2 gamma];

%LoopGainParametersTemp = LoopGainParameters;
disp(['LG1=' num2str(LoopGainParameters(6))]);
disp(['VRA=' num2str(LoopGainParameters(8))]);

figure(2); clf(2);
plot(Time(I)-100,1+[VE(I) Vchem(I)]);
hold on
plot(Time(I)-100,1+[VchemOut]);
plot(Time(I)-100,1+[Varousal+VchemOut],'k');
plot(Time(I)-100,2+[Ensembles.ArPr(I)]) % ArSig(I)]);
try
plot(Time(I)-100,1+[(Ensembles.GGpeak(I)-mean(Ensembles.GGpeak(I)))/mean(Ensembles.GGpeak(I))]) % ArSig(I)]);

plot(Time(I)-100,1+[(Ensembles.DeltaPes(I)-mean(Ensembles.DeltaPes(I)))/mean(Ensembles.DeltaPes(I))]) % ArSig(I)]);
catch
end
xlabel('Time')
ylabel('VE  Vdrive  VChemOut   ArPr')

title(['LG1='  num2str(LoopGainParameters(6)) '; ' 'VRA=' num2str(LoopGainParameters(8))] )

set(gca,'box','off')
%legend on
%hold on
 pause

LGvalues(M,:)=LoopGainParameters
if 0
saveas(2,[newdir '\CPAPTraits' '\'  num2str(M) '_LG'],'fig');
end
%Save

    catch me
    end
end

%% Save all LG values
if 0
save([newdir '\CPAPTraits' '\'  'LGvalues'],'LGvalues');
end
%%
if 1
    
ParametersFixed = [NaN NaN NaN NaN];
Parameters = [LG1*1.2 tau1 0 0]
[SSres,Error,Rsq,VchemOut,Varousal,Penalized_Error] = TheModelX1tau2(Parameters,Data,i_1,delayedVE,polyfitorder,ParametersFixed)

figure(3)
plot(Time,[VE VchemOut]);
hold on
end

%% 
disp('All working!!')


