function [Pesdelta2,Pmusdelta2,Ipes,VIpes,Pes_baselineest,Pespeak]= PesAnalysis(Pes,Flow,Time,I,Ecw,PesArtifact)
global settings
%find Pcw

if ~exist('PesArtifact')
    PesArtifact=1
end

% options
ShowFigures = 1; % settings.plotfigure is unknown at the level this is run

dt = Time(2)-Time(1);
useactualzero=0;
filter_HFcutoff_butter0 = 2;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
Pesf = filtfilt(B_butter0,A_butter0,Pes);

voltemp1 = cumsum(Flow)*dt;
vollower = interp1(Time(I.starti),voltemp1(I.starti),Time,'linear');
vollower(1:I.starti(1)-1)=voltemp1(I.starti(1)); %extrapolation using nearest
vollower(I.starti(end):end)=voltemp1(I.starti(end)); %extrapolation using nearest
voltemp2 = voltemp1 - vollower;
Pcw = Ecw*voltemp2;

% peslower = interp1(Time(I.starti),Pesf(I.starti),Time,'linear');
%     peslower(1:I.starti(1)-1)=Pesf(I.starti(1)); %extrapolation using nearest
%     peslower(I.starti(end):end)=Pesf(I.starti(end)); %extrapolation using nearest
%
%     temp = Pesf(I.starti); Pes_baselineest = median(temp);
%     figure(100); plot(Time,Pesf,Time,peslower)


%%
pcrtiles = 0:1:100;
PesPrctiles = prctile(Pesf,pcrtiles);
IX=find(Pesf(1)<=PesPrctiles,1);
if 0
    [Vdot_intended,parameters,rsquared]=PesToVflowFit(Flow,Pesf,Time,5);
else
    parameters = [-15 -10 0 0 100-pcrtiles(IX)];%[20 10 0 0 100-pcrtiles(IX)];
end

[Vdot_intended]=PesToVflowRun(Flow,Pesf,Time,Inf,parameters);
% [Ipes,~,~,~,~,~,~,~]=Vflowanalysis2CMM(Vdot_intended,Time,dt,1,useactualzero,1);
settingsbackup=settings; 


try
settings.sqrt_scaling = 0;
settings.exponent = 1;
settings.modBB_i_start=0;
%[~,~,BB_i_start,BB_i_mid,BB_i_end,~,~,~,Ttotpes,~,~,VTpes] = ...  % old
[~,~,BB_i_start,BB_i_mid,BB_i_end,~,~,~,Ttotpes,~,~,VTpes] = ... %new
    VEfromFlow(Time,Vdot_intended); %VEfromFlow_sqrt_V16(Time,Vdot_intended);

settings = settingsbackup;
catch me
settings = settingsbackup;    
end

Ipes.starti = BB_i_start;
Ipes.midi = BB_i_mid;
Ipes.endi = BB_i_end;

VIpes = VTpes./Ttotpes;
Pes_baselineest = median(Pesf(Ipes.starti));
Pmus = Pes-Pcw-Pes_baselineest;


[Vdot_intended]=PesToVflowRun(Flow,Pmus,Time,Inf,parameters);
% [Ipes,~,~,~,~,~,~,~]=Vflowanalysis2CMM(Vdot_intended,Time,dt,1,useactualzero,1);
settingsbackup=settings;


if ShowFigures
    figure(42) %Pes breath/baseline detection using low pass filtered trace
    set(gcf,'color',[1 1 1]);
    ax42(1)=subplot(3,1,1); plot(Time,Flow); hold('on');plot(Time,Vdot_intended);
    plot(Time(Ipes.starti-1),Flow(Ipes.starti-1),'r.');
    plot(Time(Ipes.midi),Flow(Ipes.midi),'k.'); box('off'); hold('off');
    ylabel('Flow');
    ax42(2)=subplot(3,1,2); plot(Time,Vdot_intended); hold('on');
    plot(Time(Ipes.starti-1),Vdot_intended(Ipes.starti-1),'r.');
    plot(Time(Ipes.midi),Vdot_intended(Ipes.midi),'k.'); box('off');  hold('off');
    ylabel('Drive->Flow');
    ax42(3)=subplot(3,1,3); plot(Time,Pesf); hold('on');
    plot(Time(Ipes.starti-1),Pesf(Ipes.starti-1),'r.');
    plot(Time(Ipes.midi),Pesf(Ipes.midi),'k.'); box('off');  hold('off');
    ylabel('Drive');
    linkaxes(ax42,'x');
end


%find baseline Pes
%figure(456);
%plot(Time,Pes,Time,0*Time+Pes_baselineest);

   
clear minPesB2 maxPesB2 Pesdelta2 Pmusdelta2 Pespeak
for i=1:length(Ipes.starti)
    minPesB2(i,1) = min(Pes(Ipes.starti(i):Ipes.endi(i))-Pes_baselineest);
    maxPesB2(i,1) = max(Pes(Ipes.starti(i):Ipes.endi(i))-Pes_baselineest);
    Pesdelta2(i,1) = Pes(Ipes.starti(i))-min(Pes(Ipes.starti(i):Ipes.endi(i)));
    Pmusdelta2(i,1) = Pmus(Ipes.starti(i))-min(Pmus(Ipes.starti(i):Ipes.endi(i)));    
    Pespeak(i,1) = min(Pes(Ipes.starti(i):Ipes.endi(i)));
end
%% Pes artifact removal
if PesArtifact == 1
    
PesUpperLim = 0.5*abs(median(minPesB2));
temp2 = maxPesB2>PesUpperLim; % changed from 0.5*abs(minPesB2); %test for artifact: upwards swing from "baseline" is larger than the downwards swing from the "baseline"
Pesdelta2(temp2)=NaN;
Pmusdelta2(temp2)=NaN;
Pespeak(temp2)=NaN;
VIpes(temp2)=NaN;
end

%%
if ShowFigures
    ax42(3)=subplot(3,1,3); plot(Time,Pes-Pes_baselineest); hold('on');
    plot(Time(Ipes.starti-1),Pes(Ipes.starti-1)-Pes_baselineest,'r.');
    plot(Time(Ipes.midi),Pes(Ipes.midi)-Pes_baselineest,'k.');
    plot(Time,Pcw,'-','color',[0.5 0.5 0.5]); box('off');
    hold('off');
    ylabel('Pes');
    linkaxes(ax42,'x');
end
