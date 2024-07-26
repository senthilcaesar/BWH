function [I,Vdot,VT,Ti,Te,VTi,VTe]=Vflowanalysis_knownI(Vflow,time,dt,minimum_figs,I)

%% Respiratory Trace Analyse

%transpose if needed:
if size(Vflow,2)<size(Vflow,1)
    Vflow=Vflow';
end
%transpose if needed:
if size(time,2)<size(time,1)
    time=time';
end

vol1=cumsum(Vflow)*dt;
leak = mean(Vflow);
voldetrend1 = leak*(time-time(1));
vol2 = vol1-voldetrend1;

if ~minimum_figs
    figure(); 
    ax1(1)=subplot(3,1,1); plot(time,vol1,time,voldetrend1);
    ax1(2)=subplot(3,1,2); plot(time,vol2);
end

%Filter volume gently for use in analysis
filter_HFcutoff_butter1 = 10;
filter_LFcutoff_butter1 = 1/15;
filter_order1 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));
vol2_filtered1 = filtfilt(B_butterHcut,A_butterHcut,vol2);
Vflow_filtered1 = filtfilt(B_butterHcut,A_butterHcut,Vflow);

if ~minimum_figs
     
    ax1(2)=subplot(3,1,2); plot(time,vol2_filtered1);
    linkaxes(ax1,'x');
end

%%
BB_i_start = I.starti;
BB_i_mid = I.midi;
BB_i_end = I.endi;

%% Calculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear','extrap');
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear','extrap');  
VTmean = nanmean(volupper-vollower);
if ~minimum_figs
     
    ax1(2)=subplot(3,1,2); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(3)=subplot(3,1,3); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;

%% Find closer local minima/maxima
searchlength=0.5; %sec
thresholdformove = 0.01*VTmean;
for i=1:length(BB_i_start)
    lefti=BB_i_start(i)-round(searchlength/dt);
    righti=BB_i_start(i)+round(searchlength/dt);
    if lefti<1,lefti=1;end
    if i>1
        if lefti<BB_i_mid(i-1),lefti=BB_i_mid(i-1)+1;end
    end
    if righti>length(Vflow),righti=length(Vflow);end
    if righti>=BB_i_mid(i),righti=BB_i_mid(i)-1;end
    Irange=lefti:righti;
    [minvol,index]=min(vol2_filtered1(Irange))
    currentvol=vol2_filtered1(BB_i_start(i));
    if (currentvol-minvol)>thresholdformove
        BB_i_start(i)=lefti+index-1;
        if i>1
        BB_i_end(i-1)=BB_i_start(i);
        end
    end
end

for i=1:length(BB_i_mid)
    lefti=BB_i_mid(i)-round(searchlength/dt);
    righti=BB_i_mid(i)+round(searchlength/dt);
    if lefti<1,lefti=1;end
    if lefti<BB_i_start(i),lefti=BB_i_start(i);end
    if righti>length(Vflow),righti=length(Vflow);end
    if righti>BB_i_end(i),righti=BB_i_end(i);end
    Irange=lefti:righti;
    [maxvol,index]=max(vol2_filtered1(Irange))
    currentvol=vol2_filtered1(BB_i_mid(i));
    if (maxvol-currentvol)>thresholdformove
        BB_i_mid(i)=lefti+index-1;
    end
end

%% Recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear','extrap');
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear','extrap');  
VTmean = nanmean(volupper-vollower);
if ~minimum_figs
     
    ax1(2)=subplot(3,1,2); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(3)=subplot(3,1,3); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;

%% Force V positive, just in case.

VT(VT<0)=0;
Vdot(Vdot<0)=0;
VTi(VTi<0)=0;
VTe(VTe<0)=0;


