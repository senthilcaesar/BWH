function [I,Vdot,VT,Ti,Te,leak,baselinedriftstd]=Vflowanalysisleak(Vflow,time,dt,minimum_figs)
warning('off');
%% Respiratory Trace Analyse Using Volume as the primary method of breath ID
%[Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs)
% Respiratory Trace Analyse
faster=1; %not used yet
Nloops=1;

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

if 1 %detrend and highpass filtering is ok as this is just to find the end-expiratory points for vollower, the trend is later added back
vol2 = vol1-voldetrend1;
else
vol2 = vol1;    
end
   
if ~minimum_figs
    figure(); 
    ax1(1)=subplot(Nloops+1,1,1); plot(time,vol1-voldetrend1);
    ax1(2)=subplot(Nloops+1,1,2); plot(time,vol2);
    linkaxes(ax1,'x');
end

%Filter volume gently for use in analysis
if 0
filter_HFcutoff_butter1 = 5;
filter_order1 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
else
filter_HFcutoff_butter1 = 5;
filter_LFcutoff_butter1 = 1/15;
filter_order1 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));    
end

vol2_filtered1 = filtfilt(B_butterHcut,A_butterHcut,vol2);

voldetrend1 = voldetrend1 + (vol2-vol2_filtered1); %keep track of detrending


%Vflow_filtered1 = filtfilt(B_butterHcut,A_butterHcut,Vflow);
Vflow_filtered1 = [diff(vol2_filtered1)/dt]; Vflow_filtered1=[Vflow_filtered1(1) Vflow_filtered1];
if ~minimum_figs
    ax1(2)=subplot(Nloops+1,1,2); plot(time,vol2_filtered1);
    linkaxes(ax1,'x');
end

%% loop -- estimate start insp values, then detrend vol, then re-estimate
%start insp times...

for XX=1:Nloops
%% detecting peaks

% accompanying filtered flow trace
%Vflow_filtered1 = [diff(vol2_filtered1)/dt]; Vflow_filtered1=[Vflow_filtered1(1) Vflow_filtered1];

[min_list,max_list] = peakdet(-vol2_filtered1,0.001*std(vol2_filtered1)); %in this configuration, max_list(i,1) > minlist(i,1).
BB_i_start=min_list(:,1);
BB_i_mid=max_list(:,1);
if BB_i_start(1)==1;
    BB_i_start(1)=[];
end
if BB_i_start(end)==length(vol2_filtered1);
    BB_i_start(end)=[];
end
if BB_i_start(1)>BB_i_mid(1)
    BB_i_mid(1)=[];
end
while BB_i_mid(end)>BB_i_start(end)
    BB_i_mid(end)=[];
end

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.-',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.-');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end

%estimate normal tidal volume from the distance between upper and lower
%vols.

vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
    
    
if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end
volupperminusvollower = volupper-vollower;
VTmean = nanmean(volupper-vollower);

BB_i_end = BB_i_start(2:end);
BB_i_start = BB_i_start(1:end-1);

VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


%% Remove extreme mini-breaths from SS and JB 

VTi_thres = VTmean*0.05;
VTe_thres = VTmean*0.05;
Ti_thres = normalTtot/10;
Te_thres = normalTtot/10;

%If BOTH inspiration and expiration are either short or small, remove the breath
criteria=find(((VTi<VTi_thres)|(Ti<Ti_thres))&((VTe<VTe_thres)|(Te<Te_thres)));

if ~isempty(criteria)
firstbreathisbad=criteria(1)==1;
while firstbreathisbad
   criteria(1)=[];
   BB_i_start(1)=[];
   BB_i_mid(1)=[];
   BB_i_end(1)=[];   
   criteria=criteria-1;   
   if isempty(criteria)
       break
   end
   firstbreathisbad=criteria(1)==1;
end
   criteria2=criteria-1;
   BB_i_start(criteria)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria2)=[]; 
end

%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
  
if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
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
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));


%% If we find a  short or small inspiration, presume it is part of the previous expiration:
criteria=find((Ti<Ti_thres)|(VTi<VTi_thres)|((Ti<Ti_thres*2)&(VTi<VTi_thres*2)));
%firstbreathisbad=sum(criteria==1);
if ~isempty(criteria)
firstbreathisbad=criteria(1)==1;
while firstbreathisbad
   criteria(1)=[];
   BB_i_start(1)=[];
   BB_i_mid(1)=[];
   BB_i_end(1)=[];   
   criteria=criteria-1;   
   if isempty(criteria)
       break
   end
   firstbreathisbad=criteria(1)==1;
end
   criteria2=criteria-1;
   BB_i_start(criteria)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria2)=[]; 
end

%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
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
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));

%% If we find a small expiration presume the breath is part of one continuous inspiration:
% if 0 %issues 
% criteria=find(VTe<VTe_thres);
% 
% if ~isempty(criteria)
% NN=length(BB_i_start);
% lastbreathisbad=criteria(end)==NN;
% while lastbreathisbad
%         criteria(end)=[];
%         BB_i_start(end)=[];
%         BB_i_mid(end)=[];
%         BB_i_end(end)=[];
%         NN=NN-1;
%        if isempty(criteria)
%             break
%        end
%         lastbreathisbad=criteria(end)==NN;        
% end
%    criteria2=criteria+1;
%    BB_i_start(criteria2)=[];
%    BB_i_mid(criteria)=[];
%    BB_i_end(criteria)=[];
% end
% 
% %% recalculate breath info:
% vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
%     vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
%     vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
% volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
%     volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
%     volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
% 
% if ~minimum_figs
%      
%     ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
%     ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
%     linkaxes(ax1,'x');
% end
% volupperminusvollower = volupper-vollower;
% VT = volupperminusvollower(BB_i_mid);
% VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
% VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
% Ti = (BB_i_mid-BB_i_start)'*dt;
% Te = (BB_i_end-BB_i_mid)'*dt;
% Ttot = Ti+Te;
% 
% criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
% Vdot = VT./Ttot;
% normalTtot=median(Ttot(criteriafornormal));
% normalTi=median(Ti(criteriafornormal));
% normalTe=median(Te(criteriafornormal));
% 
% 
% chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));


%% If we find a small(er) expiration above FRC presume the breath is part of one continuous inspiration:

VolStartInsp = (VTi-VTe)/VTmean; 
FVTe = VTe./VTi;
criteria=find((VolStartInsp>0.33&(FVTe<0.2)&Te<normalTe/2));

if ~isempty(criteria)
NN=length(BB_i_start);
lastbreathisbad=criteria(end)==NN;
while lastbreathisbad
        criteria(end)=[];
        BB_i_start(end)=[];
        BB_i_mid(end)=[];
        BB_i_end(end)=[];
        NN=NN-1;
       if isempty(criteria)
            break
       end
        lastbreathisbad=criteria(end)==NN;        
end
   criteria2=criteria+1;
   BB_i_start(criteria2)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria)=[];
end

%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
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
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));




%% If there is small(er) inspiration after previous incomplete exhalation (mid-inspiratory pause), merge breath with previous expiration:

VolStartInsp = (VTi-VTe)/VTmean; VolStartInsp = [NaN VolStartInsp(1:end-1)];
criteria=find((VolStartInsp>0.33&VTi<(VTi_thres*4)));


if ~isempty(criteria)
firstbreathisbad=criteria(1)==1;
while firstbreathisbad
   criteria(1)=[];
   BB_i_start(1)=[];
   BB_i_mid(1)=[];
   BB_i_end(1)=[];   
   criteria=criteria-1;   
   if isempty(criteria)
       break
   end
   firstbreathisbad=criteria(1)==1;
end
   criteria2=criteria-1;
   BB_i_start(criteria)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria2)=[]; 
end

%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
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
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));





%% More
VE = VTe./Ttot;
VI = VTi./Ttot;

%% Long expirations are apneas, which are then broken up into smaller pieces
    %Need to add criteria for negative flow
if 1    
    FtotIsAnApnea = 2;
    VEaverage = VTmean/normalTtot;
    TeThresholdForApnea = (0.5*normalTe+FtotIsAnApnea*normalTtot);
    Apnea_B = Te'>TeThresholdForApnea;
    Iapnea = find(Apnea_B);
    for i=length(Iapnea):-1:1
        di=BB_i_end(Iapnea(i))-(BB_i_mid(Iapnea(i))+round(normalTe/dt));
        BB_i_start = [BB_i_start(1:Iapnea(i));BB_i_mid(Iapnea(i))+round(normalTe/dt);BB_i_start(Iapnea(i)+1:end)];
        BB_i_end = [BB_i_end(1:Iapnea(i)-1);BB_i_mid(Iapnea(i))+round(normalTe/dt);BB_i_end(Iapnea(i):end)];
        BB_i_mid = [BB_i_mid(1:Iapnea(i));BB_i_start(Iapnea(i)+1)+round(normalTi/normalTtot*di);BB_i_mid(Iapnea(i)+1:end)];
        Apnea_B = [Apnea_B(1:Iapnea(i)-1);0;1;Apnea_B(Iapnea(i)+1:end)];
    end

%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
Vdot = VT./Ttot;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;


%% Divide up apneas into Ttot-sized parcels
Iapnea = find(Apnea_B);
Nbreathsapnea=round(Ttot(Iapnea)/normalTtot);
for i=length(Iapnea):-1:1
    di=round((BB_i_end(Iapnea(i))-BB_i_start(Iapnea(i)))/Nbreathsapnea(i));
%     time(BB_i_start(Iapnea(i)))
    temp = (BB_i_start(Iapnea(i)):di:((BB_i_end(Iapnea(i)))+di/5))'; %di/5 just handles rounding issues
    BB_i_start = [BB_i_start(1:Iapnea(i));temp(2:(end-1));BB_i_start((Iapnea(i)+1):end)];
    BB_i_end = [BB_i_end(1:Iapnea(i)-1);temp(2:(end-1));BB_i_end(Iapnea(i):end)];
    BB_i_mid = [BB_i_mid(1:Iapnea(i)-1);temp(1:(end-1))+round((normalTi/normalTtot*di));BB_i_mid((Iapnea(i)+1):end)];
    Apnea_B = [Apnea_B(1:Iapnea(i));1+0*temp(2:(end-1));Apnea_B((Iapnea(i)+1):end)];
end


%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
VE = VTe./Ttot;
VI = VTi./Ttot;
Vdot = VT./Ttot;
end

%% Additional information
baselinedriftstd = std(vollower+voldetrend1);
% Endexpvol = vol2_filtered1(BB_i_start);
% Endexpvol_t = time(BB_i_start);
temp = polyfit(time,vollower,1);
leak = temp(1) + leak;
if ~minimum_figs
figure(99);
plot(time,voldetrend1+vol2_filtered1,time(BB_i_start),voldetrend1(BB_i_start)+vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),voldetrend1(BB_i_mid)+vol2_filtered1(BB_i_mid),'k.',time,voldetrend1+volupper,'k:',time,voldetrend1+vollower,'r:');
title(['leak=' num2str(leak) ', baselinedriftstd=' num2str(baselinedriftstd)]);
set(gcf,'position',[50 100 500 500]);
end
%% detrend and re-run
if XX<Nloops
    %vol2_filtered1_backup = vol2_filtered1;
    %vol2_filtered1
    volarray = vol2_filtered1([BB_i_start;BB_i_end(end)]);
    %volarray = [volarray(1) volarray volarray(end)];
    %medfilt1
    %volarray2 = medfilt1(volarray,1)
    %volarray2 = volarray2(2:end-1);
    voltrend = interp1(time([BB_i_start;BB_i_end(end)]),volarray,time,'linear');        
    voltrend(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1));
    voltrend(BB_i_end(end):end)=vol2_filtered1(BB_i_end(end));
    %figure(102);plot(time,[vol2_filtered1;voltrend]);
    vol2_filtered1=vol2_filtered1-voltrend;
end


end
%% Force V positive, just in case.

VT(VT<0)=0;
Vdot(Vdot<0)=0;
VTi(VTi<0)=0;
VTe(VTe<0)=0;
VE(VE<0)=0;
VI(VI<0)=0;

%% Start/End breath indices in original indices
clear I;
I.starti = BB_i_start;
I.midi = BB_i_mid;
I.endi = BB_i_end;


