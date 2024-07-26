function [X2,X2_O2,FVAL,FVAL_O2,meanPCO2,minPCO2,maxPCO2,stdPCO2,meanPO2,minPO2,maxPO2,stdPO2,PCO2est_win,PO2est_win,time_dt_win,figurehandle] = PlantgainO2andCO2(Vflow,PETCO2,PETO2,VI,VTe,time,BB,BB_i_start,BB_t_PETCO2fix,PETCO2_Bfix,PETO2_Bfix,BB_i_mid);


%% Process VE and PCO2 data from CO2responsesSpike
%clear VI_rs VTe_rs VAI_rs
timestart=max([time(BB_i_start(1)) BB_t_PETCO2fix(1)]);
time0=timestart;
timeend=time(BB_i_start(end))+BB(end)-timestart;

% if timeend>300
%     timeend=300
% end
time_B=time(BB_i_start)-timestart;

dt=0.05;
time_dt=0:dt:timeend; %new time

% VI=detrend(VI)+mean(VI);
VDE=0.15./BB;

for i=1:length(time_dt);
    %yi = y(find(x<=xi,1,'last'));
    VI_rs(i) = VI(find(time_B<=time_dt(i),1,'last')); %Using VE here since this value is directly measured...
    VTe_rs(i) = VTe(find(time_B<=time_dt(i),1,'last')); 
    VDE_rs(i) = VDE(find(time_B<=time_dt(i),1,'last')); 
    %VAI_rs(i) = VAI(find(time_B<=time_dt(i),1,'last')); 
end
% figure(102); ax102(1)=subplot(211);plot(time_B,60*VI,'.',time_dt,60*VI_rs,'r'); set(gca,'fontname','arial narrow');
% figure(101); plot(time_B,VTe,'.',time_dt,VTe_rs,'r');

if 0
    VI_rs=VAI_rs; %overwrite
end


timePCO2=BB_t_PETCO2fix-timestart;
P_CO2=PETCO2_Bfix;
P_O2=PETO2_Bfix;
%Expired volume trace to predict reliability of signal
%VIatPCO2times=interp1(time_dt,VI_rs,timePCO2-dt,'linear','extrap'); %unused currently
VTeatPCO2times=interp1(time_dt,VTe_rs,timePCO2-dt,'linear','extrap');

% figure(101);
% ax101(1)=subplot(411);plot(time-timestart,Vflow,time(BB_i_start)-timestart,Vflow(BB_i_start),'r.',time(BB_i_mid)-timestart,Vflow(BB_i_mid),'k.'); 
% ax101(2)=subplot(412);plot(time_B,VTe,'.',time_dt,VTe_rs,'k',time_dt,VI_rs,'k:',timePCO2,VTeatPCO2times,'r.');
% ax101(3)=subplot(413);plot(time-timestart,PETCO2,timePCO2,P_CO2,'.-');
% if length(PETO2)>0
% ax101(4)=subplot(414);plot(time-timestart,PETO2,timePCO2,P_O2,'.-');
% end
% linkaxes(ax101,'x');

%% Divide into windows and analyze...
T0=timeend;
Twin=T0; %T0 %window width
Fwin=0.1;%Fshift (1-overlap)
stepsize=Twin*Fwin;
%Twin+Nwin*stepsize<T0
Nwin=1+floor((T0-Twin)/stepsize);
Tstart=0:stepsize:Nwin*stepsize;
Tend=Tstart+Twin;
clear X2_win X2_win_O2 FVAL FVAL_O2 PGO2_30 PGCO2_30 PGO2_60 PGCO2_60 startdelay COV_VI PETCO2est_winX
F=round((1/Twin)*1000)/1000:0.001:0.05;
s=-1i*2*pi*F;
Mwin=round(mean((Tend-Tstart))/dt);
for i=1:Nwin
    i
time_dt_win=time_dt(time_dt>=Tstart(i)&time_dt<=Tend(i));
VI_rs_win=VI_rs(time_dt>=Tstart(i)&time_dt<=Tend(i));
VDE_rs_win=VDE_rs(time_dt>=Tstart(i)&time_dt<=Tend(i));
% 
% if 1
%    VI_rs_win(time_dt_win>((Tend(i)-Tstart(i))/2+Tstart(i)))=0.5; 
% end

timePCO2_win=timePCO2(timePCO2>=Tstart(i)&timePCO2<=Tend(i));
P_CO2_win=P_CO2(timePCO2>=Tstart(i)&timePCO2<=Tend(i));
if length(P_O2)>0
P_O2_win=P_O2(timePCO2>=Tstart(i)&timePCO2<=Tend(i));
else
P_O2_win=-P_CO2_win;    
end
VTeatPCO2times_win=VTeatPCO2times(timePCO2>=Tstart(i)&timePCO2<=Tend(i));

time_temp=time-timestart;
time_win=time_temp(time_temp>=Tstart(i)&time_temp<=Tend(i));
PETCO2_win=PETCO2(time_temp>=Tstart(i)&time_temp<=Tend(i));
if length(P_O2)>0
PETO2_win=PETO2(time_temp>=Tstart(i)&time_temp<=Tend(i));
else
PETO2_win=-PETCO2_win;    
end

Vflow_win=Vflow(time_temp>=Tstart(i)&time_temp<=Tend(i));
order=1;
MaxFunEvals=1000;
[PCO2est_win,X2,FVAL(i),EXITFLAG,OUTPUT,Q] = plant_gain_CO2(time_dt_win,VI_rs_win,VDE_rs_win,P_CO2_win,timePCO2_win,VTeatPCO2times_win,order,MaxFunEvals);
[PO2est_win,X2_O2,FVAL_O2(i),EXITFLAG_O2,OUTPUT_O2,Q_O2] = plant_gain_CO2(time_dt_win,VI_rs_win,VDE_rs_win,150-P_O2_win,timePCO2_win,VTeatPCO2times_win,order,MaxFunEvals);
PO2est_win=150-PO2est_win;
% figurehandle=figure(102); 
figurehandle=NaN; 
% ax102(1)=subplot(411);plot(time_win,60*Vflow_win-60*mean(Vflow_win),'k'); set(gca,'fontname','arial narrow');
% set(gca,'Xtick',[],'box','off','Xcolor',[1 1 1])
% ax102(2)=subplot(412);plot(time_dt_win,60*VI_rs_win,'k'); set(gca,'fontname','arial narrow');
% set(gca,'Xtick',[],'box','off','Xcolor',[1 1 1])
% %figure(102); ax102(2)=subplot(212);plot(timePCO2_win,P_CO2_win,'.-',time_dt_win,PCO2est_win,'r'); linkaxes(ax102,'x');
% set(ax102(2),'ylim',[-0.2 60*max(VI_rs_win)*1.1])
X2
X2_O2
X2_win(i,:)=X2;
X2_win_O2(i,:)=X2_O2;
startdelay(i)=max([timePCO2_win(find(Q>0,1))-Tstart(i),timePCO2_win(find(Q_O2>0,1))-Tstart(i)]);

EXITFLAG_win(i)=EXITFLAG;
% ax102(3)=subplot(413); plot(time_win,PETCO2_win,'color',[0.5 0.5 0.5]); hold('on'); set(gca,'fontname','arial narrow');
% plot(timePCO2_win,P_CO2_win,'r.',time_dt_win,PCO2est_win,'r--'); hold('off');
% set(gca,'Xtick',[],'box','off','Xcolor',[1 1 1])
% ax102(4)=subplot(414); plot(time_win,PETO2_win,'color',[0.5 0.5 0.5]); hold('on'); set(gca,'fontname','arial narrow');
% plot(timePCO2_win,P_O2_win,'r.',time_dt_win,PO2est_win,'r--'); hold('off');
% set(gca,'Xtick',[],'box','off','Xcolor',[1 1 1])
% linkaxes(ax102,'x');
COV_VI(i)=std(VI_rs_win)/mean(VI_rs_win);

% figure(104); 
if length(X2_O2)<7
    X2_O2(7:8)=[0,0];
    X2(7:8)=[0,0];
PGO2_30(i)=abs(X2_O2(1)/(j*2*pi/30+1/X2_O2(2))); %for 1 tau only...
PGCO2_30(i)=abs(X2(1)/(j*2*pi/30+1/X2(2))); 
PGO2_60(i)=abs(X2_O2(1)/(j*2*pi/60+1/X2_O2(2))); 
PGCO2_60(i)=abs(X2(1)/(j*2*pi/60+1/X2(2))); 
PGO2f=X2_O2(1)./(s+1/X2_O2(2)); 
PGCO2f=X2(1)./(s+1/X2(2)); 
else
PGO2_30(i)=abs(X2_O2(1)/(j*2*pi/30+1/X2_O2(2))+X2_O2(7)/(j*2*pi/30+1/X2_O2(8))); 
PGCO2_30(i)=abs(X2(1)/(j*2*pi/30+1/X2(2))+X2(7)/(j*2*pi/30+1/X2(8))); 
PGO2_60(i)=abs(X2_O2(1)/(j*2*pi/60+1/X2_O2(2))+X2_O2(7)/(j*2*pi/60+1/X2_O2(8))); 
PGCO2_60(i)=abs(X2(1)/(j*2*pi/60+1/X2(2))+X2(7)/(j*2*pi/60+1/X2(8))); 
PGO2f=X2_O2(1)./(s+1/X2_O2(2))+X2_O2(7)./(s+1/X2_O2(8)); 
PGCO2f=X2(1)./(s+1/X2(2))+X2(7)./(s+1/X2(8)); 
end


% plot(F,abs(PGCO2f),F,abs(PGO2f))

PETCO2est_winX(i,:)=PCO2est_win(1:Mwin);

end

%% Sliding window estimation of true PCO2...
if 0
Nend=round(Tend(end)/dt);
iskip=round(stepsize/dt);

PETCO2est_winNaNs = NaN*zeros(Nwin,Nend);
istart=1;

for i=1:Nwin
    iend=istart+Mwin-1;
    if 0
    delta=round(startdelay(i)/dt);
    else
    delta=0;
    end
    PETCO2est_winNaNs(i,(istart+delta):iend)=PETCO2est_winX(i,(1+delta):end); %round(startdelay(i)/dt)
    %delta=0;
    %PETCO2est_winNaNs(i,istart:iend)=PETCO2est_winX(i,:); %round(startdelay(i)/dt)
    istart=istart+iskip;
end
PETCO2est_winNaNmedian=nanmedian(PETCO2est_winNaNs); %nanmedian
PETCO2est_winNaNmean=nanmean(PETCO2est_winNaNs); %nanmedian
%figure(105),plot(time-timestart,PETCO2,timePCO2,P_CO2,'.-',[1:Nend]*dt,PETCO2est_winNaNmean,[1:Nend]*dt,PETCO2est_winNaNmedian)

%plot(time-timestart,PETCO2,timePCO2,P_CO2,'.-');
end
%%
mean(X2_win_O2)
std(X2_win_O2)/Nwin^0.5
mean(X2_win_O2)
std(X2_win_O2)/Nwin^0.5

PGO2_30_median=median(PGO2_30)/60;
PGCO2_30_median=median(PGCO2_30)/60;
PGO2_60_median=median(PGO2_60)/60;
PGCO2_60_median=median(PGCO2_60)/60;
PG60ratio1=PGO2_60_median/PGCO2_60_median;
PG30ratio1=PGO2_30_median/PGCO2_30_median;
PG60ratio2=median(PGO2_60./PGCO2_60);
PG30ratio2=median(PGO2_30./PGCO2_30);

% figure(103),
% subplot(231),plot(FVAL,PGCO2_30/60,'.');
% subplot(232),plot(FVAL,PGO2_30/60,'.');
% subplot(233),plot(FVAL,PGO2_30./PGCO2_30,'.');
% subplot(234),plot(FVAL,PGCO2_60/60,'.');
% subplot(235),plot(FVAL,PGO2_60/60,'.');
% subplot(236),plot(FVAL,PGO2_60./PGCO2_60,'.');
% 
% figure(104),
% subplot(231),plot(COV_VI,PGCO2_30/60,'.');
% subplot(232),plot(COV_VI,PGO2_30/60,'.');
% subplot(233),plot(COV_VI,PGO2_30./PGCO2_30,'.');
% subplot(234),plot(COV_VI,PGCO2_60/60,'.');
% subplot(235),plot(COV_VI,PGO2_60/60,'.');
% subplot(236),plot(COV_VI,PGO2_60./PGCO2_60,'.');

%%
meanPCO2=mean(PCO2est_win);
minPCO2=min(PCO2est_win);
maxPCO2=max(PCO2est_win);
stdPCO2=std(PCO2est_win);
meanPO2=mean(PO2est_win);
minPO2=min(PO2est_win);
maxPO2=max(PO2est_win);
stdPO2=std(PO2est_win);


time_dt_win = time_dt_win + timestart;


