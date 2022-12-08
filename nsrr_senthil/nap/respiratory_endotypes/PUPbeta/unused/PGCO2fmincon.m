function [F,PCO2est,Q] = PGCO2fmincon(x,time,VI,VDE,PCO2,timePCO2,VTeatPCO2times)
%start=[Gtot tau1 delay start_P tau2 PCO2_];

Gplant=x(1);
tau1=x(2);
start_P=x(3);
PCO2_=x(4);
VI_=x(5); %VI_=mean(VI);
VDcorrection=x(6);

PCO2est=zeros(length(time),1); %start_P
PCO2est(1)=start_P;
dt=time(2)-time(1);

%VI=VI-VDcorrection*VDE;

if length(x)<7; %tau2
    %'first order'
    for i=1:length(time)-1
        PCO2est(i+1)=PCO2est(i)+dt*(-Gplant*(VI(i)-VI_)-(PCO2est(i)-PCO2_)/tau1);
    end
else
    tau2=x(8);
    taua=1/(1/tau1+1/tau2);
    taub=tau1*tau2;
    G1=Gplant;
    G2=x(7);
    PCO2est(2)=start_P;
    for i=2:length(time)-1
        PCO2est(i+1)=PCO2_+1/(1/dt^2+1/dt*taua)*((VI(i)-VI(i-1))/dt*(-(G1+G2))+(VI(i)-VI_)*(-G1/tau2-G2/tau1)-(PCO2est(i)-PCO2_)*(-1/dt*taua-2/dt^2+1/taub)-(PCO2est(i-1)-PCO2_)*1/dt^2);
    end
    %PCO2est(i+1)=PCO2_+1/(1/dt^2+1/dt*taua)*((VI(i)-VI(i-1))/dt*(-(G1+G2))+(VI(i)-VI_)*(-G1/tau2-G2/tau1)-(PCO2est(i)-PCO2_)*(-1/dt*taua-2/dt^2+1/taub)-(PCO2est(i-1)-PCO2_)*1/dt^2);
end

PCO2est_rs=interp1(time,PCO2est,timePCO2,'linear','extrap');
    
    Err0=PCO2est_rs-PCO2;
    Err0_XHz=interp1(timePCO2,Err0,time,'linear','extrap');
    Err0_XHz_detrended=detrend(Err0_XHz)+mean(Err0_XHz);
    Err=interp1(time,Err0_XHz_detrended,timePCO2,'linear','extrap');
      
    %Err=detrend(Err0)+mean(Err0);
    PCO2est = PCO2est+Err0_XHz_detrended'-Err0_XHz';


% Err=Err1./VIatPCO2times;
if 0
    Q=max(VTeatPCO2times)./VTeatPCO2times;
elseif 0
    upper=0.67*mean(VTeatPCO2times); %VTE>upper is penalised fully | VTE between lower and upper is penalized on a linear scale
    lower=0.33*mean(VTeatPCO2times); %VTE<lower is not penalised at all.
    Q=0*VTeatPCO2times+1; %default=penalise Q=1;
    m=1/(upper-lower);
    c=-lower*m;  
    Q(VTeatPCO2times<upper)=m*VTeatPCO2times(VTeatPCO2times<upper)+c; %y=mx+c
    Q(VTeatPCO2times<=lower)=0;

%     Q(VIatPCO2times==1)=0; %if apnea breath
%     Q(1)=0;
    Thres=-5;
    Q(Err<Thres)=1;
else
    Q=0*Err+1;
end

Npenalisedbreaths=sum(Q);

if 1
cost=Err.*Q;

threscost=prctile(cost(Q>0),90);

if 1
cost_max=(PCO2-mean(PCO2)).*Q;
SSE=(sum((cost(cost<threscost).^2))); % sum square error for penalised breaths
SStot=(sum((cost_max(cost<threscost).^2)));
F=SSE/SStot;
else
F=(sum((cost(cost<threscost).^2)))/Npenalisedbreaths; % sum square error for penalised breaths    
end

Q(cost>threscost)=-0.2;
else
F=median(abs(Err.*Q)); % sum square error for penalised breaths    
end

%%
%PCO2est(i+1)=PCO2est(i)+dt*(-Gplant*(VI(i)-VI_)-(PCO2est(i)-PCO2_)/tau1);
%dP/dt=-G*VI-P/tau;
%sP+P/tau=-G*VI;
%P/VI=-G/(s+1/tau);

%P/VI=-G1/(s+1/tau1)-G2/(s+1/tau2);
%P/VI=[-G1*(s+1/tau2)-G2(s+1/tau1)]/[(s+1/tau1)(s+1/tau2)];
%(s+1/tau1)(s+1/tau2) = s^2+s*(1/tau1+1/tau2)+1/(tau1*tau2);
%taua=1/(1/tau1+1/tau2);
%taub=tau1*tau2;
%P*[s^2+s*1/taua+1/taub]=VI*[-G1*s-G1*1/tau2-G2*s-G2*1/tau1]
%P*[s^2+s*1/taua+1/taub]=VI*[-(G1+G2)*s-G1/tau2-G2/tau1]
%P*s^2+P*s/taua+P/taub=VI*s[-(G1+G2)]+VI*[-G1/tau2-G2/tau1]

%VI*s<->(VI(i)-VI(i-1))/dt
%P*s^2=(dP/dt(i)-dP/dt(i-1))/dt...
%dP/dt(i)=(P(i+1)-P(i))/dt
%dP/dt(i-1)=(P(i)-P(i-1))/dt
%P*s^2=(P(i+1)-2*P(i)+P(i-1))/dt^2

%(P(i+1)-2*P(i)+P(i-1))/dt^2
%+(P(i+1)-P(i))/dt/taua
%+P(i)/taub
%=(VI(i)-VI(i-1))/dt*[-(G1+G2)]+(VI(i)-VI_)*[-G1/tau2-G2/tau1]

%P(i+1)/dt^2+P(i+1)/dt*taua
%-P(i)/dt*taua-2*P(i)/dt^2+P(i)/taub
%+P(i-1)/dt^2
%=(VI(i)-VI(i-1))/dt*[-(G1+G2)]+(VI(i)-VI_)*[-G1/tau2-G2/tau1]

%P(i+1)*[1/dt^2+1/dt*taua]
%=(VI(i)-VI(i-1))/dt*[-(G1+G2)]+(VI(i)-VI_)*[-G1/tau2-G2/tau1]
%-P(i)*[-1/dt*taua-2/dt^2+1/taub]
%-P(i-1)*1/dt^2









