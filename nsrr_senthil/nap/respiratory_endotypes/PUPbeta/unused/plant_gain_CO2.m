function [PCO2est,X2,FVAL,EXITFLAG,OUTPUT,Q] = plant_gain_CO2(time,VI,VDE,PCO2,timePCO2,VTeatPCO2times,order,MaxFunEvals) 

Gplant=5;
G2=20;
tau1=8; % =60
PCO2_=mean(PCO2);
VI_=mean(VI);
start_P=PCO2_;
start=[Gplant tau1 start_P PCO2_ VI_ 1];

tau2=0;
if tau2>0
    start=[Gplant tau1 start_P PCO2_ VI_ 1 G2 tau2];
end
x=start;

%% FminCON
upper=[ 400    60    PCO2_*2   PCO2_+10 VI_*1.25 0 20   200]; %1=100% max drift in VE over the window
lower=[ 1     2   0   PCO2_-10 VI_*0.8   0 1    60];
if tau2==0
    upper(7:8)=[];
    lower(7:8)=[];
end
OPTIONS = optimset('Display','off','TolX',1e-17,'TolFun',1e-17,'MaxIter',1e17,'MaxFunEvals',MaxFunEvals,'Algorithm','interior-point');
for i=1:4
    [X2,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) PGCO2fmincon(x,time,VI,VDE,PCO2,timePCO2,VTeatPCO2times),start,[],[],[],[],lower,upper,[],OPTIONS);
    start=X2; % This is the vector of parameters that minimise the function in the fmincon function
end

[F,PCO2est,Q]=PGCO2fmincon(X2,time,VI,VDE,PCO2,timePCO2,VTeatPCO2times);




