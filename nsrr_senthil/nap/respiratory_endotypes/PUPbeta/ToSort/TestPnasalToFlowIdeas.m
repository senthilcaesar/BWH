
[time,Vflow,index,time1,VI1,BB_i_start,BB_i_mid,BB_i_end,BB_t,VTi,VTe,VI,VE,Ttot_B,Ti,Te,FLindex1,PeakInspFlow,PeakInspFlow_t,MidInspFlow] = VEfromFlow_sqrt_V14_1(Time,Flow,settings.sqrt_scaling);

sqrt_scaling=1;
FlowSqrt=Flow;
if sqrt_scaling
    FlowSqrt(FlowSqrt>0)=FlowSqrt(FlowSqrt>0).^0.5;
    FlowSqrt(FlowSqrt<0)=-((-FlowSqrt(FlowSqrt<0)).^0.5);
end

figure(100)
plot(Time,FlowSqrt,'k',Time,0*FlowSqrt,'k--'); hold('on');
plot(time,Vflow,'b');
plot(time(BB_i_start),Vflow(BB_i_start),'r.')
plot(time(BB_i_mid),Vflow(BB_i_mid),'g.');

%%
tic
%try a range of leaks, and adjust the I-E amplitude ratio, to acheive the smallest SD in the volume trend signal / smallest calculated leak
FlowSignal=Flow-mean(Flow);

Time;


plotfigs=1;
leaklimit=1;
expinspcorrectionlimit=2;

voltemp = cumsum(FlowSignal);

Ttotest = 10;
 
addpath('C:\Users\szs88\Dropbox\PT+SS\PUPbeta\PUPbeta 20160602 O2psg\env_secant')
view = round(Ttotest*1.5/dt);
[env] = env_secant(Time,voltemp,view,'bottom');

figure(1000)
subplot(1,1,1);plot(Time,voltemp,Time,env);

SDvollower = nanstd(env);

N=11;
M=11;
XX=linspace(-leaklimit,leaklimit,M);
YY=10.^linspace(log10(1/expinspcorrectionlimit),log10(expinspcorrectionlimit),N); %%%%%% beta = E/I magnitude coefficient 

clear SDvollower
for j=1:N %leak
    for i=1:M %amplituderatio
        exponent=0.5;
        alpha = XX(i);
        beta = YY(j);
        iqr=prctile(FlowSignal,75)-prctile(FlowSignal,25);
        FlowSignal2=FlowSignal + alpha*iqr;
        FlowSignal2(FlowSignal2>0)=(FlowSignal2(FlowSignal2>0).^(exponent))/(beta^0.5);
        FlowSignal2(FlowSignal2<0)=(-(-FlowSignal2(FlowSignal2<0)).^(exponent))*(beta^0.5);
        VolSignal = cumsum(FlowSignal2)*dt; 
            VolSignal = detrend(VolSignal);
        VolEnvelope = env_secant(Time,VolSignal,view,'bottom');
        figure(1000)
        subplot(2,1,1);plot(Time,VolSignal,Time,VolEnvelope);
        subplot(2,1,2);plot(Time,FlowSignal2,Time,0*FlowSignal,'k--');
        SDvollower(i,j) = nanstd(VolEnvelope);
    end
end

[temp,mini_] = min(SDvollower);
[~,minj] = min(temp);
mini=mini_(minj);
        exponent=0.5;
        alpha = XX(mini);
        beta = YY(minj);
        iqr=prctile(FlowSignal,75)-prctile(FlowSignal,25);
        FlowSignal2=FlowSignal + alpha*iqr;
        FlowSignal2(FlowSignal2>0)=(FlowSignal2(FlowSignal2>0).^(exponent))/(beta^0.5);
        FlowSignal2(FlowSignal2<0)=(-(-FlowSignal2(FlowSignal2<0)).^(exponent))*(beta^0.5);
        VolSignal = cumsum(FlowSignal2)*dt; 
            VolSignal = detrend(VolSignal);
        VolEnvelope = env_secant(Time,VolSignal,view,'bottom');
        figure(1000)
        subplot(2,1,1);plot(Time,VolSignal,Time,VolEnvelope);
        subplot(2,1,2);plot(Time,FlowSignal2,Time,0*FlowSignal,'k--');

toc
%%
clear ParametersArray Fres1
tic
OPTIONS = optimset('Display','off','TolX',1e-6,'TolFun',1e-6,'MaxIter',50,'MaxFunEvals',50,'Algorithm','interior-point');
%'trust-region-reflective','active-set','sqp','interior-point'
lower = [-leaklimit 0.5];
upper = [leaklimit 2];
ParametersGuess = [0 1;lower;upper;lower(1) upper(2);lower(2) upper(1)]
for i=1:size(ParametersGuess,1)
    Parameters = ParametersGuess(i,:);
    [ParametersArray(:,i),Fres1(i),~,~] = fmincon(@(Parameters) PnasaltoFlow(Parameters,FlowSignal,Time,dt,0),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
end
Fres1
toc

%%
    plotfigs=0;
    leaklimit=1;
    expinspcorrectionlimit=2;

tic
clear final_parameters final_parameters1 Fres Fres1 FresXY
    OPTIONS = optimset('Display','off','TolX',1e-3,'TolFun',1e-3,'MaxIter',10,'MaxFunEvals',10,'Algorithm','interior-point');
    
    N=3;
    lower = [-leaklimit]; %alpha = leak (in units of flow/flowIQR)
    upper = [leaklimit];
    YY=10.^linspace(log10(1/expinspcorrectionlimit),log10(expinspcorrectionlimit),N); %%%%%% beta = E/I magnitude coefficient
    
    for j=1:N
        Parameters=0; %estimate leaks
        [final_parameters1(j),Fres1(j),~,~] = fmincon(@(Parameters) PnasaltoFlow2(Parameters,YY(j),FlowSignal,Time,dt,0),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
        [Fres(j),Parameters_out_manual,~]=PnasaltoFlow([final_parameters1(j),YY(j)],FlowSignal,Time,dt,plotfigs);
    end
    figure(34)
    semilogx(YY,Fres,'.-') %Fres is baseline drift variance
    
    [minFres,j]=min(Fres);

    %OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',20,'MaxFunEvals',20,'Algorithm','interior-point');

    leftj=j-1;
    rightj=j+1;
    
if leftj==0||rightj==N+1
    %do nothing
    final_parameters_=[final_parameters1(j),YY(j)];
else
    for xx=1:2
    xdatarange = YY(leftj:rightj);
    [xpeak,yvalue] = PeakFitQuadratic(xdatarange,-Fres(leftj:rightj));
    if xpeak<xdatarange(2)
        j=leftj+1;
    else
        j=leftj+2;
    end
    YY=[YY(1:j-1) xpeak YY(j:end)];        
    final_parameters1=[final_parameters1(1:j-1) NaN final_parameters1(j:end)];
    Fres1=[Fres1(1:j-1) NaN Fres1(j:end)];
    Fres=[Fres(1:j-1) NaN Fres(j:end)];
    Parameters=0; 
    [final_parameters1(j),Fres1(j),~,~] = fmincon(@(Parameters) PnasaltoFlow2(Parameters,YY(j),FlowSignal,Time,dt,0),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
    [Fres(j),Parameters_out_manual,~]=PnasaltoFlow([final_parameters1(j),YY(j)],FlowSignal,Time,dt,plotfigs);
    [minFres,j]=min(Fres);
    leftj=j-1;
    rightj=j+1;    
    end
end

final_parameters_=[final_parameters1(j),YY(j)];
[F,Parameters_out,Ydata2]=PnasaltoFlow(final_parameters_,FlowSignal,Time,dt,plotfigs);

if 1
    minFreslast=minFres
    OPTIONS = optimset('Display','off','TolX',1e-9,'TolFun',1e-9,'MaxIter',50,'MaxFunEvals',50,'Algorithm','sqp');
    %'trust-region-reflective','active-set','sqp','interior-point'
    lower = [final_parameters_(1)-0.2 final_parameters_(2)/1.2];
    upper = [final_parameters_(1)+0.2 final_parameters_(2)*1.2];
    Parameters = final_parameters_;
    [final_parameters_,minFres,~,~] = fmincon(@(Parameters) PnasaltoFlow(Parameters,FlowSignal,Time,dt,0),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
end

figure(34)
semilogx(YY,Fres,'.-')

% [F_manual,Parameters_out_manual,Ydata3]=PnasaltoFlow([0 0.57],Ydata1,Time(I),dt,plotfigs);

leak=-final_parameters_(1)*(prctile(FlowSignal,75)-prctile(FlowSignal,25));

figure(3);
a3(1) = subplot(2,1,1); plot(Time,FlowSignal,Time,leak+0*Time,'k--');
a3(2) = subplot(2,1,2); plot(Time,Ydata2);
linkaxes(a3,'x')
minFres
    toc
    
    