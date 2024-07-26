function [Vdot_intended,parameters,rsquared]=PesToVflowFit(Flow,Pes,Time,Ccw)

%Flow
%Pes, esophageal pressure
%Time
%Ccw is (estimated) chest wall compliance, should be ~0.2L/cmH2O [set to Inf for no effect].

%Code is written by Scott Sands (May 2016) based on the idea by Peter Catcheside (presented Feb 2013) to transform Pes into an intended ventilation waveform

dt = Time(2)-Time(1);
    
%Calculate Pcw, and filter it (optional). 
Pcw = -1/Ccw*cumsum(Flow)*dt;
    %low pass, 8Hz
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    Pcw = filtfilt(B_butter0,A_butter0,Pcw); 
    if 1
    %high pass, baseline removal = 30 s
        Pcw=detrend(Pcw);
        filter_LFcutoff_butterX = 1/60;
        filter_orderX = 2;
        [B_butterX,A_butterX] = butter(filter_orderX,filter_LFcutoff_butterX/(1/dt/2),'high');
        Pcw = filtfilt(B_butterX,A_butterX,Pcw);
    else
        Pcw=detrend(Pcw);
    end

%Add chest wall pressure to esophageal pressure to get inspiratory muscle pressure:    
Pmus = -Pes - Pcw; %upwards is negative pressure.

%create a structure variable to pass to the PesToVflowModel() function
xdata.data = Pmus; 
xdata.time = Time;
xdata.dt = dt;
ydata = Flow;
    
%setup curve fitting: 
    lsqoptions=optimset('display','off');
    parameters=[6 6 0 0 10]; %guess starting values: R E initialVL Vleak 
    lower=[0.1 0.1 -0.5 -0.1 0]; %lower and upper bounds, edit these as needed to avoid floor/ceiling effects
    upper=[60 60 0.5 0.2 100]; 
    
    %run in a loop, each time "parameters" contains the initial starting points. 
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@PesToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [Vdot_intended] = PesToVflowModel(parameters,xdata);
        
    %goodness of fit
    rsquared = 1 - sum((Vdot_intended-ydata).^2)/sum((ydata-mean(ydata)).^2);
    
    %"parameters" contains the final parameters needed to convert new Pmus data into new intended ventilation data.
    parameters(4)=0; %set leak to zero in the intended vol trace
     
    baselinePmus = prctile(Pmus,parameters(5));
    
% figure();
% set(gcf,'color',[1 1 1]);
% prettyfig = @()set(gca,'box','off','fontname','arial narrow');
% ax(1)=subplot(2,2,1);
% plot(xdata.time,ydata,'k'); 
% hold('on');
% plot(xdata.time,Vdot_intended,'r');
% hold('off');
% prettyfig();
% ax(3)=subplot(2,2,3); 
% plot(xdata.time,-Pmus+baselinePmus,'k',xdata.time,-Pmus-Pcw+baselinePmus,'b',xdata.time,-Pcw,'g');
% prettyfig();
% linkaxes(ax,'x');
% subplot(1,2,2); 
% plot(ydata,Vdot_intended,'k.','markersize',4);
% prettyfig();