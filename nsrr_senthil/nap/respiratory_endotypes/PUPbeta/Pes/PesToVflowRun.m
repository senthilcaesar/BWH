function [Vdot_intended]=PesToVflowRun(Flow,Pes,Time,Ccw,parameters)

%Flow
%Pes, esophageal pressure
%Time
%Ccw is (estimated) chest wall compliance, should be ~0.2L/cmH2O [set to Inf for no effect].
%parameters is obtained from prior use of "PesToVflow()"

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
Pmus = -(-Pes - Pcw); %upwards is negative pressure.

%create a structure variable to pass to the PesToVflowModel() function
xdata.data = Pmus; 
xdata.time = Time;
xdata.dt = dt;
ydata = Flow;
        
    [Vdot_intended] = PesToVflowModel(parameters,xdata);
        
    baselinePmus = prctile(Pmus,parameters(5));
    
% figure();
% set(gcf,'color',[1 1 1]);
% prettyfig = @()set(gca,'box','off','fontname','arial narrow');
% ax(1)=subplot(2,1,1);
% plot(xdata.time,ydata,'k'); 
% hold('on');
% plot(xdata.time,Vdot_intended,'r');
% hold('off');
% prettyfig();
% ax(2)=subplot(2,1,2); 
% plot(xdata.time,-Pmus+baselinePmus,'k',xdata.time,-Pmus-Pcw+baselinePmus,'b',xdata.time,-Pcw,'g');
% prettyfig();
% linkaxes(ax,'x');
