function [Fserror,delay] = FixFsDriftC1A1f(Flow,Pmask,Channels_Fs,ChannelsList) 

chFlow = find(strcmp(ChannelsList,'Flow')==1);
chPmask = find(strcmp(ChannelsList,'Pmask')==1);

Channels_Fs(chFlow);
Channels_Fs(chPmask);

TimeFlow=(0:(1/Channels_Fs(chFlow)):0+(length(Flow)-1)*(1/Channels_Fs(chFlow)))'; % This is the time vector associated with the 100Hz Flow data.  
TimePmask=(0:(1/Channels_Fs(chPmask)):0+(length(Pmask)-1)*(1/Channels_Fs(chPmask)))'; % This is the time vector associated with the 100Hz Flow data.  

Pmask_rs = interp1(TimePmask,Pmask,TimeFlow,'linear');

Lagaverage2=SyncFlowExact(Pmask_rs,-Flow,TimeFlow,0.9);

dt = TimeFlow(2)-TimeFlow(1);
secperhrlag = median(diff(Lagaverage2))/dt*3600;

delay = Lagaverage2(1);
%%
%figure(1)
Time2 = TimeFlow + Lagaverage2;
dt2 = Time2(2)-Time2(1);

%Fserror = (dt/dt2) - 1;
Fserror = (dt/median(diff(Time2))) - 1;
% 
% ax(1)=subplot(2,1,1);   plot(Time2,FlowX*10,Time,PmaskX);
% ax(2)=subplot(2,1,2);   plot(Time2,FlowX);
% linkaxes(ax,'x')

%%
% 
% DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow'))) = interp1(Time2,DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow'))),Time);
% DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pes'))) = interp1(Time2,DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pes'))),Time);

