
Time = DataEventHypnog_Mat(:,1);
dt = Time(2)-Time(1)
PmaskX = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pmask')));
FlowX = -DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')));
Lagaverage2=SyncFlowExact(PmaskX,FlowX,Time,0.9);


secperhrlag = median(diff(Lagaverage2))/dt*3600;



%%
figure(1)
Time2 = Time + Lagaverage2;
dt2 = Time2(2)-Time2(1);

Fserror = (dt/dt2) - 1;
ax(1)=subplot(2,1,1);   plot(Time2,FlowX*10,Time,PmaskX);
ax(2)=subplot(2,1,2);   plot(Time2,FlowX);
linkaxes(ax,'x')

%%

DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow'))) = interp1(Time2,DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow'))),Time);
DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pes'))) = interp1(Time2,DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pes'))),Time);

