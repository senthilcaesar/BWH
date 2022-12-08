%% Plot FOT from DataEventHypnog_Mat

DataT = array2table(DataEventHypnog_Mat);
DataT.Properties.VariableNames = ChannelsList;
%%
% DataT.Time = DataT.Time;
% Vflow = DataT.FlowOrig;
% Pmask = DataT.PmaskOrig;
% DataT.Flow = DataT.Flow;
% DataT.Pmask = DataT.Pmask;
% DataT.PflowFOT = DataT.DataT.PflowFOT;
% DataT.PpmaskFOT = DataT.DataT.PpmaskFOT;
% DataT.YFOT = DataT.DataT.YFOT;
% DataT.YFOTclean = DataT.YFOTclean;
% DataT.CohFOT = DataT.DataT.CohFOT;

%%
    figure(10); clf(10); set(gcf,'color',[1 1 1]);
    
    ax3(1)=subplot(4,1,1);
    plot(DataT.Time,DataT.FlowOrig,DataT.Time,2+DataT.PmaskOrig/nanstd(DataT.Pmask)*nanstd(DataT.Flow)); ylabel('Flow/Pmask clean');
    
    hold on
    plot(DataT.Time,DataT.Flow,'k',DataT.Time,2+DataT.Pmask/nanstd(DataT.Pmask)*nanstd(DataT.Flow),'k'); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    
    %ax3(2)=subplot(5,1,2); plot(DataT.Time,FlowHigh,DataT.Time,0.5+PmaskHigh/nanstd(PmaskHigh)*nanstd(FlowHigh)); ylabel('Flow/Pmask HF'); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    ax3(3)=subplot(4,1,2); semilogy(DataT.Time,[(abs(DataT.PflowFOT)) (abs(DataT.PpmaskFOT))]); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    ylabel('Pxx, Pyy');
    ax3(4)=subplot(4,1,3); plot(DataT.Time,DataT.YFOT); ylim([0 1.5]); set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
    hold on
    plot(DataT.Time,[DataT.YFOTclean],'k','linewidth',2.5);
    hold off
    ylabel('|Y=V/P|');
    %ax3(3)=subplot(4,1,3); plot(time2,volume1,time2,detrend(volume1));
    ax3(5)=subplot(4,1,4); plot(DataT.Time,DataT.CohFOT); set(gca,'box','off');
    ylabel('Coh');
    linkaxes(ax3,'x')

