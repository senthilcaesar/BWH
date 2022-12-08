%Idea

%Find delay time that Ve(airway open) is best explained by SpO2 (non-linear plot is ok);
run EupneaFromSpO2

W=length(BreathDataTable{1});
plotfig=0;
clear bSlope Rsq
kval=-3:10;
for kk=1:length(kval)
    k=kval(kk);
    run SpO2VEanalysisScript;
end
plotfig=settings.plotfigs;
if plotfig
    figure(13)
    subplot(2,1,1); plot(kval,Rsq,'.-');
    subplot(2,1,2); plot(kval,bSlope,'.-');
end

[VESpO2slope,kkmax]=max(bSlope);

k=kval(kkmax);

run SpO2VEanalysisScript;
if b(2)<0
    b(2)=0; b(1)=nanmean(yval);
end
VESpO2slope = b(2);
SpO2VRA=b(2)*(100-spo2wake75th)+b(1)-1;% prctile(xval,2.5) spo2wakemean10to90 spo2wake75th
VEarthres = 100*(1+VESpO2slope*(spo2wake75th-SpO2arthres));
if plotfig
    figure(11)
    subplot(3,1,[1 2]);
    text(90,6,['Sensitivity:' num2str(100*VESpO2slope,2) '%eupnea/%SpO2'],'fontsize',8);
    text(90,5.2,['ArThres:' num2str(round(SpO2arthres)) '%SpO2'],'fontsize',8);
    text(90,4.4,['VRA:' num2str(round(100*SpO2VRA)) '%eupnea'],'fontsize',8);
    %SpO2VRA_x = spo2wakemean10to90;
    linkaxes(ax11,'x')
    
    ylim([0 8]);
    xlims=get(gca,'xlim');
    xlims;
    if xlims(1)>89
        xlims(1)=89;
    end
    set(gca,'xlim',xlims);
end
SpO2VEtraits = [VESpO2slope SpO2VRA SpO2arthres k spo2wake75th VEarthres];



