function PlotEventData(EvtVE)

PtName_wSpecialChars=EvtVE.FileName;
PtName=char(regexprep(PtName_wSpecialChars,'[^a-z A-Z 0-9]',' ')); % DLM removed special chars to stop format affect in written strings
EvMatrix= EvtVE.EventVE';
VeMatrix= EvtVE.VE';
DsatMatrix=EvtVE.Dsat';
HRMatrix= EvtVE.HR';
ArousalMatrix= EvtVE.Arousal';
ArIntMatrix= EvtVE.ArousalIntensity';
EnTime=EvtVE.EnTime;

  
%% plot p-level data
plotnum=1;

partitions=cell(1,plotnum);
[partitions{1,1:plotnum}] = deal(1/plotnum);

scaled=1;
edgalphaind=0.05;
edgalphaavg=1;

lnwidthind=0.3;
lnwidthavg=3;

figure(13)
set(gcf,'Position',[1020 95 560 850])
clf

% create panel
p = panel();
p.pack('h', partitions)
p.de.margin = 2;
for pl=1:plotnum
    p(pl).pack(6,1);
    p(pl).margin = [10 0 10 0];
end
p.margin = [20 18 5 6];
p.fontsize = 12;

ylabRot=90;
ylabPos=[-0.08, 0.5, 0];

x2=EnTime;
x=EnTime;
if size(VeMatrix,2)>=5
    
    Veplotlims=[x(1) x(end) -0.1 1.1];
    Evplotlims=[x(1) x(end) 0 400];
    Dsplotlims=[x(1) x(end) 80 100];
    Hrplotlims=[x(1) x(end) 40 80];
    Arplotlims=[x(1) x(end) -0.1 1.1];
    ArScoreplotlims=[x(1) x(end) -0.1 9.1];
    
    nn=size(VeMatrix,1)*size(VeMatrix,2);
    if scaled
        if length(find(isnan(VeMatrix)))<nn
            Evplotlims=[x(1) x(end) 100*prctile(prctile(VeMatrix,10),10) 100*prctile(prctile(VeMatrix,90),90)];
        end
        if length(find(isnan(DsatMatrix)))<nn
            if prctile(prctile(DsatMatrix,10),10)==prctile(prctile(DsatMatrix,90),90)
                Dsplotlims=[x(1) x(end) 90 101];
            else
                 Dsplotlims=[x(1) x(end) prctile(prctile(DsatMatrix,10),10) prctile(prctile(DsatMatrix,90),90)];
            end
        end
        if length(find(isnan(HRMatrix)))<nn
            Hrplotlims=[x(1) x(end) prctile(prctile(HRMatrix,10),10) prctile(prctile(HRMatrix,90),90)];
        end
        if length(find(isnan(ArIntMatrix)))<nn
            ArScoreplotlims=[x(1) x(end) prctile(prctile(ArIntMatrix,10),10) prctile(prctile(ArIntMatrix,90),90)];
        end
        
    end
    
    Avgplot= p(1);
    edgealpha=edgalphaind;
    lnwidth=lnwidthind;
    lnstyle='-';
    edgeclr='k';
    
    
    %% Plot AR
    for ii=1:size(ArousalMatrix,2) 
       Avgplot(1,1).select();
       hold on;
       patchline(x,round(ArousalMatrix(:,ii)),'linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);  
    end

    edgealpha=edgalphaavg;
    lnwidth=lnwidthavg;
    patchline(x,nanmean(ArousalMatrix')','linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);

    set(gca, 'xtick', [0]);
    set(gca, 'ytick', [0.5],'YGrid', 'on');
    ylab=ylabel('Arousal','rot',ylabRot);
    set(ylab, 'Units', 'Normalized', 'Position', ylabPos);
    axis(Arplotlims);   
    title(PtName);


    %% Plot Events
    edgealpha=edgalphaind;
    lnwidth=lnwidthind;
    for ii=1:size(EvMatrix,2) 
       Avgplot(2,1).select();
       hold on;
       patchline(x,EvMatrix(:,ii),'linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);  
    end

    edgealpha=edgalphaavg;
    lnwidth=lnwidthavg;
    patchline(x,nanmean(EvMatrix')','linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);
    
    set(gca, 'xtick', [0]);
    set(gca, 'ytick', [0.5],'YGrid', 'on');
    axis(Veplotlims);
    ylab=ylabel('Resp. Events','rot',ylabRot);
    set(ylab, 'Units', 'Normalized', 'Position', ylabPos);
    
   

    %% Plot VE
    edgealpha=edgalphaind;
    lnwidth=lnwidthind;
    for ii=1:size(VeMatrix,2) 
       Avgplot(3,1).select();
       hold on;
       patchline(x2,100*VeMatrix(:,ii),'linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);  
    end

    edgealpha=edgalphaavg;
    lnwidth=lnwidthavg;
    patchline(x2,100*nanmean(VeMatrix')','linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);
    set(gca, 'xtick', [0]);
    set(gca, 'ytick', [100 200],'YGrid', 'on');

    axis(Evplotlims);
    ylab=ylabel('VE(%eup)','rot',ylabRot);
    set(ylab, 'Units', 'Normalized', 'Position', ylabPos);
%     disp([EvtVE.VeFtrs.FsinEst EvtVE.VeFtrs.FexpEst EvtVE.VeFtrs.FsquEst ])
    

    %% Plot Sat
    edgealpha=edgalphaind;
    lnwidth=lnwidthind;
    for ii=1:size(DsatMatrix,2) 
       Avgplot(4,1).select();
       hold on;
       patchline(x,DsatMatrix(:,ii),'linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);  
    end 

    edgealpha=edgalphaavg;
    lnwidth=lnwidthavg;
    patchline(x,nanmean(DsatMatrix')','linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);
    set(gca, 'xtick', [0]);
    set(gca, 'ytick', [90 98],'YGrid', 'on');

    ylab=ylabel('SpO2(%)','rot',ylabRot);
    set(ylab, 'Units', 'Normalized', 'Position', ylabPos);
    axis(Dsplotlims);
    
    %% Plot ARIntensity
    edgealpha=edgalphaind;
    lnwidth=lnwidthind;
    for ii=1:size(ArIntMatrix,2) 
       Avgplot(5,1).select();
       hold on;
       patchline(x,ArIntMatrix(:,ii),'linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);  
    end

    edgealpha=edgalphaavg;
    lnwidth=lnwidthavg;
    patchline(x,nanmean(ArIntMatrix')','linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);
    set(gca, 'xtick', [0]);
    set(gca, 'ytick', [3 6],'YGrid', 'on');
    ylab=ylabel('ArInt','rot',ylabRot);
    set(ylab, 'Units', 'Normalized', 'Position', ylabPos);
    axis(ArScoreplotlims);

    %% Plot HR
    edgealpha=edgalphaind;
    lnwidth=lnwidthind;
    for ii=1:size(HRMatrix,2)  
       Avgplot(6,1).select();
       hold on;
       patchline(x,HRMatrix(:,ii),'linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);  
    end

    edgealpha=edgalphaavg;
    lnwidth=lnwidthavg;
    patchline(x,nanmean(HRMatrix')','linestyle',lnstyle,'edgecolor',edgeclr,'linewidth',lnwidth,'edgealpha',edgealpha);
    set(gca, 'ytick', [50 70],'YGrid', 'on');
    ylab=ylabel('HR(beats/min)','rot',ylabRot);
    set(ylab, 'Units', 'Normalized', 'Position', ylabPos);
    axis(Hrplotlims); 
    
    %% X axis (common for all subplots)
    xlabel('Time(s)')

end

