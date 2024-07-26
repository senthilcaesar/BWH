
    
if 1
        [file,dir] = uigetfile();
        filedir = [dir, file];
        matObj = matfile(filedir);
        varlist = who(matObj);
        load(filedir);   
             
    [~,subj,~] = fileparts(filedir);
    if strcmp(subj(end-2:end),'XHz')
        subj = subj(1:end-4); %strip off _XHz
    end
    if strcmp(subj(end-2:end),'nea')
        subj = subj(1:end-8); %strip off _XHz
    end
end
        [~,~,raw]=xlsread('E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018\PUPstart\AMasterSpreadsheet','S4:T63')
        raw = string(raw);
        newdir='E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018'
  %idx=(strcmp(subj,raw(:,2)))
%   
%   for m=1:60
%     subj = char(raw(m,2));
%     try
%     load([newdir '\Converted\' subj '_XHz.mat'])
%    close all
 
%     
%     idcs   = strfind(dir,'\');
%     newdir = dir(1:idcs(end-1)-1);
%     load([newdir '\CPAPTraits\' subj '_Veupnea.mat'])
    

    oldformat=1
    
        %%

%select periods of normal eupneic ventilation on CPAP during NREM sleep without arousals
%avoid flow limitation, post-arousal overshoot/undershoot, REM

if oldformat


Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1)); %EEG selection uses scoring
Flow=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %EEG selection uses scoring
Pmask=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pmask')==1)); %EEG selection uses scoring
Epochs=DataEventHypnog_Mat(:,strcmp(ChannelsList,'Epochs')==1); %EEG selection uses scoring
EventsAr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1)); %EEG selection uses scoring
EventsResp=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsResp')==1)); %EEG selection uses scoring
Edi=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Edi')==1)); %EEG selection uses scoring
WPr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'WPr')==1)); %EEG selection uses scoring
ArPr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'ArPr')==1)); %EEG selection uses scoring
Pes=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pes')==1)); %EEG selection uses scoring
GGpmax=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'GGpmax')==1)); %EEG selection uses scoring

Position=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Position')==1)); %EEG selection uses scoring

signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','Edi','[WPr ArPr]'};
else
    Time=SigT.Time;
    Flow=SigT.Flow;
    Pmask=SigT.Pmask;
    Edi=SigT.Edi;
    signallist={'[SigT.Epochs SigT.EventsAr 1+SigT.Position]','SigT.EventsResp','SigT.Pmask','SigT.Flow','SigT.Edi','[SigT.WPr SigT.ArPr]'};

end
dt=Time(2)-Time(1);





% 
% if Veupnea.ranges(1,1) < Time(1)
% Veupnea.ranges = Veupnea.ranges +Time(1)
% end

%signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','Edi','[Pes Pmus]','GGpmax','[WPr ArPr]'};

global xvalues yvalues range ax2
clearvalues=0;
range=10*60; 
 xvalues=[]
 yvalues=[]
PlotAndSelectData(signallist,Time,clearvalues);


%%
if 1
Veupnea.ranges=xvalues;

saveontherun=0;
if saveontherun
    save(filenameanddir,'Veupnea','-append');
end
end
%% Calculate Veupnea

%Edi=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Edi')==1)); %EEG selection uses scoring
Drive = Edi;
% Time=(0:length(Flow.values)).*Flow.interval
% dt=Flow.interval

for i=size(Veupnea.ranges,1):-1:1
    if Veupnea.ranges(i,2)<Veupnea.ranges(i,1)+10
        Veupnea.ranges(i,:)=[];
    end
end
minimum_figs=0;
clear T0 Veupnea_n Veupnea_Ttot_n Veupnea_deltaEdi_n
for i=1:size(Veupnea.ranges,1)
    irange=(Time>Veupnea.ranges(i,1)-2)&(Time<=Veupnea.ranges(i,2)+2);
    %[Ix,VI,VTi,Ti,Te]=Vflowanalysis1(Flow(irange),Time(irange),dt,minimum_figs);
    [Ix,VI,VTi,Ti,Te]=Vflowanalysis2(Flow(irange),Time(irange),dt,minimum_figs);
    Ttot=Ti+Te;
    T0(i)=sum(Ttot);
    Veupnea_n(i)=sum(VI.*Ttot)/T0(i);
    Veupnea_Ttot_n(i)=mean(Ttot);
    CPAP_baseline(i)=mean(Pmask(irange))
    
%     
%         close(3)
%       
%         figure(3)
%         try
%         if i>1
%         subplot(3,1,1); plot(Time(irange),Flow(irange),'k',[Veupnea.ranges(i,1) Veupnea.ranges(i,1)],[min(Flow(irange)) max(Flow(irange))],'r',[Veupnea.ranges(i,2) Veupnea.ranges(i,2)],[min(Flow(irange)) max(Flow(irange))],'r');
%     %    subplot(3,1,2); plot(Time(irange),vol1,'k',[Veupnea.ranges(i,1) Veupnea.ranges(i,1)],[min(vol1(irange)) max(vol1(irange))],'r',[Veupnea.ranges(i,2) Veupnea.ranges(i,2)],[min(vol1(irange)) max(vol1(irange))],'r');
%         subplot(3,1,3); plot(Time(irange),Drive(irange),'k');
%         title(i)'
%         end
%         catch me1
%         end
    
            DriveTemp = Drive(irange);
            DriveTemp(Ix.starti(:,1));
            clear DriveMaxB DriveMaxBi DriveMinB
            for ii=1:size(Ix.starti,1)-1
                jrange = Ix.starti(ii,1):Ix.starti(ii+1,1);
                [DriveMaxB(ii,1),idx] = max(DriveTemp(jrange))
                DriveMaxBi(ii,1) = idx + jrange(1) - 1;
            end
            
            for ii=1:size(Ix.starti,1)-1
                if ii==1
                    li = 1;
                else
                    li = DriveMaxBi(ii-1);
                end
                irange2 = li:DriveMaxBi(ii);
                [DriveMinB(ii,1)] = min(DriveTemp(irange2));
            end
            
            
            deltaEdi=DriveMaxB-DriveMinB
            Veupnea_deltaEdi_n(i)= mean(deltaEdi)
    
    
    
    
    
end

Veupnea.mean=sum(Veupnea_n.*T0)/sum(T0)*60;
Veupnea.duration=sum(T0)/60;
Veupnea.cov=std(Veupnea_n)/mean(Veupnea_n);
Veupnea.deltaEdi = mean(Veupnea_deltaEdi_n);

%% Save
%save(filenameanddir,'Pcrit','-append');
if 0
idcs   = strfind(dir,'\');
 newdir = dir(1:idcs(end-1)-1);
 if ~(exist([newdir '\CPAPTraits'], 'dir') == 7)
        mkdir([newdir '\CPAPTraits']);
 end
end

save([newdir '\CPAPTraits' '\' subj '_Veupnea'],'Veupnea');
% saveas(12,[newdir '\CPAPTraits' '\' subj '_Pcrit'],'fig');
% saveas(13,[newdir '\CPAPTraits' '\' subj '_PcritMedian'],'fig');
% saveas(12,[newdir '\CPAPTraits' '\' subj '_Pcrit'],'png');
% saveas(13,[newdir '\CPAPTraits' '\' subj '_PcritMedian'],'png');
    

