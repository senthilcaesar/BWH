% Produce Group ensemble plots

function [GroupEnsembles,GroupEnsembleTime,EAinfoGroup, Mrange,nremcount]=plotGroupEnsembles(MrangeOverride,preEvtArCrit,selectState)
%%

global settings AMasterSpreadsheet

if ~exist('preEvtArCrit')
    preEvtArCrit=0;
end
if ~exist('selectState')
    selectState=9; %all states
end
preEvtArWin=15;
thres=10;



%% Load AMasterSpreadsheet
path={'E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018\Analyzed\'}

%[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');
[num,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,1),'UniformOutput',false)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,1:2);

NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));

settings = ImportSettings(settings,AMasterSpreadsheet);

MaxM = size(num,1);

try
    M = max(settings.Mrange);
catch
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end
if exist('MrangeOverride')
    Mrange=MrangeOverride(:)';
else
    Mrange=settings.Mrange;
end

%%
k=1
GroupEnsembles=[];
EAinfoGroup=[];
GroupEnsemblesInfo=[];
success=1;
%% Load analyzed data
for i=Mrange
    i
    try
      
        if 1
        filedir = [path{:} 'EventAnalyzed\' settings.savename '_' num2str(i) '.mat'];
        
        if exist(filedir)==2
            disp(['Processing: ' num2str(i) ': ' settings.savename '_' num2str(i) '.mat']);
            load (filedir,'Ensembles','Boxes','EAinfo', 'Evts');
        end
        end
        
        
        
        if preEvtArCrit
        for j=1:size(Boxes.EventsAr,1)
            ri=round(100-Evts.RespT.EventDuration(j));
            li=ri-preEvtArWin+1;
            if li<1
                continue
            end
            criteria(j,1)=(sum(Boxes.EventsAr(j,li:ri),2)==0);
        end
        else
            criteria=true(size(Boxes.EventsAr,1),1);
        end
        
          if ~exist('Evts.RespT.Epochs')
            try
                Evts.RespT.Epochs=Evts.RespT.state
            catch
            end 
          end
            
        if selectState==1
            hypok=[ 2]
        elseif  selectState==2
            hypok=[ 1 ]
        elseif  selectState==3
            hypok=[0 ]
        elseif  selectState==4
            hypok=[0 1 2]
        elseif  selectState==5
            hypok=[3]
        else
            hypok=[0 1 2 3 4]
        end
            
            
           criteria= criteria & sum(Evts.RespT.Epochs==hypok,2)>0
           EventsUsed(i)= sum(sum(Evts.RespT.Epochs==hypok,2)>0)/length(criteria);
      
       
        disp(['using '   num2str(sum(criteria)) ' out of '  num2str(length(criteria)) ' events'])
    
        if 1  
        [Ensembles,EAinfo2,success] = subsetBoxes(Boxes, Ensembles,EAinfo,criteria,thres)

        end

%%%     load(['PhenoDrive2020_' num2str(i)], 'EvtsData') %analyzed
        if success
        temptimeerr = EAinfo2.EnTime(1) - Ensembles.Time(1);
        stretch = 20/abs(EAinfo2.deltaT2);
        eventTime(k)=abs(EAinfo2.deltaT2);
        TimeNew = (Ensembles.Time + temptimeerr)*stretch;
        
        BoxesList=fieldnames(Boxes);
        for j=1:length(BoxesList)
            temp=Ensembles.(BoxesList{j}); %SignalsT.Time,EvtSig,Boxes.Time,'nearest');
            tempinterpd = interp1(TimeNew,temp,Ensembles.Time);
            
            if 0 %debug plots
                figure(89); plot(EAinfo2.EnTime,temp); hold('on');
                plot(Ensembles.Time,tempinterpd); hold('on');
            end
            
            GroupEnsembles.(BoxesList{j})(k,:) = tempinterpd(:);
        end
        
      
    GGexponent=1;
    GGcalvalue(k) = GroupEnsembles.GGpeak(k,101-20).^GGexponent;
    GroupEnsembles.GGpeakPmax(k,:) =GroupEnsembles.GGpeak(k,:);
    GroupEnsembles.GGtonicPmax(k,:) = GroupEnsembles.GGtonic(k,:);
     GroupEnsembles.GGpeak(k,:) = 100*GroupEnsembles.GGpeakPmax(k,:).^GGexponent/GGcalvalue(k);
    GroupEnsembles.GGtonic(k,:) = 100*GroupEnsembles.GGtonicPmax(k,:).^GGexponent/GGcalvalue(k);
    SpO2calvalue(k) = max(GroupEnsembles.SpO2(k,101-20:101));
    GroupEnsembles.SpO2Calibrated(k,:) = 100*GroupEnsembles.SpO2(k,:)./SpO2calvalue(k)
   
    InclSubj(k)=i;
    k=k+1
        end
        
        
        
    catch
    end
    
end

EAinfoGroup = EAinfo2;
EAinfoGroup = struct;
EAinfoGroup.meanEvtTimes=eventTime;
EAinfoGroup.EnTime = [-100:100];
EAinfoGroup.InclSubj=InclSubj;
timebase=mean(EAinfoGroup.meanEvtTimes)
GroupEnsembleTime.Time=EAinfoGroup.EnTime/20*timebase
EAinfoGroup.PlotEnTime=EAinfoGroup.EnTime/20*timebase
EAinfoGroup.deltaT2 = 20;

settings.PlotEventData2Option=1;
EAinfoTemp = PlotEventData2(GroupEnsembles,GroupEnsembleTime,EAinfoGroup)

end


%%
% %%
% 
% settings.PlotEventData2Option=3;
% %%
% criteria= P==1
% timebase=mean(eventTime(criteria))
% %Ensembles2=Ensembles;
% Ensembles2.Time=Ensembles.Time/20*timebase
% 
% 
% EAinfoTemp = PlotEventData2(GroupEnsembles2,Ensembles2,EAinfoGroup,logical(criteria));
% ax=get(gcf,'children')
% temp=get(ax(3),'ytick');
% %scalefactor=nanmean(SpO2calvalue(criteria))
% %set(ax(3),'yticklabels',round(temp/100*scalefactor,2));
% saveas(15,'temp');
% uiopen('temp.fig',1);
% delete('temp.fig');
% 
% %%
% criteria= P==0; %drive early, GG very late -- compares well visually against group A
% timebase=mean(eventTime(criteria))
% %Ensembles2=Ensembles;
% Ensembles2.Time=Ensembles.Time/20*timebase
% 
% EAinfoTemp = PlotEventData2(GroupEnsembles2,Ensembles2,EAinfoGroup,logical(criteria));
% ax=get(gcf,'children')
% %temp=get(ax(3),'ytick');
% scalefactor=nanmean(SpO2calvalue(criteria))
% %set(ax(3),'yticklabels',round(temp/100*scalefactor,2));
% %%
% criteria= P==-1
% timebase=mean(eventTime(criteria))
% %Ensembles2=Ensembles;
% Ensembles2.Time=Ensembles.Time/20*timebase
% 
% EAinfoTemp = PlotEventData2(GroupEnsembles2,Ensembles2,EAinfoGroup,logical(criteria));
% ax=get(gcf,'children')
% %temp=get(ax(3),'ytick');
% scalefactor=nanmean(SpO2calvalue(criteria))
% %set(ax(3),'yticklabels',round(temp/100*scalefactor,2));
% 
% %%
% criteria= P==-1 | P==0 | P==-2
% timebase=mean(eventTime(criteria))
% Ensembles2=Ensembles;
% Ensembles2.Time=Ensembles.Time/20*timebase
% 
% EAinfoTemp = PlotEventData2(GroupEnsembles2,Ensembles2,EAinfoGroup,logical(criteria));
% ax=get(gcf,'children')
% temp=get(ax(3),'ytick');
% scalefactor=nanmean(SpO2calvalue(criteria))
% %set(ax(3),'yticklabels',round(temp/100*scalefactor,2));
% 
% 
% 
% %%
% criteria = ~isnan(P)
% timebase=mean(eventTime(criteria))
% %Ensembles2=Ensembles;
% Ensembles2.Time=Ensembles.Time/20*timebase
% 
% EAinfoTemp = PlotEventData2(GroupEnsembles2,Ensembles2,EAinfoGroup,logical(criteria));
% ax=get(gcf,'children')
% temp=get(ax(3),'ytick');
% scalefactor=nanmean(SpO2calvalue(criteria))
% %set(ax(3),'yticklabels',round(temp/100*scalefactor,2));