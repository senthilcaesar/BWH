%% Find ArthresEdi
oldformat=0

%clear all;close all
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
       
    idcs   = strfind(dir,'\');
    newdir = dir(1:idcs(end-1)-1);
    
    %%
load([newdir '\CPAPTraits\' subj '_Pcrit.mat'])
load([newdir '\CPAPTraits\' subj '_Veupnea.mat'])

load([newdir '\CPAPTraits\' subj '_Vactive.mat'])
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


Thor=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Thorax')==1)); %EEG selection uses scoring
Abdo=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Abdomen')==1)); %EEG selection uses scoring
Thor=Thor-mean(Thor);
Abdo=Abdo-mean(Abdo);

signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','Edi','[WPr ArPr]'}
else 
    Time=SigT.Time;
    Flow=SigT.Flow;
    Pmask=SigT.Pmask;
    Edi=SigT.Edi;
    signallist={'[SigT.Epochs SigT.EventsAr 1+SigT.Position]','SigT.EventsResp','SigT.Pmask','SigT.Flow','SigT.Edi','[SigT.WPr SigT.ArPr]'};
end



dt=Time(2)-Time(1);

%%
global xvalues yvalues range ax2
clearvalues=0;
range=5*60;
xvalues=[];
yvalues=[];
% It = Time(find(diff(EventsAr)==1));
PlotAndSelectData(signallist,Time,clearvalues);

%%

Vactive.ranges=xvalues;
%xvalues=Vactive.ranges;
%xvalues=Vactive.ranges+Time(1);

%Vactive.ranges=Vactive.ranges+Time(1);
%yvalues=0*xvalues;
%Vactive.ranges=xvalues(I,:);

%% Analyze Vactive v2

%right click if zero flow
%left click otherwise
%there is no way to exclude data from here as yet
%Vactive.ranges=Pcritruntimes
deltaX=0.5;%2.5;
Vactive.Vactive_list=[];
Vactive.CPAP_list=[];
Vactive.State=[];
Vactive.Drive=[];
Vactive.DriveDelta=[];

i=1;
exclude=[];
clear button
figure(3);
Vactive.Vactive_list=[];

Drive = Edi;


%Flow=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %EEG selection uses scoring
%     filter_HFcutoff_butter0 = 6;
%     filter_order0 = 2;
%     [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
%     Flow = nanfilter(Flow,B_butter0,A_butter0,1); %filtfilt, otherwise flow signal is right-shifted

    Fthres = 1/3;
    
while i<=size(Vactive.ranges,1)
        close(3)
        try
        figure(3)
        try
        if i>1
        subplot(3,1,1); plot(Time(I),Flow(I)+1,'k',[Vactive.ranges(i,1) Vactive.ranges(i,1)],[min(Flow(I)) max(Flow(I))],'r',[Vactive.ranges(i,2) Vactive.ranges(i,2)],[min(Flow(I)) max(Flow(I))],'r');
        subplot(3,1,2); plot(Time(I),vol1,'k',[Vactive.ranges(i,1) Vactive.ranges(i,1)],[min(vol1(I)) max(vol1(I))],'r',[Vactive.ranges(i,2) Vactive.ranges(i,2)],[min(vol1(I)) max(vol1(I))],'r');
        subplot(3,1,3); plot(Time(I),Drive(I),'k');
        title(i)'
        end
        catch me1
        end
        
        lowerX=Vactive.ranges(i,1)-deltaX;
        upperX=Vactive.ranges(i,2)+deltaX;
        
        I=find(Time>lowerX&Time<upperX);
        I2=find(Time>Vactive.ranges(i,1)&Time<Vactive.ranges(i,2));
        subplot(3,1,1); plot(Time(I),Flow(I)+1,'k',[Vactive.ranges(i,1) Vactive.ranges(i,1)],[min(Flow(I)) max(Flow(I))],'r',[Vactive.ranges(i,2) Vactive.ranges(i,2)],[min(Flow(I)) max(Flow(I))],'r');
        vol1=cumsum(Flow(I)-mean(Flow(I)))*dt;
        vol2=cumsum(Flow(I2)-mean(Flow(I2)))*dt;
        thresh=std(vol2)*Fthres;
        [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
        %[indVI, indVE] = CheckVind(max_list, min_list);
        %TidalVol = V(indVI)-V(indVE);
        subplot(3,1,2); plot(Time(I),vol1,'k',[Vactive.ranges(i,1) Vactive.ranges(i,1)],[min(vol1) max(vol1)],'r',[Vactive.ranges(i,2) Vactive.ranges(i,2)],[min(vol1) max(vol1)],'r');
        hold('on')
        subplot(3,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'k.');
        subplot(3,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'r.');
        %[ymin,temptemp]=get(ax2(4),'ylim')
        
        subplot(3,1,3); plot(Time(I),Drive(I),'k');
        
        lefti=find(Time(I(1)-1+min_list(:,1))<(Vactive.ranges(i,1)+0.5),1,'last');
        righti=find(Time(I(1)-1+min_list(:,1))>Vactive.ranges(i,2),1,'first')-1;
        if isempty(lefti)
            lefti=1;
        end
        if isempty(righti)||righti==0
            righti=length(min_list)-1;
        end
        if size(min_list,1)>1
            subplot(3,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[min(vol1) max(vol1)],'k:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[min(vol1) max(vol1)],'k:');
            
            xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
            ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
            ypoly = vol1(min_list(lefti:(righti+1),1));
            ppoly = polyfit(xpoly,ypoly,1);
            Yvoldrift = polyval(ppoly,Time(I));
            
            %repeat find end-exp and end-insp now that we have done some leak correction:
            vol3=cumsum(Flow(I)-mean(Flow(I))-ppoly(1))*dt;
            thresh=std(vol3)*Fthres;
            [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
            
            subplot(3,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'ko');
            subplot(3,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'ro');
            
            lefti=find(Time(I(1)-1+min_list(:,1))<Vactive.ranges(i,1),1,'last');
            righti=find(Time(I(1)-1+min_list(:,1))>Vactive.ranges(i,2),1,'first')-1;
            if isempty(lefti)
               lefti=1;
               
            end
            if isempty(righti)||righti==0
                righti=size(min_list,1)-1;
            end
            subplot(3,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[min(vol1) max(vol1)],'r:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[min(vol1) max(vol1)],'r:');
            
            xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
            ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
            ypoly = vol1(min_list(lefti:(righti+1),1));
            ppoly = polyfit(xpoly,ypoly,1);
            Yvoldrift = polyval(ppoly,Time(I));
            
            subplot(3,1,1); plot(Time(I),Flow(I)-mean(Flow(I))-ppoly(1),'k',[Vactive.ranges(i,1) Vactive.ranges(i,1)],[min(Flow(I)-mean(Flow(I))-ppoly(1)) max(Flow(I)-mean(Flow(I))-ppoly(1))],'r',[Vactive.ranges(i,2) Vactive.ranges(i,2)],[min(Flow(I)-mean(Flow(I))-ppoly(1)) max(Flow(I)-mean(Flow(I))-ppoly(1))],'r');
            
            Yvol1_exp = ypoly(1:end-1); %end exp
            Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
            Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
            Yvoldrift_insp = polyval(ppoly,Time(I(1)-1+max_list(lefti:(righti),1)));
            
            subplot(3,1,2); plot(Time(I),Yvoldrift);
            endexpvol=Yvol1_exp-Yvoldrift_exp;
            endinspvol=Yvol1_insp-Yvoldrift_insp;
            volbreath=endinspvol-endexpvol;
            %     Flow1_corrected = Flow(I)-mean(Flow(I))-ppoly(1);
            %     for n=1:(length(xpoly)-1);
            %         I3=find(Time>xpoly(n)&Time<xpoly(n+1))-I(1);
            %         Pcrit.Peakflow_list{i}(n)=max(Flow1_corrected(I3))
            %     end
            Vactive.Vactive_list(i)=mean(volbreath)/ttot*60;
            
            DriveTemp = Drive(I);
            DriveTemp(min_list(:,1));
            clear DriveMaxB DriveMaxBi DriveMinB
            for ii=1:size(min_list,1)-1
                irange = min_list(ii,1):min_list(ii+1,1);
                [DriveMaxB(ii,1),idx] = max(DriveTemp(irange))
                DriveMaxBi(ii,1) = idx + irange(1) - 1;
            end
            for ii=1:size(min_list,1)-1
                if ii==1
                    li = 1;
                else
                    li = DriveMaxBi(ii-1);
                end
                irange2 = li:DriveMaxBi(ii);
                [DriveMinB(ii,1)] = min(DriveTemp(irange2));
            end
            Vactive.Drive(i) = mean(DriveMaxB);
            Vactive.DriveDelta(i) = mean(DriveMaxB-DriveMinB);
            
            subplot(3,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'ro');
        else
            volbreath=NaN;
            Vactive.Vactive_list(i)=mean(volbreath)/ttot*60;
            
        end
        Vactive.CPAP_list(i)=mean(Pmask(I(1)-1+min_list(lefti:(righti+1),1))); %zero flow values
        Vactive.State(i)=Epochs(I2(1));
        hold('off');
        
        Vactive.Vactive_list(i)=mean(volbreath)/ttot*60;
        %Vactive.CPAP_list(i)=mean(Pmask(I(1)-1+min_list(lefti:(righti+1),1)));
        
    catch me
    end
    
    disp('right click if zero flow | right click left of screen to exclude');
    [x1,~,button(i)] = ginput(1);
    
    if x1<Time(I(1))&&button(i)==3 %exclude
        exclude(i)=1;
        disp('excluding selected');
    else
        exclude(i)=0;
    end
    
    if x1<Time(I(1))&&button(i)==1&&i>1 %go back
        i=i-2;
    end
    i=i+1;
end

%%
%button(length(Vactive.Vactive_list)+1:end)=[];
Vactive.Vactive_list(button==3)=0;

I = button==3 & Vactive.CPAP_list==0;
Vactive.CPAP_list(I)=Pmask(round((Vactive.ranges(I,1)-Time(1))/dt)+1); %fix to replace zeros in Pmask when zero flow is selected.

I=exclude==1;
Vactive.Vactive_list(I)=[];
Vactive.CPAP_list(I)=[];
Vactive.State(I)=[];
Vactive.Drive(I)=[];
Vactive.DriveDelta(I)=[];


%% Vactive analysis and plots

ImportPreviousManualSelection=0; %ImportPreviousManualSelection
plotwithnumbers=0;
removezerozeros=1; %remove erroneous double 0 entries
Vactive.minPmask = Inf;
%Manual settings as needed:
if ~ImportPreviousManualSelection
  %  Vactive.minPmask = Inf;
    removelist1 = []; % remove outlier data for better fit quality: e.g.  [10 18 4 5 12 26] 
end

%% Draw active and passive pressure-VE plots
global fixedslope upperlimit
upperlimit.on=0;
figure(14); clf(14);
set(14,'color',[1 1 1]);
plotpoints=0;
ax2(1)=subplot(1,1,1);


if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(2);
else
    fixedslope=NaN;
end

if 1
    Pmaskdata=Pcrit.CPAP(Pcrit.include_data);
    VdotEdata=Pcrit.VE(Pcrit.include_data);
end
[~,Vcritvalue,~,Pcritvalue,PVSlope]=pcrit_model_run(Pmaskdata,VdotEdata,0,plotpoints) %Pmaskdata VdotMaxdata VdotEdata
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);


            temp=get(ax2(1),'children');
            for i=1:length(temp)-1
                set(temp(i),'Color',[0 0 0]);
            end
            set(temp(end),'FaceAlpha',0);

%Start plotting Vactive points



if removezerozeros
    I0 = find(Vactive.CPAP_list==0 & Vactive.Vactive_list==0);
    Vactive.ranges(I0,:)=[];
    Vactive.CPAP_list(I0)=[];
    Vactive.Vactive_list(I0)=[];
    Vactive.State(I0)=[];
end

%Pcrit.forcedslopes = [3.6519 0.9120];
if ~ImportPreviousManualSelection || ~isfield(Vactive,'include_data') %skip this if loaded older data
    Vactive.include_data=1:length(Vactive.CPAP_list); %1:length(Pcrit.CPAP) [20:45];
    removelist1 = [removelist1 find(Vactive.CPAP_list>Vactive.minPmask)];
    Vactive.include_data(removelist1)=[];
else %use existing include_data
    removelist1 = find(any(Vactive.include_data==(1:length(Vactive.CPAP_list))',2)==0)';
end

Pdata=Vactive.CPAP_list(Vactive.include_data);
Vdata=Vactive.Vactive_list(Vactive.include_data);
% Vdata(Pdata>maxPdata)=[];
% Pdata(Pdata>maxPdata)=[];

fixedslope=Pcrit.PVSlope_VE;
%fixedslope=NaN;
[VcritvalueSEMA,VcritvalueA,PcritvalueSEMA,PcritvalueA,PVSlopeA]=pcrit_model_run(Pdata,Vdata,0,1);

xlim([min([0 Vactive.CPAP_list(Vactive.include_data) Pcrit.CPAP(Pcrit.include_data)])-0.2 max([0 Vactive.CPAP_list(Vactive.include_data) Pcrit.CPAP(Pcrit.include_data)])+0.2]);

if plotwithnumbers
    for i=1:length(Vactive.include_data)
        h=text(Pdata(i)+0.1,Vdata(i)+0.1,num2str(Vactive.include_data(i)),'fontname','arial narrow','fontsize',9);
    end
end

set(gcf,'position',[1120          96         374         440]);

Vactive.VactiveSampleAbovePassiveLine = Vactive.Vactive_list - (Pcrit.PVSlope_VE.*Vactive.CPAP_list+Vcritvalue);
Vactive.VactiveSampleAboveActiveLine = Vactive.Vactive_list - (Pcrit.PVSlope_VE.*Vactive.CPAP_list+VcritvalueA);

%% Drive Effects
figure(15); clf(15); set(gcf,'color',[1 1 1]);

fixedslope=NaN; %fixed slope zero not working
upperlimit.on=0;
TempDrive = Vactive.Drive(Vactive.include_data);
    DriveMean = nanmean(TempDrive); TempDrive = TempDrive-nanmean(TempDrive);
[~,Vcomp_atMeanDrive,~,~,Vcomp_Slope,~]=pcrit_model_run(TempDrive,Vactive.VactiveSampleAbovePassiveLine(Vactive.include_data),0,1,1);

% figure(16); clf(16); set(gcf,'color',[1 1 1]);
% plot(Vactive.Drive(Vactive.include_data),Vactive.VactiveSampleAbovePassiveLine(Vactive.include_data),'k.','markersize',12);
% [~,temp,~,~,Vcomp_Slope,~]=pcrit_model_run(Vactive.Drive(Vactive.include_data),Vactive.VactiveSampleAbovePassiveLine(Vactive.include_data),0,1,1);

if plotwithnumbers
    for i=1:length(Vactive.include_data)
        h=text(TempDrive(i)+0.1,Vactive.VactiveSampleAbovePassiveLine(i)+0.1,num2str(Vactive.include_data(i)),'fontname','arial narrow','fontsize',9);
    end
end

temp = get(gca,'xtick');
set(gca,'xticklabels',round(temp+DriveMean,1));
% % 
% temp = get(gca,'ytick');
% set(gca,'yticklabels',round(temp+VcritvalueA-Vcritvalue,2));

box('off');
ylabel('Vcomp, L/min');
xlabel('Drive');

set(gcf,'position',[  754   523   366   425]);

%% Summary
% Vcomp data are correct since 5pm 2020-05-18; load rerun and resave older data

%load Pcrit and Veupnea
Vpassive.Feupnea = Pcrit.Vcritvalue/Veupnea.mean;
Vpassive.Direct_Feupnea_median = median(Pcrit.VE(Pcrit.CPAP>-1&Pcrit.CPAP<1))/Veupnea.mean;
Vpassive.Direct_Feupnea_N = length(Pcrit.CPAP>-1&Pcrit.CPAP<1);

%load Veupnea:
Vactive.Vcrit_Feupnea = VcritvalueA/Veupnea.mean;
Vactive.Vcrit_Feupnea_SEM = VcritvalueSEMA/Veupnea.mean;
Vactive.Active_Pcrit = PcritvalueA;
Vactive.Active_Pcrit_SEM = PcritvalueSEMA;
Vactive.Direct_Feupnea_mean = mean(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1))/Veupnea.mean;
Vactive.Direct_Feupnea_median = median(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1))/Veupnea.mean;
Vactive.Direct_Feupnea_N = length(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1));
Vactive.Direct_Feupnea_SEM = std(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1))/(Vactive.Direct_Feupnea_N^0.5)/Veupnea.mean;
Vactive.DriveMean=DriveMean;
Vactive.Vcomp_atMeanDrive=Vcomp_atMeanDrive;
Vactive.Vcomp_Slope=Vcomp_Slope;

%% Save
%save(filenameanddir,'Pcrit','-append');


 if ~(exist([newdir '\CPAPTraits'], 'dir') == 7)
        mkdir([newdir '\CPAPTraits']);
 end
    
save([newdir '\CPAPTraits' '\' subj '_Vactive'],'Vactive');
saveas(14,[newdir '\CPAPTraits' '\' subj '_Vactive'],'fig');
saveas(14,[newdir '\CPAPTraits' '\' subj '_Vactive'],'png');
saveas(15,[newdir '\CPAPTraits' '\' subj '_VcompVsDrive'],'fig');
saveas(15,[newdir '\CPAPTraits' '\' subj '_VcompVsDrive'],'png');
save([newdir '\CPAPTraits' '\' subj '_Vpassive'],'Vpassive');

