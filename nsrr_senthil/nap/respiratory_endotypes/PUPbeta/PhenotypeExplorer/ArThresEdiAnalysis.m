%% Find ArthresEdi (not ready)


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
        
        %%

%select periods of normal eupneic ventilation on CPAP during NREM sleep without arousals
%avoid flow limitation, post-arousal overshoot/undershoot, REM


Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1)); %EEG selection uses scoring
Flow=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %EEG selection uses scoring
Pmask=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pmask')==1)); %EEG selection uses scoring
Epochs=DataEventHypnog_Mat(:,strcmp(ChannelsList,'Epochs')==1); %EEG selection uses scoring
EventsAr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1)); %EEG selection uses scoring
EventsResp=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsResp')==1)); %EEG selection uses scoring
Edi=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Edi')==1)); %EEG selection uses scoring
WPr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'WPr')==1)); %EEG selection uses scoring
ArPr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'ArPr')==1)); %EEG selection uses scoring

Position=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Position')==1)); %EEG selection uses scoring

dt=Time(2)-Time(1);
signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','Edi','[WPr ArPr]'};

global xvalues yvalues range ax2
clearvalues=0;
range=5*60;
It = Time(find(diff(EventsAr)==1));
PlotAndSelectData(signallist,Time,clearvalues);

%%

%Right-click to exclude and move on
%Click left and right WITHIN the 5 breaths prior to arousal to include these breaths

clear VARx1 VARy1
figure(2)
Xminute=3;
Yminute=1;
dsf=5;
X=4;


ax2(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),'b',downsample(Time,dsf),downsample(EventsArXHz,dsf),'r');
if exist('WakeSleep'), hold('on'); plot(downsample(Time,dsf),downsample(WakeSleep,dsf),'color',[1 0.3 0.3]); hold('off'); end
ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Flow,dsf)); ylim([-1.2 1.2]);
%ax2(5)=subplot(X,1,5); plot(downsample(Time,dsf),downsample(Edi,dsf));
ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(FlowEdi,dsf)); ylim([-1.2 1.2]);
linkaxes(ax2,'x')

%arousals are at temp==1:
temp=diff(EventsArXHz);
I=find(temp==1);
i=1;
while i<=length(I)
    lowerX=Time(I(i))-Yminute*60;
    upperX=Time(I(i))+Yminute*60;
    xlim([lowerX upperX])
    
    lefti=(lowerX-Time(1))/dt+1;
    righti=(upperX-Time(1))/dt+1;
    
    xdata.data = Edi(lefti:righti); %downsample(Pes(I),downsamplefactor);
    xdata.time = Time(lefti:righti);
    xdata.dt = dt;
    
    linkaxes(ax2,'x')
    %[ymin,temptemp]=get(ax2(4),'ylim')
    %set(ax2(4),'ylim',[min([Pepi1trend(Time>lowerX&Time<upperX)+ArThresEpi.mean]-5) prctile(Pepi(Time>lowerX&Time<upperX),99.5)]);
    [VARx1(i,1),~,button(i)]=ginput(1);
    %     if button(i)~=3
    %         [VARx1(i,2),VARy1(i,2),button(i)]=ginput(1);
    %     end
    i=i+1;
end

VARx1(button==3)=[]; %keep 'approved' x values for left and right of flow trace.
VARx1 = [VARx1 VARx1+0.1];
ArthresEdi.ranges=VARx1;
ArthresEdi.button = button;
%% --save
if saveontherun
    save(filenameanddir,'ArthresEdi','-append');
end

%% Analyze ArThres using Edi and model
%load(filenameanddir,'ArthresEdi');
button=[];
chooseArs=0;
%ArthresEdi.comments = 'variable calibration';
ArthresEdi.Frescale = 1;
for i=1:size(ArthresEdi.ranges,1)
    
    deltaX=8;
    lowerX=ArthresEdi.ranges(i,1)-deltaX;
    upperX=ArthresEdi.ranges(i,2)+deltaX;
    lefti=(lowerX-Time(1))/dt+1;
    righti=(upperX-Time(1))/dt+1;
    
    modelY = ArthresEdi.Frescale*FlowEdi(lefti:righti); %downsample(Pes(I),downsamplefactor);
    
    xdata.time = Time(lefti:righti);
    xdata.dt = dt;
%     if 0
%         modelY = EdiToVflowModel(EdiToFlow.Parameters(end,:),xdata); % parameters_x(end,:)
%     else
%         modelY = FlowEdi(lefti:righti);   
%     end
    
    figure(3)
    %I=find(Time>lowerX&Time<upperX);
    subplot(2,1,1); plot(xdata.time,modelY,'k',[ArthresEdi.ranges(i,1) ArthresEdi.ranges(i,1)],[-1 1],'r',[ArthresEdi.ranges(i,2) ArthresEdi.ranges(i,2)],[-1 1],'r');
    vol1=cumsum(modelY-mean(modelY))*dt;
    thresh=std(vol1)/10;
    [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    %[indVI, indVE] = CheckVind(max_list, min_list);
    %TidalVol = V(indVI)-V(indVE);
    subplot(2,1,2); plot(xdata.time,vol1,'k',[ArthresEdi.ranges(i,1) ArthresEdi.ranges(i,1)],[-1 1],'r',[ArthresEdi.ranges(i,2) ArthresEdi.ranges(i,2)],[-1 1],'r');
    hold('on')
    subplot(2,1,2); plot(xdata.time(max_list(:,1)),vol1(max_list(:,1)),'k.');
    subplot(2,1,2); plot(xdata.time(min_list(:,1)),vol1(min_list(:,1)),'r.');
    %[ymin,temptemp]=get(ax2(4),'ylim')
    
    lefti=find(xdata.time(min_list(:,1))<ArthresEdi.ranges(i,1),1,'last');
    righti=find(xdata.time(min_list(:,1))>ArthresEdi.ranges(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)
        righti=lefti;
    end
    if righti<lefti
        righti=lefti;
    end
    subplot(2,1,2); plot([xdata.time(min_list(lefti,1)) xdata.time(min_list(lefti,1))],[-1 1],'r:',[xdata.time(min_list(righti+1,1)) xdata.time(min_list(righti+1,1))],[-1 1],'r:');
    
    xpoly = xdata.time(min_list(lefti:(righti+1),1));
    ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
    ypoly = vol1(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,xdata.time);
    
    %repeat find end-exp and end-insp now that we have done some leak correction:
    vol3=cumsum(modelY-mean(modelY)-ppoly(1))*dt;
    thresh=std(vol3)/10;
    [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    
    subplot(2,1,2); plot(xdata.time(max_list(:,1)),vol1(max_list(:,1)),'ko');
    subplot(2,1,2); plot(xdata.time(min_list(:,1)),vol1(min_list(:,1)),'ro');
    
    lefti=find(xdata.time(min_list(:,1))<ArthresEdi.ranges(i,1),1,'last');
    righti=find(xdata.time(min_list(:,1))>ArthresEdi.ranges(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)||righti==0
        righti=size(min_list,1)-1;
    end
    subplot(2,1,2); plot([xdata.time(min_list(lefti,1)) xdata.time(min_list(lefti,1))],[-1 1],'r:',[xdata.time(min_list(righti+1,1)) xdata.time(min_list(righti+1,1))],[-1 1],'r:');
    
    xpoly = xdata.time(min_list(lefti:(righti+1),1));
    ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
    ypoly = vol1(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,xdata.time);
    
    Yvol1_exp = ypoly(1:end-1); %end exp
    Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
    Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
    Yvoldrift_insp = polyval(ppoly,xdata.time(max_list(lefti:(righti),1)));
    
    subplot(2,1,2); plot(xdata.time,Yvoldrift);
    hold('off');
    
    endexpvol=Yvol1_exp-Yvoldrift_exp;
    endinspvol=Yvol1_insp-Yvoldrift_insp;
    volbreath=endinspvol-endexpvol;
    ArthresEdi.data(i)=mean(volbreath)/ttot*60;
    ArthresEdi.CPAP(i)=mean(Pmask(I(1)-1+min_list(lefti:(righti+1),1)));
    if chooseArs
    [~,~,button(i)]=ginput(1);
    else
        button(i)=1;
    end
end

II = ArthresEdi.CPAP>-1&ArthresEdi.CPAP<1;

ArthresEdi.mean = mean(ArthresEdi.data(button==1));
ArthresEdi.median = median(ArthresEdi.data(button==1));
ArthresEdi.Feupnea = ArthresEdi.median/Veupnea.mean;
ArthresEdi.F_sem = std(ArthresEdi.data(button==1))/length(button==3)^0.5/Veupnea.mean;

%%








%%
ArThresPes.data=[];
Pesnadirabs=[];
ArThresPes.CPAP=[]; %note this is the Pmask measured at the Pes nadir rather than during zero flow.
ArThresPes.data_t=[];
ArThresPes.data_i=[];

for i=1:length(ArThresPes.x_select)
    [Pesnadir_uncorrected,index]=min(Pepi(Time>ArThresPes.x_select(i)-1&Time<ArThresPes.x_select(i)+1));
    ArThresPes.data_t(i)=ArThresPes.x_select(i)-1+(index-1)*dt;
    ArThresPes.data_i(i)=find(Time>ArThresPes.data_t(i),1,'first');
    %plot
    lowerX=ArThresPes.x_select(i)-15;
    upperX=ArThresPes.x_select(i)+15;
    ArThresPes.data(i)=Pesnadir_uncorrected-ArThresPes.y_select(i);
    %ArThresPes.CPAP(i)=Pmask(find(Time>ArThresPes.data_t(i),1));
    %Pesnadirabs(i)=ArThresPes.data(i)+Pmask(ArThresPes.data_i(i));
    %plot
    %ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),downsample(Time,dsf),downsample(Time,dsf)*0+estimatedzeroPepi,'k:',ArThresEpi.data_t(i),Pepi(ArThresEpi.data_i(i)),'r.');
    ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),ArThresPes.data_t(i),Pepi(ArThresPes.data_i(i)),'r.');
    xlim([lowerX upperX]);
    set(ax2(3),'ylim',[min(Pepi(Time>lowerX&Time<upperX)) prctile(Pepi(Time>lowerX&Time<upperX),95)]);
    hold('off');
end

%% Summary analysis
ArThresPes.state=NaN*ArThresPes.data_i;
for i=1:length(ArThresPes.data_i)
    try
        ArThresPes.state(i)=max([EpochsXHz(ArThresPes.data_i(i)-15/dt) EpochsXHz(ArThresPes.data_i(i)-0/dt)]);
    catch me
    end
end
criteria = (ArThresPes.state<5);
ArThresPes.mean = nanmean(ArThresPes.data(criteria));
ArThresPes.median = nanmedian(ArThresPes.data(criteria));
ArThresPes.N = sum(~isnan(ArThresPes.data(criteria)));
ArThresPes.SEM = nanstd(ArThresPes.data(criteria))/ArThresPes.N^0.5;

ArThresPes.mean_abs = nanmean(Pesnadirabs(criteria));
ArThresPes.median_abs = nanmedian(Pesnadirabs(criteria));
ArThresPes.SEM_abs = nanstd(Pesnadirabs(criteria))/ArThresPes.N^0.5;

ArThresPes.cov = nanstd(ArThresPes.data(criteria))/-ArThresPes.mean;

ArThresPes.onCPAP = nanmean(ArThresPes.data(criteria&(ArThresPes.CPAP<-1|ArThresPes.CPAP>1)));
ArThresPes.offCPAP = nanmean(ArThresPes.data(criteria&ArThresPes.CPAP>=-1&ArThresPes.CPAP<=1));
ArThresPes.offCPAPN = sum(criteria&~isnan(ArThresPes.data)&ArThresPes.CPAP>=-1&ArThresPes.CPAP<=1);

%% --save
save(filenameanddir,'ArThresPes','-append');
