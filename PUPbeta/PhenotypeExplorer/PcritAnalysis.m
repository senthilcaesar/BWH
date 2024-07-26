    clear all;close all;
    [file,dir] = uigetfile('E:\Dropbox (Partners HealthCare)\PAtO\Traits\Converted');
        filedir = [dir, file];
        matObj = matfile(filedir);
        varlist = who(matObj);
        load(filedir);
        oldformat=0
        
    [~,subj,~] = fileparts(filedir);
    if strcmp(subj(end-2:end),'XHz')
        subj = subj(1:end-4); %strip off _XHz
    end
    if strcmp(subj(end-2:end),'rit')
        subj = subj(1:end-6); %strip off _XHz
    end
        
%     [~,~,raw]=xlsread('C:\Users\lg373\Dropbox (Partners HealthCare)\PhenotypeDrive2018\PUPstart\AMasterSpreadsheet','S4:T63')
%     raw = string(raw);
%     idx=(strcmp(subj,raw(:,2)))
%     subj = char(raw(idx,1))
        %%


%right click = finished completely
%click left of screen to move back Ymin
%click right of screen to move forward Ymin
%click left and right within screen to select a range to analyze later

if oldformat
Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1)); %
Flow=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)); %
Pmask=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pmask')==1)); %

Epochs=DataEventHypnog_Mat(:,strcmp(ChannelsList,'Epochs')==1); %
EventsAr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1)); 
EventsResp=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsResp')==1)); 
Edi=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Edi')==1)); 
WPr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'WPr')==1)); 
ArPr=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'ArPr')==1)); 
Position=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Position')==1)); 

signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','Edi','[WPr ArPr]'}
else
    Time=SigT.Time;
    Flow=SigT.Flow;
    Pmask=SigT.Pmask;   
    Edi=SigT.Edi;
    signallist={'[SigT.Epochs SigT.EventsAr 1+SigT.Position]','SigT.EventsResp','SigT.Pmask','SigT.Flow','SigT.Edi','[SigT.WPr SigT.ArPr]'};

end
dt=Time(2)-Time(1);

%signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','Edi','[WPr ArPr]'};
%signallist={'Pmask','Flow'};
%signallist={'[Epochs EventsAr 1+Position]','EventsResp','Pmask','Flow','[WPr ArPr]'};

global xvalues yvalues range ax2
clearvalues=1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range=5*60;

PlotAndSelectData(signallist,Time,clearvalues);

%%
Pcrit.runtimes=xvalues;

%% Analyze Pcrit

%1. %left click if flow is ok, or unsure; right click if zero flow,
%2. left click (inside plot) ok | right click (inside plot) if zero flow | ...
    % ...  right click left of screen to exclude | left click left of screen to go back 1

%Pcrit.runtimes=Pcritruntimes
figure(3);
set(gcf,'color',[1 1 1]);
ax3(1)=subplot(2,1,1); ylabel('Flow'); box('off');
ax3(2)=subplot(2,1,2); ylabel('Vol'); box('off');
i=1;
exclude=[];
Pcrit.VE=[];
Pcrit.Peakflow=[];
Pcrit.CPAP=[];
deltaX=0.5%1.2;

while i<=size(Pcrit.runtimes,1) 
    try
        delete(ax3(2));
    catch me1
    end
    
    %try
    
    lowerX=Pcrit.runtimes(i,1)-deltaX;
    upperX=Pcrit.runtimes(i,2)+deltaX;
    figure(3)
    I=find(Time>lowerX&Time<upperX);
    I2=find(Time>Pcrit.runtimes(i,1)&Time<Pcrit.runtimes(i,2));
    ax3(1)=subplot(2,1,1); plot(Time(I),Flow(I)+1,'k',[Pcrit.runtimes(i,1) Pcrit.runtimes(i,1)],[-1 1],'r',[Pcrit.runtimes(i,2) Pcrit.runtimes(i,2)],[-1 1],'r');
    ylabel('Flow'); box('off');
    %right click if zero flow
    [x1,~,button(i)] = ginput(1);
    if button(i)==3
        Pcrit.Peakflow_list{i}=0;
        Pcrit.Peakflow(i)=0;
        Pcrit.VE_list{i}=0;
        Pcrit.VE(i)=0;
        Pcrit.CPAP_list{i}=Pmask(find(Time>Pcrit.runtimes(i,1),1));
        Pcrit.CPAP(i)=Pcrit.CPAP_list{i};
        title(num2str(Pcrit.CPAP(i)));
        exclude(i)=0;
        Pcrit
        i=i+1;
        continue
    end
    vol1=cumsum(Flow(I)-mean(Flow(I)))*dt;
    vol2=cumsum(Flow(I2)-mean(Flow(I2)))*dt;
    thresh=std(vol2)/10;
    [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    %[indVI, indVE] = CheckVind(max_list, min_list);
    %TidalVol = V(indVI)-V(indVE);
    ax3(2)=subplot(2,1,2); plot(Time(I),vol1,'k',[Pcrit.runtimes(i,1) Pcrit.runtimes(i,1)],[-1 1],'r',[Pcrit.runtimes(i,2) Pcrit.runtimes(i,2)],[-1 1],'r');
    hold('on');
    
    plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'k.');
    plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'r.');
    %[ymin,temptemp]=get(ax2(4),'ylim')
    ylabel('Vol'); box('off');
    lefti=find(Time(I(1)-1+min_list(:,1))<Pcrit.runtimes(i,1),1,'last');
    righti=find(Time(I(1)-1+min_list(:,1))>Pcrit.runtimes(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)||righti==0
        righti=length(min_list)-1;
    end
    if size(min_list,1)>1
        ax3(2)=subplot(2,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[-1 1],'k:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[-1 1],'k:');
        
        xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
        ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
        ypoly = vol1(min_list(lefti:(righti+1),1));
        ppoly = polyfit(xpoly,ypoly,1);
        Yvoldrift = polyval(ppoly,Time(I));
        
        %repeat find end-exp and end-insp now that we have done some leak correction:
        vol3=cumsum(Flow(I)-mean(Flow(I))-ppoly(1))*dt;
        thresh=std(vol3)/10;
        [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
        
        ax3(2)=subplot(2,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'ko');
        ax3(2)=subplot(2,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'ro');
        ylabel('Vol'); box('off');
        lefti=find(Time(I(1)-1+min_list(:,1))<Pcrit.runtimes(i,1),1,'last');
        righti=find(Time(I(1)-1+min_list(:,1))>Pcrit.runtimes(i,2),1,'first')-1;
        if isempty(lefti)
            lefti=1;
        end
        if isempty(righti)||righti==0
            righti=size(min_list,1)-1;
        end
        ax3(2)=subplot(2,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[-1 1],'r:');
        
        xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
        ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
        ypoly = vol1(min_list(lefti:(righti+1),1));
        ppoly = polyfit(xpoly,ypoly,1);
        Yvoldrift = polyval(ppoly,Time(I));
        hold('off');
        ax3(1)=subplot(2,1,1); plot(Time(I),Flow(I)-mean(Flow(I))-ppoly(1),'k',[Pcrit.runtimes(i,1) Pcrit.runtimes(i,1)],[-1 1],'r',[Pcrit.runtimes(i,2) Pcrit.runtimes(i,2)],[-1 1],'r');
        
        Yvol1_exp = ypoly(1:end-1); %end exp
        Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
        Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
        Yvoldrift_insp = polyval(ppoly,Time(I(1)-1+max_list(lefti:(righti),1)));
        
        ax3(2)=subplot(2,1,2); hold('on'); plot(Time(I),Yvoldrift);
        endexpvol=Yvol1_exp-Yvoldrift_exp;
        endinspvol=Yvol1_insp-Yvoldrift_insp;
        volbreath=endinspvol-endexpvol;
        Flow1_corrected = Flow(I)-mean(Flow(I))-ppoly(1);
        for n=1:(length(xpoly)-1)
            I3=find(Time>xpoly(n)&Time<xpoly(n+1))-I(1);
            Pcrit.Peakflow_list{i}(n)=max(Flow1_corrected(I3));
        end
        Pcrit.Peakflow(i)=median(Pcrit.Peakflow_list{i});
        Pcrit.VE_list{i}=volbreath/ttot*60;
        Pcrit.VE(i)=median(volbreath)/ttot*60;
        Pcrit.CPAP_list{i}=Pmask(I(1)-1+min_list(lefti:(righti+1),1));
        Pcrit.CPAP(i)=median(Pmask(I(1)-1+min_list(lefti:(righti+1),1)));
        
    else
        volbreath=NaN;
        Flow1_corrected = Flow(I)-mean(Flow(I))-ppoly(1);
        Pcrit.Peakflow_list{i}=NaN;
        Pcrit.Peakflow(i)=median(Pcrit.Peakflow_list{i})
        Pcrit.VE_list{i}=volbreath/ttot*60;
        Pcrit.VE(i)=median(volbreath)/ttot*60;
        Pcrit.CPAP_list{i}=Pmask(I(1)-1+min_list(lefti));
        Pcrit.CPAP(i)=median(Pcrit.CPAP_list{i});
    end
    hold('off');
    
    %right click if zero flow | right click left of screen to exclude |
    [x1,~,button(i)] = ginput(1);
    
    if x1<Time(I(1))&&button(i)==3 %exclude
        exclude(i)=1;
    else
        exclude(i)=0;
    end
    
    if x1<Time(I(1))&&button(i)==1&&i>1 %go back
        i=i-2;
    end
    
   
    Pcrit
    i=i+1;
    
end

Pcrit.VE(button==3)=0;
Pcrit.Peakflow(button==3)=0;

Pcrit.VE(exclude==1|isnan(Pcrit.CPAP))=[];
Pcrit.Peakflow(exclude==1|isnan(Pcrit.CPAP))=[];
Pcrit.CPAP(exclude==1|isnan(Pcrit.CPAP))=[];


%% Pcrit analysis and plots

ImportPreviousManualSelection=0; %if loading/rerunning already-analyzed Pcrit structure variable data from file: 1
plotwithnumbers=0; %1 for review, 0 for presentation

if ~exist('Pcrit.minPmask')
Pcrit.minPmask = Inf;
end
%Manual settings as needed:
if ~ImportPreviousManualSelection
  %  Pcrit.minPmask= Inf;
    removelist1 = [1 2 3 4 5 6 7 8];%[25:38 1 2 3 9 10 17 18]; % remove outlier data for better fit quality: e.g.  [10 18 4 5 12 26] 
end

%%
%First run: set plotwithnumbers to 1, set removelist1 to [], set Pcrit.minPmask=Inf;
%Set Pcrit.minPmask to level needed for a linear fit (exclude high CPAP data that starts to be on plateau)
%Add outliers to removelist1
%Last run: set plotwithnumbers to 0;

fh=figure(12); clf(12);
global fixedslope upperlimit;
fixedslope=NaN;
upperlimit.on=0;
Pcrit.forcedslopes = [NaN NaN];

%Pcrit.forcedslopes = [3.6519 0.9120];
if ~ImportPreviousManualSelection %skip this if loaded older data
    Pcrit.include_data=1:length(Pcrit.CPAP); %1:length(Pcrit.CPAP) [20:45];
    removelist1 = [removelist1 find(Pcrit.CPAP>Pcrit.minPmask)];
    Pcrit.include_data(removelist1)=[];
else %use existing include_data
    removelist1 = find(any(Pcrit.include_data==(1:length(Pcrit.CPAP))',2)==0)';
end

Pmaskdata=Pcrit.CPAP(Pcrit.include_data);
VdotMaxdata=Pcrit.Peakflow(Pcrit.include_data)*60;
VdotEdata=Pcrit.VE(Pcrit.include_data);

set(fh,'color',[1 1 1]);
ax1(1)=subplot(2,1,1);
plotpoints=1;

if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(1);
end
%[VcritSEM,Vcritvalue,PcritSEM,Pcritvalue,PVSlope,Rsquared]=pcrit_model_run(Pmask,Vflow,exclude_zeroflowdata,plotpoints)
[Pcrit.VdotcritSEM,Pcrit.Vdotcritvalue,Pcrit.PcritSEM,Pcrit.Pcritvalue,Pcrit.PVSlope_Vmax]=pcrit_model_run(Pmaskdata,VdotMaxdata,0,plotpoints)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('Vpeak','fontname','arial narrow','fontsize',12);
ax1(2)=subplot(2,1,2);

if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(2);
end

[Pcrit.VcritSEM,Pcrit.Vcritvalue,Pcrit.PcritSEM_VE,Pcrit.Pcritvalue_VE,Pcrit.PVSlope_VE]=pcrit_model_run([Pmaskdata],[VdotEdata],0,plotpoints);
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
linkaxes(ax1,'x'); 

if plotwithnumbers
    for i=1:length(Pcrit.include_data)
        ax2(1)=subplot(2,1,2);
        h=text(Pmaskdata(i)+0.1,VdotEdata(i)+0.1,num2str(Pcrit.include_data(i)),'fontname','arial narrow','fontsize',9);
        ax2(1)=subplot(2,1,1);
        h=text(Pmaskdata(i)+0.1,VdotMaxdata(i)+0.1,num2str(Pcrit.include_data(i)),'fontname','arial narrow','fontsize',9);
    end
end

set(gcf,'position',[745    42   375   924])

%% Pcrit analysis and plots - median binned data
Pcrit.usealldataformedianbinned=0;

Ndata=[]; VdotMaxdata=[]; VdotEdata=[]; Pmaskdata=[]; 
VdotEdataU=[];
VdotEdataL=[];
VdotMaxdataU=[];
VdotMaxdataL=[];

if ~(isfield(Pcrit,'include_data2') && ~exist('removelist1')) %skip this if loaded older data
    Pcrit.include_data2=1:length(Pcrit.CPAP); 
    %removelist2 = []; 
    removelist2 = removelist1; 
    removelist2 = [removelist2 find(Pcrit.CPAP>Pcrit.minPmask)];
    Pcrit.include_data2(removelist2)=[];
end

if Pcrit.usealldataformedianbinned
Pmaskdata_=Pcrit.CPAP;
VdotMaxdata_=Pcrit.Peakflow*60;
VdotEdata_=Pcrit.VE;
else
Pmaskdata_=Pcrit.CPAP(Pcrit.include_data2);
VdotMaxdata_=Pcrit.Peakflow(Pcrit.include_data2)*60;
VdotEdata_=Pcrit.VE(Pcrit.include_data2);
end


binwidth=1;
CPAPbins1=-0.5+(floor(min(Pmaskdata_))-0):binwidth:(ceil(max(Pmaskdata_))-binwidth*0.5);

%CPAPbins1=-0.5+(floor(min(Pcrit.CPAP))-0):binwidth:(ceil(max(Pcrit.CPAP))-binwidth*0.5);

for i=1:length(CPAPbins1)
    tempI=find(Pmaskdata_>CPAPbins1(i)&Pmaskdata_<=(CPAPbins1(i)+binwidth));
    Ndata(i)=length(tempI);
    VdotMaxdata(i)=median(VdotMaxdata_(tempI));
    VdotMaxdataU(i)=prctile(VdotMaxdata_(tempI),75);
    VdotMaxdataL(i)=prctile(VdotMaxdata_(tempI),25);
    VdotEdata(i)=median(VdotEdata_(tempI));
    VdotEdataU(i)=prctile(VdotEdata_(tempI),75);
    VdotEdataL(i)=prctile(VdotEdata_(tempI),25);    
    Pmaskdata(i)=mean(Pmaskdata_(tempI));
end

plotwithnumbers=1;
fh=figure(13); clf(13);
fixedslope=NaN;

I2 = find(Pmaskdata>Pcrit.minPmask|isnan(Pmaskdata));
    Ndata(I2)=[];
    VdotMaxdata(I2)=[];
    VdotEdata(I2)=[];
    Pmaskdata(I2)=[];
    VdotEdataU(I2)=[];
    VdotEdataL(I2)=[];
    VdotMaxdataU(I2)=[];
    VdotMaxdataL(I2)=[];

set(fh,'color',[1 1 1]);
ax1(1)=subplot(2,1,1);
plotpoints=1;


if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(1);
end

[Pcrit.VdotcritSEM_medianmethod,Pcrit.Vdotcritvalue_medianmethod,Pcrit.PcritSEM_medianmethod,Pcrit.Pcritvalue_medianmethod,Pcrit.PVSlope_Vmax_medianmethod]=pcrit_model_run(Pmaskdata,VdotMaxdata,0,plotpoints)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('Vpeak','fontname','arial narrow','fontsize',12);
for i=1:length(Ndata)
    plot([Pmaskdata(i) Pmaskdata(i)],[VdotMaxdataL(i) VdotMaxdataU(i)],'k');
end

if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(2);
end

ax1(2)=subplot(2,1,2);
[Pcrit.VcritSEM_medianmethod,Pcrit.Vcritvalue_medianmethod,~,~,Pcrit.PVSlope_VE_medianmethod]=pcrit_model_run([Pmaskdata],[VdotEdata],0,plotpoints);
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
linkaxes(ax1,'x'); 
for i=1:length(Ndata)
    plot([Pmaskdata(i) Pmaskdata(i)],[VdotEdataL(i) VdotEdataU(i)],'k');
end

set(gcf,'position',[745    42   375   924])


%% Save
%save(filenameanddir,'Pcrit','-append');

idcs   = strfind(dir,'\');
 newdir = dir(1:idcs(end-1)-1);
 if ~(exist([newdir '\CPAPTraits'], 'dir') == 7)
        mkdir([newdir '\CPAPTraits']);
 end

 % newdir = 'C:\Users\SOpdeBeeck\OneDrive - uantwerpen\Pcrit feasibility\Data\Converted\'
  % newdir = 'C:\Users\sarao\OneDrive - uantwerpen\Pcrit feasibility\Data\Converted\'
 % subj = 'PcritUZA013'
 
save([newdir '\CPAPTraits' '\' subj '_Pcrit'],'Pcrit');
saveas(12,[newdir '\CPAPTraits' '\' subj '_PcritPaper'],'fig');
saveas(13,[newdir '\CPAPTraits' '\' subj '_PcritMedianPaper'],'fig');
saveas(12,[newdir '\CPAPTraits' '\' subj '_PcritPaper'],'png');
saveas(13,[newdir '\CPAPTraits' '\' subj '_PcritMedianPaper'],'png');
