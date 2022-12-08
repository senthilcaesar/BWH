%% clear workspace 
clear all
close all

%% Load matlab data 
prefix = 'C:\Users\szs88\Dropbox (Personal)\';
prefix = 'C:\Users\szs88\Dropbox (Partners HealthCare)\';
savename = 'OralApplianceInternational\OA3_merged'; %change to OA3_AB from OA3_SS

%prefix = 'C:\Users\SIU\Dropbox\';
%savename = 'OralApplianceInternational\Analyzed\OA3_merged'; %change to OA3_AB from OA3_SS
load([prefix savename],'SleepData','MiscEpochData','Ratings','AnalysisIndex','LGplusinfo','EventsInfo','LG_QualityInfo','DataOut','ArousalDat','fitQual','StoNData','AHIdata','CPAPData')
M=length(LGplusinfo);


%%
clear AHItotalX
for i=1:M
    AHItotalX(i,1)=AHIdata{i}(58);
end

%% Import clinical variables from Excel
filexls=[prefix 'OralApplianceInternational\PUPbeta 20161117\AnalyzeDataSpreadsheet.xlsx'];
row1=3; 
row2=row1+M;
col1='L';
col2='AM';
sheet=1;
loadnewvariables={'BaselineAHI','DeltaAHI','DeltaAHIp','BMI','Age','Neck','TreatmentAHI','Sex','WC','ProtrusionMax','ProtrusionFinalMM','ProtrusionFinalPmax','BaselineAI','TreatmentAI','BaselineMinSpO2','TreatmentMinSpO2','BaselineAHInrem','TreatmentAHInrem','BaselineAHIrem','TreatmentAHIrem'};
colsnewvariables=[ 1     	     4          3           7     6     10     2              5     11      12          13                  14                        15                16               17               18              19             20                  21              22                      ];

range=[col1 num2str(row1) ':' col2 num2str(row2)];
[num,txt,raw]=xlsread(filexls,sheet,range);
for i=1:length(loadnewvariables)
    eval(['clear ' loadnewvariables{i} ';'])
    col=colsnewvariables(i);    %Includes flow events
    for j=1:size(raw,1)
        try
            eval([loadnewvariables{i} '(j,1)=raw{j,col};']);
        catch me
            eval([loadnewvariables{i} '(j,1)=NaN;']);
        end
    end
    if 1 %trimming data to length of available PUPdata
        eval([loadnewvariables{i} '((M+1):end)=[];'])
    end
end

Exclude=[BaselineAHI<20];

%% Check AHI

TableX=[BaselineAHI AHItotalX Exclude]

plot(AHItotalX(65:88),BaselineAHI(65:88),'.');

%% Plot specific signals
if 0
figure(2)
ax1(1)=subplot(4,1,1); plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,[9 14]));
ax1(2)=subplot(4,1,2); plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,5));
ax1(3)=subplot(4,1,3); plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,[2 ]));
ax1(4)=subplot(4,1,4); plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,[10:12 15]));

linkaxes(ax1,'x')
end
%% Analyse Loop gain and arousal threshold data 

minNevents=1;
maxwakethres=30;
maxFREM=0;

clear LG1 LG2 LGn temp temp2 FVAL VRA1 VRA2 ArThres MeanEx delay

% events, arousals, longestwake, position, FremMAX
x=[minNevents      0         maxwakethres   NaN  0];

%position

%
usefirstX=0;
uselastX=0;
usefirstpercent=0;
%usefirstpercent=0;
% 10.0000    0.4756   -0.9000    1.0000    0.6000    0.0873

usebest=0; Nbest=10;
usemediannotmeanLG1=1;
rem_subjects=1; %ie

removethese=find(Exclude==1);%find((BaselineAHI<20)|(BMI>37.5)) %43 cm / 17 inches is commonly used.

%clear Tn LG1 delay

    for i=1:M
        temp=LGplusinfo{i};
        try
            tempLG1=temp(:,7); 
        catch me
            tempLG1=NaN; 
        end
        FisNaN(i)=sum(isnan(tempLG1))/length(tempLG1);
    end
    
    for i=1:M
        try
        temp=LGplusinfo{i};
        if isempty(LGplusinfo{i})
            LG1(i,1)=NaN;
            LGn(i,1)=NaN;
            LG2(i,1)=NaN;
            delay(i,1)=NaN;
            MeanEx(i,1)=NaN;
            LGn(i,1)=NaN;
            Tn(i,1)=NaN;
            VRA2(i,1)=NaN;
            VRA1(i,1)=NaN;
            ArThres(i,1)=NaN;
            LG3min(i,1)=NaN;
            LG6min(i,1)=NaN;
            LG90s(i,1)=NaN;
        continue
        end
        
        LGdataallzeros=(0==sum(LGplusinfo{i}'))';
        
        tempLG1=temp(:,7); tempLGn=temp(:,5);  tempLG2=temp(:,8); tempTn=temp(:,6);
        tempLG0=temp(:,2); temptau=temp(:,3);  tempdelay=temp(:,9); 
        tempVRA1=temp(:,10); tempVRA2=temp(:,11); tempArThres=temp(:,12);
        temp2=LG_QualityInfo{i};
        temppositiondata=LG_QualityInfo{i}(:,5);
        templongestwakedata=SleepData{i}(:,7);
        tempFREM=SleepData{i}(:,6);
            templongestwakedata(length(tempLG1)+1:end)=[];
            tempFREM(length(tempLG1)+1:end)=[];
            %SleepData(winNum+1,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
        
        alpha=(tempLG2./tempLG1).^2;
        beta=(1-alpha)./(4*alpha-1);
        tempLG3min=tempLG1.*((1+beta)./(1+beta*(1/3)^2)).^0.5; %error, needs fixing...
        tempLG6min=tempLG1.*((1+beta)./(1+beta*(1/6)^2)).^0.5; %error, needs fixing...
        tempLG90s=tempLG1.*((1+beta)./(1+beta*(2/3)^2)).^0.5; %error, needs fixing...
             
        %GetRsquared:
        temp3=fitQual{i};
        clear OneMinusRsq
        for ii=1:size(temp3,2)
            temp3pt1=temp3{ii};
            if length(temp3pt1)>1
                OneMinusRsq(ii,:)=1-temp3pt1(2); %Element #2 of FitQual
            else
                OneMinusRsq(ii,:)=NaN;
            end
        end
        % Rsq=temp3(:,1);
        FVAL=temp(:,13);
        MeanE=temp2(:,3); %is actual scored events. LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position EXITFLAG];
        SStot=FVAL./OneMinusRsq;
        
        N_events=temp2(:,2); tempN_arousals=temp2(:,1);
        
        if 1
            %criteria=N_events>=x(4)&FVAL<x(1)&MeanE>x(2)&tempN_arousals>=x(5)&tempLG1>x(3)&MeanE<x(7)&tempLG1<x(8);
            criteria=templongestwakedata<=x(3)&tempN_arousals>=x(2)&N_events>=x(1)&tempFREM<=x(5)&~LGdataallzeros; %(temppositiondata==0|temppositiondata==x(4))&
        else
            criteria=ones(length(tempLG1),1)&~LGdataallzeros;
        end
        
        %Nbest=10; %usebest=0; %greatly improves arthres estimate. r>0.8!! used:x=[0.8^2    0.30    1    0.90    1.0000   0.22^2];
        if usebest %added criteria - take best 10 measures based on FVAL
            if 0
                FVAL2=FVAL; FVAL2(criteria==0)=NaN;
                tempsorted=sort(FVAL2); thres=tempsorted(Nbest+1); if isnan(thres),thres=Inf; end
                criteria=criteria&(FVAL2<thres);
            else
                OneMinusRsq2=OneMinusRsq; %OneMinusRsq2(criteria==0)=NaN;
                tempsorted=sort(OneMinusRsq2); thres=tempsorted(Nbest+1); if isnan(thres),thres=Inf; end
                if 0
                    criteria=criteria&(OneMinusRsq2<thres);
                else
                    criteria=criteria|(OneMinusRsq2<thres);
                end
            end
        end
        
        if usefirstX>0 %added criteria - take first 10 measures across the night
            tempfind=find(cumsum(criteria)>usefirstX,1);
            if length(tempfind)>0
                criteria(tempfind:end)=0;
            end
        end
        
        if usefirstX==0
        if uselastX>0 %added criteria - take first 10 measures across the night
            tempfind=find(cumsum(flipud(criteria))>=uselastX,1);
            if length(tempfind)>0
                criteria(1:end-tempfind)=0;
            end
        end
        end
        
        if usefirstpercent>0
            tempfind=find((100*cumsum(criteria)/sum(criteria))>usefirstpercent,1);
            if length(tempfind)>0
                criteria(tempfind:end)=0;
            end
        end
        
        if usemediannotmeanLG1
            LG1(i,1)=prctile(tempLG1(~isnan(tempLG1)&criteria),50);
        else
            LG1(i,1)=mean(tempLG1(~isnan(tempLG1)&criteria));
        end
        
        LG1_SD(i,1)=std(tempLG1(~isnan(tempLG1)&criteria));
        LG1_N(i,1)=sum(~isnan(tempLG1)&criteria); % no_of_LGmeasures(i,1)=sum(~isnan(tempLG1)&criteria);
        LG1_SE(i,1)=LG1_SD(i,1)/LG1_N(i,1)^0.5;
        LG1_N_nocriteria(i,1)=sum(~isnan(tempLG1)); % no_of_LGmeasures(i,1)=sum(~isnan(tempLG1)&criteria);
        
%         meanSE=mean(LG1_SE);
%         F_incl=LG1_N./LG1_N_nocriteria;
%         meanFincl=mean(F_incl);
%         minLGN=min(LG1_N);
        
        %minmeanmedianNmeasures=[min(no_of_LGmeasures) mean(no_of_LGmeasures) median(no_of_LGmeasures)]
        % LGn_=tempLGn(criteria);
        % VRA1_=tempVRA1(criteria);
        % VRA2_=tempVRA2(criteria);
        MeanEx(i,1)=prctile(MeanE(~isnan(MeanE)&criteria),50);
        % LG1(i,1)=prctile(tempLG1(criteria),50)
        LG2(i,1)=prctile(tempLG2(~isnan(tempLG2)&criteria),50);
        LGn(i,1)=prctile(tempLGn(~isnan(tempLGn)&criteria),50);
        Tn(i,1)=prctile(tempTn(~isnan(tempTn)&criteria),50);
        VRA1(i,1)=prctile(tempVRA1(~isnan(tempVRA1)&criteria&tempN_arousals>0),50);
        VRA2(i,1)=prctile(tempVRA2(~isnan(tempVRA1)&criteria&tempN_arousals>0),50);
        if tempN_arousals<1
                VRA1(i,1)=NaN;
                VRA2(i,1)=NaN;
        end
        ArThres(i,1)=prctile(tempArThres(~isnan(tempArThres)&criteria),50);
        % Rsq_median(i,1)=prctile(Rsq(criteria),50)
        LG0direct(i,1)=prctile(tempLG0(~isnan(tempLG0)&criteria),50);
        tau(i,1)=prctile(temptau(~isnan(temptau)&criteria),50);
        delay(i,1)=prctile(tempdelay(~isnan(tempdelay)&criteria),50);
        LG3min(i,1)=prctile(tempLG3min(~isnan(tempLG3min)&criteria),50);
        LG6min(i,1)=prctile(tempLG6min(~isnan(tempLG6min)&criteria),50);
        LG90s(i,1)=prctile(tempLG90s(~isnan(tempLG90s)&criteria),50);
        
        %limitVRA1=1.59;
        %VRA1_boundaryN(1,i,1)=sum(tempVRA1(~isnan(tempVRA1)&criteria)>limitVRA1);
        
        %within subject correlation
        if 0
        figure(101);
        set(gcf,'color',[1 1 1]);
        subplot(ceil((M-1)/10),11,i);
        plot(N_events(~isnan(tempLGn)&criteria),tempLGn(~isnan(tempLGn)&criteria),'.',0:20,0.7:0.7,':k');
        xlim([0 20]);
        ylim([0 2]);
        end
%         figure(102);subplot(ceil((M-1)/10),11,i);plot(N_events(~isnan(tempTn)&criteria),tempTn(~isnan(tempTn)&criteria),'.',0:20,40:40,':k');
%         xlim([0 20]);
%         ylim([0 100]);
        
        if length(tempLG1(~isnan(tempLG1)&criteria))>1
            [temptemp1,temptemp2]=corrcoef(tempLG1(~isnan(tempLG1)&criteria)',N_events(~isnan(tempLG1)&criteria)');
            R_LGvsNevents(1,i)=temptemp1(1,2);
            P_LGvsNevents(1,i)=temptemp2(1,2);
        else
            R_LGvsNevents(1,i)=NaN;
            P_LGvsNevents(1,i)=NaN;
        end
            mean_R_LGvsNevents(1)= mean(R_LGvsNevents(1,:));
            sem_R_LGvsNevents(1)= std(R_LGvsNevents(1,:))/length(R_LGvsNevents(1,:))^0.5;
        %try plotting all on top of each other
%         if 1==2
%             figure(102); hold('on');
%             plot(N_events(~isnan(tempLG1)&criteria),tempLG1(~isnan(tempLG1)&criteria),'.');
%             hold('off');
%         end
        
        catch me
                LG1(i,1)=NaN;
        LG2(i,1)=NaN;
        MeanEx(i,1)=NaN;
        LGn(i,1)=NaN;
        Tn(i,1)=NaN;
        VRA2(i,1)=NaN;
        VRA1(i,1)=NaN;
        ArThres(i,1)=NaN;
        LG3min(i,1)=NaN;
        LG6min(i,1)=NaN;
        LG90s(i,1)=NaN;
        delay(i,1)=NaN;
        end
    end %end for i... (subject loop)
    
    if rem_subjects&&length(removethese)>0
LG1(removethese)=NaN;
        LG2(removethese)=NaN;
        MeanEx(removethese)=NaN;
        LGn(removethese)=NaN;
        Tn(removethese)=NaN;
        VRA2(removethese)=NaN;
        VRA1(removethese)=NaN;
        ArThres(removethese)=NaN;
        LG3min(removethese)=NaN;
        LG6min(removethese)=NaN;
        LG90s(removethese)=NaN;
        delay(removethese)=NaN;
    end

    VRA1=100*VRA1;
    
%% Analyze Collapsibility and responsiveness data
 
 minNevents=0;
 maxwakethres=360;
 maxFREM=0;

ploton=0;
    minNevents = 0;
    unnormalize = 0;
    minwakethres = 330;
    NN = M;
    Vpassive=NaN*zeros(NN,1);
    Vactive=NaN*zeros(NN,1);
    medianV=NaN*zeros(NN,1);
    arthres=NaN*zeros(NN,1);
    for n=[1:NN]
        if ploton
            figure(1)
            tempX=5;
            subplot(tempX,ceil(NN/tempX),n);

            hold('on');
        end
    N_events = LG_QualityInfo{n}(:,2);
    Pos = LG_QualityInfo{n}(:,5);
    FREM=SleepData{n}(:,6);
    minwake = SleepData{n}(:,7);
    JJ=length(FREM);
    
    if length(N_events)<JJ, N_events(end:JJ)=NaN; end
    if length(Pos)<JJ, Pos(end:JJ)=NaN; end
        
        %criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=minwakethres;
        criteria = N_events>=minNevents & minwake<=minwakethres;
        if sum(criteria)==0
            continue
        end
        
    [t,x,y,a,win,nota,veup,hyp,e] = VEVdriveArray(DataOut,n,criteria);

    if unnormalize
        y = y.*veup;
        x = x.*veup;
        ymean = mean(y); %re-normalize with mean over whole night
        y = y/ymean;
        x = x/ymean;       
    end
    
    %arthres(n) = ArThresFromPSG(1,x,win);
    arthres(n) = ArThres(n);
    %arthres(n) = 1.6;   
    if 0
            x = x(nota==1);
            y = y(nota==1);
    else
            x = x(nota==1&hyp~=4);
            y = y(nota==1&hyp~=4);
    end
    
    if 1 %set minimum for arthres PER SE
            if arthres(n)<1.00
                arthres(n)=1.00;
            end
    end
    
    arthres2(n) = arthres(n);
    if arthres2(n)<1.05
        arthres2(n)=1.05;
    end
        
    medianV(n) = median(y);
    
    [Vpassive(n),Vactive(n)]=VEVdriveAnalysis(x,y,arthres2(n),ploton,1,10);
    if ploton
        title([num2str(BaselineAHI(n),2) '-' num2str(BaselineAHI(n)-DeltaAHI(n),2)]);
    end
    end  
    
    Vcomp = Vactive-Vpassive;
    
    Vcomp = 100*Vcomp;
    arthres = 100*arthres;
    Vactive = 100*Vactive;
    Vpassive = 100*Vpassive;

%% Best to predict Percent, absolute or treatment AHI?

Exclude(isnan(TreatmentAHI))=1;

figure(9)
XX=BaselineAHI(Exclude==0);
YY=TreatmentAHI(Exclude==0);
[r,p]=corrcoef(XX,YY); 
rp=[r(1,2) p(1,2)]
scatter(XX,YY,20,[0.2 0.2 0.5],'filled','markerfacealpha',0.5) 

XX=BaselineAHI(Exclude==0);
YY=DeltaAHI(Exclude==0);
[r,p]=corrcoef(XX,YY); 
rp=[r(1,2) p(1,2)]
scatter(XX,YY,20,[0.2 0.2 0.5],'filled','markerfacealpha',0.5) 

XX=BaselineAHI(Exclude==0);
YY=DeltaAHIp(Exclude==0);
[r,p]=corrcoef(XX,YY); 
rp=[r(1,2) p(1,2)]
scatter(XX,YY,20,[0.2 0.2 0.5],'filled','markerfacealpha',0.5) 

DeltaAHIf = (BaselineAHI-TreatmentAHI)./(BaselineAHI+TreatmentAHI);

XX=BaselineAHI(Exclude==0);
YY=DeltaAHIf(Exclude==0);
[r,p]=corrcoef(XX,YY); 
rp=[r(1,2) p(1,2)]
scatter(XX,YY,20,[0.2 0.2 0.5],'filled','markerfacealpha',0.5) 


    %% Predicting outcomes using single traits 

addpath([cd '\MatlabFunctionDatabase']);

xvalueslist = {'Vpassive','LGn','VRA1','arthres0p5','Vcomp'}; %
fx = 0.1+0*[1:length(xvalueslist)]; fxw = 0.05; %for bar locations
yvaluestr = {'DeltaAHIf'}; ythres=(1-0.3)/(1+0.3);
%yvaluestr = {'TreatmentAHI'}; ythres=10; 
%yvaluestr = {'DeltaAHI'}; ythres=20;

corrtype = 'Pearson';% '%Pearson' 'Spearman'

Nrows=5;

figure(82);
set(gcf,'color',[1 1 1]);

msize=15;
    markertype = '.';
    
AA=ceil(length(xvalueslist)/Nrows);
BB=ceil(length(xvalueslist)/AA);
clear thresopt
Y=eval(yvaluestr{1});
Y(Exclude==1)=NaN;
for i=1:length(xvalueslist)
    xvalue = eval(xvalueslist{i});
    
    Y=eval(yvaluestr{1});
    Y(isnan(xvalue))=[];
    xvalue(isnan(xvalue))=[];
    criteriaR = Y>ythres; 
    
    plotrange=1:length(Y);
    plotdata1=zeros(length(Y),1); plotdata1(plotrange)=1;
    criteriaR2 = Y>ythres;
    criteriaMid = Y>Inf;
    criteriaMid2 = Y>Inf;
    criteriaN = Y<=ythres%(~criteriaR2)&(~criteriaMid)&(~criteriaMid2)&~isnan(Y);
    
    subplot(AA,BB,i); 
    
    X=xvalue;
    [R,P]=corr(X(~isnan(Y)&~isnan(X)),Y(~isnan(Y)&~isnan(X)),'type',corrtype);

        funcsensspecstr = 'y+(1-x)';

[PPVNPVSensSpec,thresoptX,YvariableByPredictedOutcome,PredT,thresopt,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,AUCdata]=univariateleave1out(X,criteriaR,funcsensspecstr,Y);

X_ = [nanmean(X(PredT==1)) nanmean(X(PredT==0))];
if X_(1)<X_(2)
    pos=0;
else
    pos=1;
end
X_ = [nanmean(X(PredT==1-pos)) nanmean(X(PredT==pos))];
Y_ = [nanmean(Y(PredT==1-pos)) nanmean(Y(PredT==pos))];
Errors = [nanstd(Y(PredT==1-pos)) nanstd(Y(PredT==pos))]./[nansum(PredT==1-pos) nansum(PredT==pos)].^0.5;
[~,phighvslow]=ttest2(Y(PredT==1),Y(PredT==0)); 


thresopt_X(i)=thresopt;
plot(thresopt*[1 1],[min(Y),max(Y)],'k:');
    hold('on');
    temp=1;
    
    
    plot(xvalue,Y,markertype,'markersize',msize,'color',[1 1 1]) %dummy plot
    set(gca,'fontname','arial narrow','fontsize',15);
    %title([num2str(R,2),',',num2str(P,3),',',num2str(phighvslow,2)],'fontsize',8);
    
    if sum((X==1|X==0))==length(X) %binary
        xlim([-0.5 1.5]);
    end

    xlimtemp = get(gca,'xlim');
    xrangedelta = max(X)-min(X);
    
    stepsize = 10^round(log10(xrangedelta))/5;
    
    xlimmax = max(X)+0.1*xrangedelta;
        xlimmax = ceil(xlimmax/stepsize)*stepsize;
    xlimmin = min(X)-0.1*xrangedelta;
        xlimmin = floor(xlimmin/stepsize)*stepsize;
        
        if xlimtemp(1)<xlimmin
            xlimtemp(1)=xlimmin;
        end
        
        if xlimtemp(2)>xlimmax
            xlimtemp(2)=xlimmax;
        end
        
        xlim(xlimtemp);
        
    if sum((X==1|X==0))==length(X) %binary
        xlim([-0.5 1.5]);
        xlimtemp = get(gca,'xlim');
    end
    
Temp = xlimtemp(1)+fx(i)*(xlimtemp(2)-xlimtemp(1)); 
X_ = Temp + [-1 1]*fxw*(xlimtemp(2)-xlimtemp(1)); 
      
hold('on');
% BarColor = [0 0 0; 0 0 0];
% xlimtemp = get(gca,'xlim');
% barwidth = diff(xlimtemp)*0.08;
% hatchwidth = barwidth*0.25;
% plotbarwitherrors(Y_,X_,Errors,BarColor,0,barwidth,hatchwidth);
% 
    plot(xvalue(criteriaMid2&plotdata1==temp),Y(criteriaMid2&plotdata1==temp),markertype,'markersize',msize,'color',[1 0.5 0]) %partial nonresponders
    hold('on');
    plot(xvalue(criteriaN&plotdata1==temp),Y(criteriaN&plotdata1==temp),markertype,'markersize',msize,'color',[1 0 0]);
    plot(xvalue(criteriaMid&plotdata1==temp),Y(criteriaMid&plotdata1==temp),markertype,'markersize',msize,'color',[0.9 0.75 0]); %partial responders light green: 0.4 0.8 0
    plot(xvalue(criteriaR2&plotdata1==temp),Y(criteriaR2&plotdata1==temp),markertype,'markersize',msize,'color',[0 0.5 0]); %M, on top
    box('off');
    hold('off');
    set(gca,'fontname','arial narrow','fontsize',15);
%     
    %title(['R:',num2str(R,2),'[',num2str(P,3),'%,',num2str(phighvslow,2),',AUC:',num2str(AUCdata(1),2),'±',num2str(AUCdata(2),2),'[',num2str(AUCdata(3),2),'],PPV/NPV:',num2str(PPVNPVSensSpec(1)*100,2),'/',num2str(PPVNPVSensSpec(2)*100,2)],'fontsize',8);
    title(['R:',num2str(R,2),'[Pcorrel=',num2str(P,3),'],Pgroup=',num2str(phighvslow,2)],'fontsize',8);

    if i==1
        ylabel(yvaluestr{1});
    end
    xlabel(xvalueslist{i});
    ylim([min(Y) max(Y)]);
hold('off');
end


%% Predicting outcomes using multiple traits 

%% Setup predictors "Amatrix"
arthres_=arthres; arthres(arthres<100)=100;
arthres0p5 = 100+(100*((arthres-100)/100).^0.5);

Vpassive_=Vpassive; Vpassive_(Vpassive>100)=100;
Vpassive0p5 = 100-(100*((100-Vpassive_)/100).^0.5);
Vpassive0p33 = 100-(100*((100-Vpassive_)/100).^(1/3));

%'Vcomp','Tn','delay',
clear PredT remlist2 table1 tempXX
ythres=70
testthreshold=ythres; 
BaselineAHI = BaselineAHI(1:length(LG1)); %,'BaselineAHI','BMI','Sex','Age','Neck'
BaselineAHI = BaselineAHI(:)
xvalueslist={'BaselineAHI','BMI','Age','Sex','Neck'};
%xvalueslist={'Vpassive','LGn','LG1','arthres','Vcomp','Tn','delay','BaselineAHI','VRA1','Vactive','LG2'};%}; %'VpassiveT','VcompT','arthresT',
xvalueslist={'Vpassive','LGn','VRA1','Vactive','arthres0p5','Sex','Neck'};%}; %'VpassiveT','VcompT','arthresT',
%xvalueslist={'Vpassive','LGn','VRA1','Vactive','arthres0p5'};%}; %'VpassiveT','VcompT','arthresT',
xvalueslist={'Vpassive','LGn','VRA1','Vactive','arthres0p5','Neck','BMI','BaselineAHI','Sex','Age'};%}; %'VpassiveT','VcompT','arthresT',
    %very good, especially without subj81, and using quadratic model
    %xvalueslist={'Vpassive','LGn','VRA1','Vactive','arthres0p5'};%}; %'VpassiveT','VcompT','arthresT',
xvalueslist={'Vpassive','LGn','VRA1','Vactive','arthres0p5'};%}; %'VpassiveT','VcompT','arthresT',

xvalueslist={'Vpassive0p5','LGn','Vcomp','arthres0p5','VRA1'};%}; %'VpassiveT','VcompT','arthresT',

if 0
xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2'};%}; %'VpassiveT','VcompT','arthresT',

Vcomp0p5=Vcomp;
Vcomp0p5(Vcomp<0) = -(-Vcomp(Vcomp<0)).^0.5;
Vcomp0p5(Vcomp>=0) = (Vcomp(Vcomp>=0)).^0.5;
%xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2','Vpassive.^0.5','LGn.^0.5','VRA1.^0.5','Vcomp0p5','arthres0p5.^0.5','LG1','delay','Vactive','BaselineAHI','Neck','BMI','Sex','Age'};%}; %'VpassiveT','VcompT','arthresT',
xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2','Vpassive.^3','LG1','LG1.^2','delay','Age'};%}; %'VpassiveT','VcompT','arthresT',
end



%xvalueslist={'Vpassive','LGn','VRA1','Vactive','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vactive.^2','arthres0p5.^2'};%}; %'VpassiveT','VcompT','arthresT',

yvaluestr = {'DeltaAHIp'}; 
Y=eval(yvaluestr{1});
%setup Amatrix
Amatrix=[];
for i=1:length(xvalueslist)
    Amatrix = [Amatrix eval(xvalueslist{i})];
end%,arthresT,VcompT,LGnT,VactiveT,LG1T]; %VactiveT,LG1T VcompT./arthresT

% 
% if 0 %import and pool data with Marques study
%     Amatrix_=Amatrix; Y_ = Y;
%     load('AmatrixPlus.mat','Amatrix','Y');
%     M1 = length(Y);
%     Amatrix = [Amatrix;Amatrix_];
%     DeltaAHIp = [Y;Y_];
%     %rewrite variables
%     Yvariable = DeltaAHIp;
%     BaselineAHI = Amatrix(:,8);
%     TreatmentAHI = Amatrix(:,8) - Amatrix(:,8).*Yvariable/100;
% end
% 
% if 1 %import and pool data with Marques study
%     Amatrix_=Amatrix; Y_ = Y;
%     load('AmatrixPlus.mat','Amatrix','Y');
%     M1 = length(Y);
%     Amatrix = [Amatrix;Amatrix_];
%     DeltaAHIp = [Y;Y_];   
% end

if 0 %makes things difficult later
Irem = sum(isnan(Amatrix),2)>1;
Amatrix(Irem==1,:)=[];
DeltaAHIp(Irem==1)=[];
end
% 
% if 1
%     arthres0p5 = 100+(100*((Amatrix(:,4)-100)/100).^0.5);
%     Amatrix(:,4) = arthres0p5;
%     
%      Yvariable = DeltaAHIp;
%     BaselineAHI = Amatrix(:,8);
%     TreatmentAHI = Amatrix(:,8) - Amatrix(:,8).*Yvariable/100;
%     DeltaAHI = BaselineAHI - TreatmentAHI;
%     Exclude(end:length(TreatmentAHI))=0;
% end

    
%% Logistic regression method (unused)
if 0
   testthreshold=70
funcsensspecstr = 'y+(1-x)';
%criteriaR = 1*(Yvariable>testthreshold); Yvariable=DeltaAHIp;
criteriaR = 1*(DeltaAHIp>50&TreatmentAHI<10); Yvariable=TreatmentAHI;
for Nfeatures_=1:2
maxNfeatures=Nfeatures_;
%tempignorevars = [2:8];
tempignorevars = [];
forcedvars = [];
breakifdeltacostlessthan=-Inf;
rangevalidate=1:M; %only works for rows 1 to N
excludefortraining=[]; %6 7 9
[PPVNPVSensSpec,table1,thresoptX,YvariableByPredictedOutcome,PredT,B,Bstats,Ilist,thresopt,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,~,AUCdata]=logregleave1out(Amatrix,criteriaR,tempignorevars,funcsensspecstr,Yvariable,maxNfeatures,breakifdeltacostlessthan,1,rangevalidate,excludefortraining,forcedvars);
logitvalue=log(thresopt./(1-thresopt));
AHIdifferentialsA = [YvariableByPredictedOutcome(3)-YvariableByPredictedOutcome(1)]
AHIdifferentialsB = [YvariableByPredictedOutcomeAll(3)-YvariableByPredictedOutcomeAll(1)]
normalized_cutoff = (B(1)-logitvalue)/B(2);
performance = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),Yvariable(rangevalidate));
% scores2 = B(2)*Xvariable + B(3)*Yvariable(rangevalidate);
% labels = 1*(criteriaR|criteriaMid);
% [X,Y,T,AUC,OPTROCPT2] = perfcurve(labels,scores2,1); %need to find the threshold value that gives the OPTROCPT!
% if AUC<0.5, AUC=1-AUC; end
% [ci,se,p] = AUCci(AUC,N1,N2);

%clear tempXX
tempXX(Nfeatures_,:) = [PPVNPVSensSpec YvariableByPredictedOutcome YvariableByPredictedOutcomeAll]
end
end
%% logreg explore linear and quadratic terms
if 0
    if 1
    %Amatrix([81],:)=NaN; 
    Amatrix([81],:)=NaN; %wass no data for 81 until the noise threshold was adjusted to allow data; should be excluded(?)
    
    end

    
    criteriaR = 1*(DeltaAHIp>ythres);
        criteriaR(isnan(DeltaAHIp))=NaN;
        rangevalidate=1:M;

%Lower loop gain and higher arthres -> responder in basic logreg analysis
%weak trends for higher Vpassive and higher Vcomp
Ilisttest = [1 2 3 4 5];
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)

%Lower loop gain and higher arthres -> responder in basic logreg analysis
%trends for higher Vpassive and higher Vcomp
Ilisttest = [1 4 6 7 10];
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)

%Lower loop gain and higher arthres -> responder in basic logreg analysis
%xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2'};%}; %'VpassiveT','VcompT','arthresT',
Ilisttest = [1 6 11];
xvalueslist(Ilisttest)
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)
XX=0:100;
Ytemp = Btemp(1)+Btemp(2)*XX+Btemp(3)*XX.^2 + Btemp(4)*XX.^3; 
Ytemp = 1./(1+exp(-(Ytemp)));
YY=criteriaR+0.05*randn(length(criteriaR),1);
scatter(Vpassive,YY,20,[0.4 0.4 1],'filled','markerfacealpha',0.7)
hold('on');
plot(XX,Ytemp);
hold('off');

Ilisttest = [1 6 11 5 2 ];
xvalueslist(Ilisttest)
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)

    %linear regr
    PredT_=NaN*BaselineAHI;
    for j=1:M+1
    Ilisttest = 1:length(xvalueslist);
        %xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2','Vpassive.^3','LGn.^3','VRA1.^3','Vcomp.^3','arthres0p5.^3','LG1','delay','Vpassive.^4','Vactive','BaselineAHI','Neck','BMI','Sex','Age'};%}; %'VpassiveT','VcompT','arthresT',
    pvals=Inf;
    Train=1:M;
    if j<M+1
        Train(j)=[];
    end
    for i=1:length(Ilisttest)
    xvalueslist(Ilisttest);
    [Btemp,~,temp] = glmfit(Amatrix(Train,Ilisttest),DeltaAHIf(Train)); %,'weights',weights
    pvals = temp.p(2:end);
    [~,ii]=max(pvals);
    %i=i+5;.
    if max(pvals)>0.1&length(Ilisttest)>2
        xvalueslist(Ilisttest(ii));
        Ilisttest(ii)=[];
    else
        break
    end
    end
    xvalueslist(Ilisttest)';
    PrPredict = Amatrix(:,Ilisttest)*Btemp(2:end)+Btemp(1);
    
    %figure(10)
    %scatter(PrPredict(Train),DeltaAHIf(Train),20,[0.4 0.4 1],'filled','markerfacealpha',0.7); hold('on')
    %scatter(PrPredict(Test),DeltaAHIf(Test),20,[1 0.4 0.4],'filled','markerfacealpha',0.7); hold('off')
    [x,y,t,AUC,~] = perfcurve(criteriaR(Train)*1,PrPredict(Train),1); %need to find the threshold value that gives the OPTROCPT!
    
             [~,I]=max(y+(1-x));
             thresopt=mean(t(I:(I+1)));
             
             if j<M+1
             PredT_(j) = (PrPredict(j)>thresopt)*1;
                if isnan(PrPredict(j))
                    PredT_(j) = NaN;
                end
             else
             PredT = (PrPredict>thresopt)*1;  
             PredT(isnan(PrPredict)) = NaN;
             xvalueslist(Ilisttest)'
             end
    end
    
    performanceP = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),DeltaAHIp(rangevalidate));         
    performanceT = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentAHI(rangevalidate));
    performanceB = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHI(rangevalidate));
    performanceD = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHI(rangevalidate)-TreatmentAHI(rangevalidate));
    
        

Ilisttest = [1 6 11 5 2];
xvalueslist(Ilisttest)
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)

Ilisttest = [20 21 23];
xvalueslist(Ilisttest)
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)

Ilisttest = [1 6 11 5 2];
xvalueslist(Ilisttest)
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
pvals = temp.p(2:end)

PrPredict = Amatrix(:,Ilisttest)*Btemp(2:end)+Btemp(1);
PrPredict = 1./(1+exp(-(PrPredict)));
figure(10)
YY=criteriaR+0.05*randn(length(criteriaR),1);
scatter(PrPredict,YY,20,[0.4 0.4 1],'filled','markerfacealpha',0.7)

funcsensspec = eval(funcsensspecstr); %y+(1-x)%y+(1-x); % min distance 1./sqrt((1-y).^2+(x).^2) % harmonic mean y.*(1-x)./(y+(1-x))
 
             [x,y,t,AUC,~] = perfcurve(criteriaR*1,PrPredict,1); %need to find the threshold value that gives the OPTROCPT!
             [~,I]=max(y+(1-x));
             thresopt=mean(t(I:(I+1)));
    PredT = PrPredict>thresopt;
    performanceP = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),DeltaAHIp(rangevalidate));         
    performanceT = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),TreatmentAHI(rangevalidate));
    performanceB = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),BaselineAHI(rangevalidate));
    performanceD = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),BaselineAHI(rangevalidate)-TreatmentAHI(rangevalidate));
    
             
%             if AUC<0.5
%                 AUC=1-AUC;
%                 y=1-y;
%                 x=1-x;
%             end
%             tempAUC(kk)=AUC; %AUC
%             tempmaxsensspec(kk)=max(eval(funcsensspecstr));
%     
    
    
end
    %%



%% Machine learning: Support vector machine 
%Takes a few min to complete backwards elimination
if 1
    %Amatrix([81],:)=NaN; 
    Amatrix([81],:)=NaN; %wass no data for 81 until the noise threshold was adjusted to allow data; should be excluded(?)
    
end

close all
%Yvariable=TreatmentAHI; testthreshold=10; criteriaR = 1*(Yvariable<testthreshold); 
%Yvariable=DeltaAHIp; testthreshold=70; criteriaR = 1*(Yvariable>testthreshold);
remove_ = (Exclude==1|isnan(TreatmentAHI)|isnan(DeltaAHIp)|isnan(DeltaAHI)|sum(isnan(Amatrix),2)>=1);

%Yvariable=TreatmentAHI; testthreshold=10; criteriaR = 1*(DeltaAHIp>50&TreatmentAHI<10);
testthreshold=70;%nanmedian(DeltaAHIp(remove_==0));
%testthreshold=nanmedian(DeltaAHIp);
%testthreshold=50;
%ythres=2/3*100;
%criteriaR = 1*(DeltaAHIp>50&TreatmentAHI<10);
criteriaR = 1*(DeltaAHIp>testthreshold); Yvariable=DeltaAHIp; %TreatmentAHI<10
%criteriaR = 1*(DeltaAHIp>70|TreatmentAHI<10);

%criteriaR = 1*(DeltaAHIp>100*2/3); Yvariable=DeltaAHIp; %TreatmentAHI<10

%criteriaR = 1*(DeltaAHIp>nanmedian(DeltaAHIp)); Yvariable=DeltaAHIp; %TreatmentAHI<10
%Yvariable=TreatmentAHI; criteriaR = 1*(Yvariable<10); 

%Yvariable=(TreatmentAHI<10)+(DeltaAHIp>75);
%criteriaR = 1*((TreatmentAHI<10)|(DeltaAHIp>75)); %TreatmentAHI<10

% criteriaR = 1*(DeltaAHIp>50&TreatmentAHI<10); 
%criteriaR = 1*(DeltaAHIp>70); Yvariable=DeltaAHIp;

Yvariable(remove_==1)=NaN;
criteriaR(remove_==1)=NaN;
Amatrix(remove_==1,:)=NaN;

M=length(Yvariable);
M_ = sum(~isnan(Yvariable)); 
M__ = sum(~isnan(Yvariable)&sum(isnan(Amatrix),2)==0);

    minN=1;
    clear tempXX2
    tempXX2=zeros(size(Amatrix,2),34)*NaN;
    for Nfeatures_= min([size(Amatrix,2) 5]):-1:2 %was 3, 2 for speed
        coststr = 'phighvslow_'; %coststr = 'phighvslow_'
        %coststr = '1-SensplusSpec_';
        sigma = 2; %Nfeatures_=2, use 1.5, Nfeatures=3 use 2
        svpmethodstr = 'QP'; %'SMO' 'QP'
        kernelfunction = 'quadratic'; %'linear' 'rbf' 'quadratic' 'polynomial'
            %quadratic model is better, and in theory we should be able to
            %get the equation out...
        %coststr = '1-SensplusSpec_';
        breakifdeltacostlessthan=-Inf;
        maxNfeatures=Nfeatures_;
        forcedvars = []; %forcedvars = [1 5 4]; [8 10]
        tempignorevars = [];
        excludefortraining=[]; rangevalidate=[1:M]; 
%         if 1 %only Cistulli data
%             excludefortraining=[1:M1]; rangevalidate=[(M1+1):M]; 
%         elseif 0 %Cistulli + Marques data
%             
%         elseif 0 %train all, validate Cistulli %%%%%%%%%%%%%%%
%             excludefortraining=[]; rangevalidate=[(M1+1):M]; 
%         elseif 0 %train all, validate Marques %%%%%%%%%%%%%%%
%             excludefortraining=[]; rangevalidate=[1:M1];       
%         elseif 0 %train all, validate Marques %%%%%%%%%%%%%%%
%             excludefortraining=[1:M1]; rangevalidate=[1:M1];           
%         else %only Marques data
%             excludefortraining=[(M1+1):M]; rangevalidate=[1:M1];     
%         end
        [PPVNPVSensSpec,table1,YvariableByPredictedOutcome,PredT,IlistA,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,PredTAll]=svmleave1out(Amatrix,criteriaR,tempignorevars,Yvariable,maxNfeatures,breakifdeltacostlessthan,rangevalidate,excludefortraining,forcedvars,coststr,svpmethodstr,kernelfunction,sigma);
        tempXX2(Nfeatures_,:) = [PPVNPVSensSpec PPVNPVSensSpecAll YvariableByPredictedOutcome YvariableByPredictedOutcomeAll]
        sensplusspec=tempXX2(:,1)+tempXX2(:,2);
        phighvlow=tempXX2(:,23);
        phighvlow_=tempXX2(:,26);
        if 0
            [~,maxi]=max(sensplusspec);
        else
            [~,maxi]=min(phighvlow_);
        end
        if maxi<Nfeatures_&&Nfeatures_>=length(forcedvars)&&Nfeatures_>minN
            break %model now deteriorating not improving, thus stop adding features.
        end
        Ilist_{Nfeatures_} = IlistA;
    end
    %%
    if 1 %find best in list and rerun
    maxNfeatures=max([maxi length(forcedvars) minN]);
    [PPVNPVSensSpec,table1,YvariableByPredictedOutcome,PredT,Ilist,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,PredTAll]=svmleave1out(Amatrix,criteriaR,tempignorevars,Yvariable,maxNfeatures,breakifdeltacostlessthan,rangevalidate,excludefortraining,forcedvars,coststr,svpmethodstr,kernelfunction,sigma);
    else %use last run's data
        
        Ilist=IlistA;
    end
    performance = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),Yvariable(rangevalidate))
    performanceAll = PredictiveValue(1*criteriaR(rangevalidate),PredTAll(rangevalidate),Yvariable(rangevalidate));
    
    chosenvariableslist = xvalueslist(Ilist)
    
    %xvalueslist(Ilist_{6})
    
    performanceT = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),TreatmentAHI(rangevalidate));
    performanceB = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),BaselineAHI(rangevalidate));
    performanceD = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),BaselineAHI(rangevalidate)-TreatmentAHI(rangevalidate));
    
    % plot SVM
     %find linear predictive equation, so far ok only with autoscale off.
        %Ilist = [1 2 4];
        %kernelfunction = 'rbf'; %'linear' 'rbf'
        SVMModel = svmtrain(Amatrix(:,Ilist),criteriaR,'autoscale','on','method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
        
        %for 'linear' method only:
        W = ((SVMModel.Alpha'*SVMModel.SupportVectors).*SVMModel.ScaleData.scaleFactor)';
        bias = SVMModel.Bias + sum(W'.*SVMModel.ScaleData.shift);
        B = [bias;W];
        
        Z = (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,SVMModel.ScaleData.scaleFactor.*(Amatrix(:,Ilist)+SVMModel.ScaleData.shift),SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias; 
        figure(22)
        Zplot=-Z; Zplot(Zplot>2)=2; Zplot(Zplot<-2)=-2;
        scatter(Zplot,DeltaAHIf,20,[1 0.2 0.2],'filled','markerfacealpha',0.7);
        PredTsvm = svmclassify(SVMModel,Amatrix(rangevalidate,Ilist));
        %performanceSVM = PredictiveValue(criteriaR(rangevalidate),PredTsvm(rangevalidate),Yvariable(rangevalidate));
        
        %% Find quadratic fit to data
        
        %create data
        clear X
        L=10;
        for i=1:5
        minX=min(Amatrix(:,Ilist(i)));
        maxX=max(Amatrix(:,Ilist(i)));
        dX=(maxX-minX)/(L-1);
        X{i}=(minX:dX:maxX)';
        end
        
        rowdata=NaN*zeros(L^5,5);
        count=1;
        for i=1:L
            for j=1:L
                for k=1:L
                    for l=1:L
                        for m=1:L
                            rowdata(count,:)=[X{1}(i) X{2}(j) X{3}(k) X{4}(l) X{5}(m)];
                            count=count+1;
                        end
                    end
                end
            end
        end
        
        PredTsvmA = svmclassify(SVMModel,rowdata);
        
        rowdata2 = [rowdata rowdata.^2 rowdata(:,1).*rowdata(:,2) rowdata(:,1).*rowdata(:,3) rowdata(:,1).*rowdata(:,4) rowdata(:,1).*rowdata(:,5) rowdata(:,2).*rowdata(:,3) rowdata(:,2).*rowdata(:,4) rowdata(:,2).*rowdata(:,5) rowdata(:,3).*rowdata(:,4) rowdata(:,3).*rowdata(:,5) rowdata(:,4).*rowdata(:,5)];
        
        [Btemp,~,temp] = glmfit(rowdata2,PredTsvmA,'binomial'); %,'weights',weights
        [Btemp temp.se temp.p]
        
        temppred0=rowdata2*Btemp(2:end)+Btemp(1);
        temppred = 1*(temppred0> 0) ;
        
        [x,y,t,AUC,~] = perfcurve(PredTsvmA*1,temppred,1); %need to find the threshold value that gives the OPTROCPT!
             [~,I]=max(y+(1-x));
             thresopt=mean(t(I:(I+1)));
    
        AmatrixX = [Amatrix(:,Ilist) Amatrix(:,Ilist).^2];
        
        PredTest = (AmatrixX*Btemp(2:end)+Btemp(1)) > 0.5 ;
        
        PredTsvm
        
        performance_est = PredictiveValue(1*criteriaR(rangevalidate),PredTest(rangevalidate),Yvariable(rangevalidate))
        
    
        %%
        
        
        
        [class,err,POSTERIOR,logp,coeff] = classify(Amatrix,rowdata,PredTsvmA,'Quadratic')
        [class,err,POSTERIOR,logp,coeff] = classify(Amatrix,Amatrix,criteriaR,'Quadratic')
        
        %Amatrix0 = (Amatrix-nanmean(Amatrix))./nanstd(Amatrix);
        Ilist=[1 2 5];
        [class,err,POSTERIOR,logp,coeff] = classify(Amatrix(:,Ilist),Amatrix(:,Ilist),criteriaR,'Quadratic');
        [X,Y] = meshgrid(linspace(0.2,0.8),linspace(100,250));
        X = X(:); Y = Y(:);
        [C,~,~,~,~] = classify([X Y],Amatrix(:,Ilist),criteriaR,'Quadratic');
        figure(34);
        scatter(X,Y,5,[C C 1-C]);
        hold('on');
        scatter(Amatrix(:,Ilist(1)),Amatrix(:,Ilist(2)),5,[criteriaR 1-criteriaR 1-criteriaR]);

        K=coeff(1,2).const;
        L=coeff(1,2).linear;
        Q=coeff(1,2).quadratic;
        
        myf = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);
        myf2 = @(A) K + [A]*L + sum(([A]*Q) .* [A], 2);
        
        myf2([80 0.35 122])
        
        h2 = ezplot(myf,[0.2 0.8 100 250]);
        set(h2,'Color','m','LineWidth',2)
        axis([4.5 8 2 4])
        
        performanceAll;
        performance_est2 = PredictiveValue(1*criteriaR(rangevalidate),class(rangevalidate),Yvariable(rangevalidate))
        
         %% separate training/validation
%         Ilist=[1 2 3 5];
%         excludefortraining=1:40;
%         rangevalidate=41:108;
%         SVMModel = svmtrain(Amatrix(excludefortraining,Ilist),criteriaR(excludefortraining),'autoscale','on','method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
%         PredT2 = svmclassify(SVMModel,Amatrix(rangevalidate,Ilist));
%         
%         performance = PredictiveValue(1*criteriaR(rangevalidate),PredT2,Yvariable(rangevalidate))
%         
%% decision tree version

%         Mdl = fitctree(Amatrix,criteriaR)
%         view(Mdl,'mode','graph')
%         PredTtree = predict(Mdl,Amatrix);
%         performance_tree = PredictiveValue(1*criteriaR(rangevalidate),PredTtree(rangevalidate),Yvariable(rangevalidate))
    Ilist=[1 2 4 5]
    tree_numsplits=10;
    tree_minLeafSize=6;
    tree_methodslist = {'interaction-curvature' 'AllSplits' 'curvature'}
    tree_method=tree_methodslist{1};
        PredTtree_=NaN*criteriaR;
        for i=1:M
            train=1:M;
            train(i)=[];
            Mdl = fitctree(Amatrix(train,Ilist),criteriaR(train),'MaxNumSplits',tree_numsplits,'MinLeafSize',tree_minLeafSize,'PredictorSelection',tree_method);
            PredTtree_(i,1) = predict(Mdl,Amatrix(i,Ilist));
        end
        performance_tree = PredictiveValue(1*criteriaR(rangevalidate),PredTtree_(rangevalidate),Yvariable(rangevalidate))
        
        performance_treeT = PredictiveValue(1*criteriaR(rangevalidate),PredTtree_(rangevalidate),TreatmentAHI(rangevalidate));
        Mdl = fitctree(Amatrix(:,Ilist),criteriaR(:),'MaxNumSplits',tree_numsplits,'MinLeafSize',tree_minLeafSize,'PredictorSelection',tree_method);
        PredTtree = predict(Mdl,Amatrix(:,Ilist));
        performance_treeall = PredictiveValue(1*criteriaR(rangevalidate),PredTtree(rangevalidate),Yvariable(rangevalidate));
        view(Mdl,'mode','graph')

        %%
        if 0 %not sure what this is in retrospect
        PredTAll_left = PredTAll;
        temp = Amatrix(:,10)<70&PredTAll==0;
        PredTAll_left(temp==1)=NaN;
        col=4
        
        [~,p]=ttest2(Amatrix(PredTAll_left==0,col),Amatrix(PredTAll_left==1,col))
        [mean(Amatrix(PredTAll_left==0,col)) mean(Amatrix(PredTAll_left==1,col)) std(Amatrix(PredTAll_left==0,col))/sum(PredTAll_left==0).^0.5 std(Amatrix(PredTAll_left==1,col))/sum(PredTAll_left==0).^0.5 p]
        end
        
        %% Force features for 3d plot
        if 1 %find best in list and rerun
        Ilist = [1 2 3]
        SVMModel = svmtrain(Amatrix(:,Ilist),criteriaR,'autoscale','on','method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
        
        %for 'linear' method only:
        W = ((SVMModel.Alpha'*SVMModel.SupportVectors).*SVMModel.ScaleData.scaleFactor)';
        bias = SVMModel.Bias + sum(W'.*SVMModel.ScaleData.shift);
        B = [bias;W];
        
        Z = (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,Amatrix(:,Ilist),SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias; 
        %figure()
        %plot(Z,Yvariable,'.');
        PredTsvm = svmclassify(SVMModel,Amatrix(rangevalidate,Ilist));
        %performanceSVM = PredictiveValue(criteriaR(rangevalidate),PredTsvm(rangevalidate),Yvariable(rangevalidate));
        end
        
%% Plot classifier

plotrange=[1:size(Amatrix,1)]; %plot Cistulli + Marques 

xcolor = 1;
        if 1%sum(Yvariable==DeltaAHIp)==M
            if 0
            criteriaRplot = Yvariable>testthreshold;
            criteriaMid = Yvariable>50&~criteriaRplot;
            criteriaMid2 = Yvariable>33.33&~criteriaMid&~criteriaRplot;
            criteriaN = (~criteriaRplot)&(~criteriaMid)&(~criteriaMid2)&~isnan(Yvariable);
            elseif 0
            criteriaRplot = Yvariable>testthreshold&TreatmentAHI<10;
            criteriaMid = Yvariable>testthreshold&~criteriaRplot;
            criteriaMid2 = Yvariable>50&TreatmentAHI<10&~criteriaMid&~criteriaRplot;
            criteriaN = (~criteriaRplot)&(~criteriaMid)&(~criteriaMid2)&~isnan(Yvariable);
            else
            criteriaRplot = Yvariable>testthreshold;
            criteriaMid = Yvariable>Inf;
            criteriaMid2 = Yvariable>Inf;
            criteriaN = Yvariable<=testthreshold;     
            end
        else
            criteriaRplot = Yvariable<10;
            criteriaMid = Yvariable<15&~criteriaRplot;
            criteriaMid2 = Yvariable<20&~criteriaMid&~criteriaRplot;
            criteriaN = (~criteriaRplot)&(~criteriaMid)&(~criteriaMid2)&~isnan(Yvariable);            
        end
        %NaN outside plotrange
            temp = criteriaRplot; criteriaRplot=criteriaRplot*NaN; criteriaRplot(plotrange) = temp(plotrange);
            temp = criteriaMid; criteriaMid=criteriaMid*NaN; criteriaMid(plotrange) = temp(plotrange);
            temp = criteriaMid2; criteriaMid2=criteriaMid2*NaN; criteriaMid2(plotrange) = temp(plotrange);
            temp = criteriaN; criteriaN=criteriaN*NaN; criteriaN(plotrange) = temp(plotrange);
            
        if length(Ilist)==1
            figure(); set(gcf,'color',[1 1 1]);
            
            xrange = [min(Amatrix(plotrange,Ilist(1))) max(Amatrix(plotrange,Ilist(1)))];
            tempA = (xrange(2)-xrange(1))/5;
            tempB = (xrange(2)+xrange(1))/2;
            temp = Amatrix(criteriaRplot==1,Ilist(1)); temp = tempB + tempA*randn(length(temp),1);
            plot(Amatrix(criteriaRplot==1,Ilist(1)),temp,'.','markersize',20,'color',[0 0.5 0]);
            set(gca,'box','off','fontname','arial narrow','fontsize',14);
            hold('on');
            temp = Amatrix(criteriaMid==1,Ilist(1)); temp = tempB + tempA*randn(length(temp),1);
                plot(Amatrix(criteriaMid==1,Ilist(1)),temp,'.','markersize',20,'color',[xcolor 0.8 0]);
            temp = Amatrix(criteriaMid2==1,Ilist(1)); temp = tempB + tempA*randn(length(temp),1);
            plot(Amatrix(criteriaMid2==1,Ilist(1)),temp,'.','markersize',20,'color',[1 0.5 0]);
            temp = Amatrix(criteriaN==1,Ilist(1)); temp = tempB + tempA*randn(length(temp),1);
            plot(Amatrix(criteriaN==1,Ilist(1)),temp,'.','markersize',20,'color',[1 0 0]);

                yrange = xrange;
                [X,Y] = meshgrid(linspace(xrange(1),xrange(2)),linspace(yrange(1),yrange(2)));
                Xorig = X; Yorig = Y;
                if ~isempty(SVMModel.ScaleData) %scale the mesh
                    X = SVMModel.ScaleData.scaleFactor(1) * (X + SVMModel.ScaleData.shift(1));
                    Y = X;
                end
                Z = (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,X(:),SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias;
                contour(Xorig,Yorig,reshape(Z,size(X)),[0 0],'k');
            xlabel(xvalueslist{Ilist(1)});
            ylim_ = get(gca,'ylim');
            set(gca,'ylim',(ylim_-mean(ylim_))*2)
            set(gca,'ycolor',[1 1 1]);
            
        elseif length(Ilist)==2
            figure(); set(gcf,'color',[1 1 1]);
            plot(Amatrix(criteriaRplot==1,Ilist(1)),Amatrix(criteriaRplot==1,Ilist(2)),'.','markersize',20,'color',[0 0.5 0]);
            set(gca,'box','off','fontname','arial narrow','fontsize',14);
            hold('on');
            plot(Amatrix(criteriaMid==1,Ilist(1)),Amatrix(criteriaMid==1,Ilist(2)),'.','markersize',20,'color',[xcolor 0.8 0]);
            plot(Amatrix(criteriaMid2==1,Ilist(1)),Amatrix(criteriaMid2==1,Ilist(2)),'.','markersize',20,'color',[1 0.5 0]);
            plot(Amatrix(criteriaN==1,Ilist(1)),Amatrix(criteriaN==1,Ilist(2)),'.','markersize',20,'color',[1 0 0]);
            if 1
                xrange = [min(Amatrix(plotrange,Ilist(1))) max(Amatrix(plotrange,Ilist(1)))];
                yrange = [min(Amatrix(plotrange,Ilist(2))) max(Amatrix(plotrange,Ilist(2)))];
                lims = [xrange yrange];
                [X,Y] = meshgrid(linspace(xrange(1),xrange(2)),linspace(yrange(1),yrange(2)));
                Xorig = X; Yorig = Y;
                if ~isempty(SVMModel.ScaleData) %scale the mesh
                    X = SVMModel.ScaleData.scaleFactor(1) * (X + SVMModel.ScaleData.shift(1));
                    Y = SVMModel.ScaleData.scaleFactor(2) * (Y + SVMModel.ScaleData.shift(2));
                end
                Z = (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,[X(:),Y(:)],SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias; 
                contour(Xorig,Yorig,reshape(Z,size(X)),[0 0],'k');
                contour(Xorig,Yorig,reshape(Z,size(X)),[-0.1 -0.1],'k');
                contour(Xorig,Yorig,reshape(Z,size(X)),[0.1 0.1],'k');
            else
                b0 = bias;
                b1 = W(1);
                b2 = W(2);
                xplot=[min(Amatrix(:,Ilist(1))) max(Amatrix(:,Ilist(1)))];
                yplot = (b0+b1*xplot)/(-b2);
                yplotB = (b0-1+b1*xplot)/(-b2);
                yplotC = (b0+1+b1*xplot)/(-b2);
                plot(xplot,yplot,'-','color',[0 0 0]);
                plot(xplot,yplotB,'-','color',[0.5 0.5 0.5]);
                plot(xplot,yplotC,'-','color',[0.5 0.5 0.5]);
                axis([xplot min(Amatrix(:,Ilist(2))) max(Amatrix(:,Ilist(2)))]);
            end
            xlabel(xvalueslist{Ilist(1)});
            ylabel(xvalueslist{Ilist(2)});
        elseif length(Ilist)==3
            % 3d figure in spheres
            compressxyz=[1 1 1]; %[1 0.5 0.5]
            figure();
            set(gcf,'color',[1 1 1]);
            %htemp=scatter3(Amatrix(:,Ilist(1)),Amatrix(:,Ilist(2)),Amatrix(:,Ilist(3)));
            htemp=scatter3(Amatrix(:,Ilist(1)),Amatrix(:,Ilist(2)),Amatrix(:,Ilist(3)),'marker','none');
            xrange = [min(Amatrix(plotrange,Ilist(1))) max(Amatrix(plotrange,Ilist(1)))];
            yrange = [min(Amatrix(plotrange,Ilist(2))) max(Amatrix(plotrange,Ilist(2)))];
            zrange = [min(Amatrix(plotrange,Ilist(3))) max(Amatrix(plotrange,Ilist(3)))];
            cubevol = diff(xrange)*diff(yrange)*diff(zrange);
            %ballsize=0.02*cubevol;
            ballsize=1*max([diff(xrange) diff(yrange) diff(zrange)]);
            %ballsize=0.02*cubevol^(1/3);
            camerapos=get(gca,'CameraPosition');
            
%             criteriaR = Yvariable>66.67;
%             criteriaMid = Yvariable>50&~criteriaR;
%             criteriaMid2 = Yvariable>33.33&~criteriaMid&~criteriaR;
%             criteriaN = (~criteriaR)&(~criteriaMid)&(~criteriaMid2)&~isnan(Yvariable);
            C=zeros(size(Amatrix,1),3)*NaN;
            for i=1:size(Amatrix,1)
                if criteriaRplot(i)==1
                    C(i,:)=[0 0.5 0];
                elseif criteriaMid(i)==1
                    C(i,:)=[1 0.8 0];
                elseif criteriaMid2(i)==1
                    C(i,:)=[1 0.5 0];
                elseif criteriaN(i)==1
                    C(i,:)=[1 0 0];
                end
            end
            
            for i=1:length(Yvariable)
                distancefromeye(i)=sqrt((Amatrix(i,Ilist(1))-camerapos(1)).^2+(Amatrix(i,Ilist(2))-camerapos(2)).^2+(Amatrix(i,Ilist(3))-camerapos(3)).^2);
            end
            distancefromeye=distancefromeye/nanmean(distancefromeye);
            distancefromeye=distancefromeye.^1.5;
            distancefromeye=1./distancefromeye;
            hold('on');
            h=scatter3sph(Amatrix(plotrange,Ilist(1)),Amatrix(plotrange,Ilist(2)),Amatrix(plotrange,Ilist(3)),'size',ballsize*distancefromeye(plotrange),'color',C(plotrange,:),'transp',1,'compress',compressxyz);
            hold('on');
            xrange = [min(Amatrix(plotrange,Ilist(1))) max(Amatrix(plotrange,Ilist(1)))];
            yrange = [min(Amatrix(plotrange,Ilist(2))) max(Amatrix(plotrange,Ilist(2)))];
            zrange = [min(Amatrix(plotrange,Ilist(3))) max(Amatrix(plotrange,Ilist(3)))];
            
            Fextraxyz=0.02;
            xlim([min(xrange)-Fextraxyz*diff(xrange) max(xrange)+Fextraxyz*diff(xrange)]);
            ylim([min(yrange)-Fextraxyz*diff(yrange) max(yrange)+Fextraxyz*diff(yrange)]);
            zlim([min(zrange)-Fextraxyz*diff(zrange) max(zrange)+Fextraxyz*diff(zrange)]);
            %surface "wall"
            if 1
                svm_3d_matlab_vis(SVMModel,Amatrix(:,Ilist),criteriaR,1);
            else
                m= -B(2)/B(3); c=(-B(1)-B(4)*zrange(1))/B(3); %at z=zmin
                xmin = xrange(1); xmax=xrange(2);
                ymin = yrange(1); ymax=yrange(2);
                [xmin,xmax,ymin2,ymax2] = xybox(xmin,xmax,ymin,ymax,m,c);
                dZ=diff(zrange);
                dX=-dZ*B(4)/B(3)/(B(2)/B(3)+1); %from 0=B(2)/B(3)*dX + dY + B(4)/B(3)*dZ; perhaps B(2)/B(3)-1 for dY=-dX
                Xrange = [xmin xmax xmax+dX xmin+dX];
                zintent = [min(Amatrix(plotrange,Ilist(3))) min(Amatrix(plotrange,Ilist(3))) max(Amatrix(plotrange,Ilist(3))) max(Amatrix(plotrange,Ilist(3)))]
                Yrange=(0-B(1)-B(2)*Xrange-B(4)*zintent)/B(3);
                h1=fill3(Xrange,Yrange,zintent,[0.5 0.5 0.5])
                set(h1,'edgecolor','none','facealpha',0.33);
            end
            
            set(gca,'xgrid','on','ygrid','on','zgrid','on')
            
            %lighting for depth
            set(h,'facelighting','phong','ambientstrength',0.5);
            h2=light('position',[nanmean(Amatrix(plotrange,Ilist(1))) nanmean(Amatrix(plotrange,Ilist(2))) nanmean(Amatrix(plotrange,Ilist(3)))],'style','local');
            set(gca,'CameraPosition',camerapos,'projection','perspective','gridlinestyle','-','TickLength',[0 0],'fontsize',20);
            set(gca,'fontname','arial narrow');
            
            xlabel(xvalueslist{Ilist(1)});
            ylabel(xvalueslist{Ilist(2)});
            zlabel(xvalueslist{Ilist(3)});
        end
        
        % 3D case
        %plot3(X(:,1), X(:,2), (-bias-w1*X(:,1)-w2*X(:,2))/w3 )
    
%% Basic logistic regression
if 0
Ilisttest = [2 8]
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest(1)),criteriaR(:),'binomial'); %,'weights',weights
    temp.p(2)
[Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
         temp.p  
         
         figure(4)
         plot(criteriaR+0.03*(randn(length(criteriaR),1)),Amatrix(:,2),'.','markersize',8); xlim([-2 3])
end

%% Results: treatment AHI<10, linear SVM: LGn<0.43 (weak predictor)
% B=[-0.4888;1.1321];
% temppredict = B(1)+B(2)*LGn;
% LGn_cutoff=-B(1)/B(2);
% plot(criteriaR,LGn,'.','markersize',10); xlim([-2 3]);
% hold('on');
% plot([-0.5 1.5],[1 1]*LGn_cutoff,'k')
% [mean(TreatmentAHI(temppredict>0)) mean(TreatmentAHI(temppredict<=0))];
% [~,p]=ttest2((TreatmentAHI(temppredict>0)),(TreatmentAHI(temppredict<=0)));
% 
% LGn_cutoff2=0.43;
% [mean(TreatmentAHI(LGn>LGn_cutoff2)) mean(TreatmentAHI(LGn<=LGn_cutoff2))];
% [~,p]=ttest2((TreatmentAHI(LGn>LGn_cutoff2)),(TreatmentAHI(LGn<=LGn_cutoff2)));
% 
%% Plot bar graphs of key variables for actual and predicted responders.
variables = {'DeltaAHIp','TreatmentAHI'} %BaselineAHI
%PredTAll=PredT
%PredT=PredT_
Nplots=2;
if 1
plotmedinstead = [1 1 1]; %plots do not work with new matlab (2017 N)
else
    plotmedinstead = [0 0 0]; %plots do not work with new matlab (2017 N)
end
figure(23); set(gcf,'color',[1 1 1]);
fontsize=8;
%criteriaR((M1+1):M)=NaN;
baseline=0;
for i=1:length(variables)
    subplot(Nplots,length(variables),i);
    X_ = [0.25 0.75];
    Y_ = [nanmean(1*eval([variables{i} '(criteriaR==1)'])) nanmean(1*eval([variables{i} '(criteriaR==0)']))];
    Errors = [nanstd(1*eval([variables{i} '(criteriaR==1)']))/sum(~isnan(1*eval([variables{i} '(criteriaR==1)'])))^0.5 nanstd(1*eval([variables{i} '(criteriaR==0)']))/sum(~isnan(1*eval([variables{i} '(criteriaR==0)'])))^0.5];
    if ~plotmedinstead(i)
        plotbarwitherrors(Y_,X_,Errors,[0 0 0;0 0 0],0,0.4,0.1);
    else
        Q1 = [prctile(1*eval([variables{i} '(criteriaR==1)']),25) prctile(1*eval([variables{i} '(criteriaR==0)']),25)];
        Q3 = [prctile(1*eval([variables{i} '(criteriaR==1)']),75) prctile(1*eval([variables{i} '(criteriaR==0)']),75)];
        plotbarwitherrorsmedian(Y_,X_,Q3,Q1,[0 0 0;0 0 0],baseline,0.4,0.5);
    end
    xlim([0 1]);
    set(gca,'xtick',[],'XColor',[1 1 1],'fontsize',fontsize);
end
% for i=1:length(variables)
%     subplot(Nplots,length(variables),i+length(variables));
%     X_ = [0.25 0.75];
%     Y_ = [nanmean(1*eval([variables{i} '(PredT==1)'])) nanmean(1*eval([variables{i} '(PredT==0)']))];
%     Errors = [nanstd(1*eval([variables{i} '(PredT==1)']))/sum(~isnan(1*eval([variables{i} '(PredT==1)'])))^0.5 nanstd(1*eval([variables{i} '(PredT==0)']))/sum(~isnan(1*eval([variables{i} '(PredT==0)'])))^0.5];
%     if ~plotmedinstead(i)
%         plotbarwitherrors(Y_,X_,Errors,[0 0 0;0 0 0],0,0.4,0.1);
%     else
%         Y_ = [nanmedian(1*eval([variables{i} '(PredT==1)'])) nanmedian(1*eval([variables{i} '(PredT==0)']))];
%         Q1 = [prctile(1*eval([variables{i} '(PredT==1)']),25) prctile(1*eval([variables{i} '(PredT==0)']),25)];
%         Q3 = [prctile(1*eval([variables{i} '(PredT==1)']),75) prctile(1*eval([variables{i} '(PredT==0)']),75)];
%         plotbarwitherrorsmedian(Y_,X_,Q3,Q1,[0 0 0;0 0 0],baseline,0.4,0.5);
%     end
%     xlim([0 1]);
%     set(gca,'xtick',[],'XColor',[1 1 1],'fontsize',fontsize);
% end
for i=1:length(variables)
    subplot(Nplots,length(variables),i+1*length(variables));
    X_ = [0.25 0.75];
    Y_ = [nanmean(1*eval([variables{i} '(PredT_==1)'])) nanmean(1*eval([variables{i} '(PredT_==0)']))];
    Errors = [nanstd(1*eval([variables{i} '(PredT_==1)']))/sum(~isnan(1*eval([variables{i} '(PredT_==1)'])))^0.5 nanstd(1*eval([variables{i} '(PredT_==0)']))/sum(~isnan(1*eval([variables{i} '(PredT_==0)'])))^0.5];
    if ~plotmedinstead(i)
        plotbarwitherrors(Y_,X_,Errors,[0 0 0;0 0 0],0,0.4,0.1);
    else
        Y_ = [nanmedian(1*eval([variables{i} '(PredT_==1)'])) nanmedian(1*eval([variables{i} '(PredT_==0)']))];
        Q1 = [prctile(1*eval([variables{i} '(PredT_==1)']),25) prctile(1*eval([variables{i} '(PredT_==0)']),25)];
        Q3 = [prctile(1*eval([variables{i} '(PredT_==1)']),75) prctile(1*eval([variables{i} '(PredT_==0)']),75)];
        plotbarwitherrorsmedian(Y_,X_,Q3,Q1,[0 0 0;0 0 0],baseline,0.4,0.5);
    end
    xlim([0 1]);
    set(gca,'xtick',[],'XColor',[1 1 1],'fontsize',fontsize);

end

%% Does physiology subgroup add to traditional measure to predict outcomes?

% criteriaR = TreatmentAHI<10;
% %criteriaR = DeltaAHIp>70;
% 
% PredX = [Age Neck Sex BaselineAHI BMI];
% PredX(Exclude==1,:)=NaN;
% [B,dev,stats] = glmfit(PredX,criteriaR(:),'binomial'); %,'weights',weights
%     stats.p
% [B1,dev1,stats1] = glmfit([PredX PredTAll],criteriaR(:),'binomial'); %,'weights',weights
%     stats1.p
    
%% Linear regression to predict outcomes (Amatrix)
if 0
PredX = Amatrix;
PredX(Exclude==1,:)=NaN;
[B,dev,stats] = glmfit(PredX,DeltaAHIp(:),'normal'); %,'weights',weights
    stats.p
    [B,dev,stats2] = glmfit(PredX,criteriaR(:),'binomial'); %,'weights',weights
    stats2.p
end

%% Does physiology subgroup add to traditional measure to predict outcomes?


PredX = [Age Neck Sex BaselineAHI BMI];
PredX(Exclude==1,:)=NaN;
[B,dev,stats] = glmfit(PredX,DeltaAHIp(:),'normal'); %,'weights',weights
    stats.p
[B1,dev1,stats1] = glmfit([PredX PredTAll],DeltaAHIp(:),'normal'); %,'weights',weights
    stats1.p
%% Does physiology subgroup add to traditional measure to predict outcomes?

PredX = [Age Neck Sex BaselineAHI BMI];
PredX(Exclude==1,:)=NaN;
[B,dev,stats] = glmfit(PredX,TreatmentAHI(:),'normal'); %,'weights',weights
    stats.p
[B1,dev1,stats1] = glmfit([PredX PredT],TreatmentAHI(:),'normal'); %,'weights',weights
[B1,dev1,stats1] = glmfit([PredX PredT_],TreatmentAHI(:),'normal'); %,'weights',weights
    stats1.p
    
    
    %% Bivariate stats, clinical vs outcome
    clear Data
    PredX = [Age Neck Sex BaselineAHI BMI BaselineAHI.^0.33];
    for i=1:size(PredX,2)
    [B1,dev1,stats1] = glmfit(PredX(Iincl,i),DeltaAHIf(Iincl),'normal'); %,'weights',weights
    Data(i,:) = [B1(2) stats1.p(2)];
    end
    
    %% Clinical plus subgroup vs outcome
    clear Data Data_
    Scoring4p = zeros(M,1);
    Scoring4p(65:88)=1;
    
    PredX = [Age*0+1 Age Neck Sex BaselineAHI BMI BaselineAHI.^0.33 Scoring4p];
    for i=1:size(PredX,2)
    [B1,dev1,stats1] = glmfit([PredX(Iincl,i) PredT(Iincl)],DeltaAHIf(Iincl),'normal'); %,'weights',weights
    Data(i,:) = [B1(2) stats1.p(2)];
    end
    for i=1:size(PredX,2)
    [B1,dev1,stats1] = glmfit([PredX(Iincl,i) PredT_(Iincl)],DeltaAHIf(Iincl),'normal'); %,'weights',weights
    Data_(i,:) = [B1(2) stats1.p(2) B1(3) stats1.p(3)];
    end
    
    %% Multivariate stats
    clear Data
    PredX = [Age Neck Sex BMI BaselineAHI.^0.33];
    [B1,dev1,stats1] = glmfit(PredX(Iincl,:),DeltaAHIf(Iincl),'normal'); %,'weights',weights
    [B1 stats1.p]
    tabledata = [PredX(Iincl,:) DeltaAHIf(Iincl)];
    
    %% Bivariate stats, traits vs absolute outcome
    clear Data
    PredX = [Vpassive0p5 LGn Vcomp arthres0p5 VRA1 Age Neck Sex BaselineAHI BMI BaselineAHI.^0.33];
    for i=1:size(PredX,2)
    [B1,dev1,stats1] = glmfit([PredX(Iincl,i) BaselineAHI(Iincl).^0.33],DeltaAHI(Iincl),'normal'); %,'weights',weights
    Data(i,:) = [B1(2) stats1.p(2)];
    end
    

    
%% Plot actual data
PredX = PredT;
figure(5); clf(5);
subplot(1,3,1);
n=length(BaselineAHI(PredX==1&Exclude==0));
scatter(0.1*randn(1,n)+1,BaselineAHI(PredX==1&Exclude==0),20,[0 0 0],'filled','markerfacealpha',0.3); 
hold('on');
n=length(BaselineAHI(PredX==0&Exclude==0));
scatter(0.1*randn(1,n)+3,BaselineAHI(PredX==0&Exclude==0),20,[0 0 0],'filled','markerfacealpha',0.3); 
xlim([0 4]);    

subplot(1,3,2);
n=length(DeltaAHIp(PredX==1&Exclude==0));
scatter(0.1*randn(1,n)+1,DeltaAHIp(PredX==1&Exclude==0),20,[0 0 0],'filled','markerfacealpha',0.3); 
hold('on');
n=length(DeltaAHIp(PredX==0&Exclude==0));
scatter(0.1*randn(1,n)+3,DeltaAHIp(PredX==0&Exclude==0),20,[0 0 0],'filled','markerfacealpha',0.3); 
xlim([0 4]);

subplot(1,3,3);
n=length(TreatmentAHI(PredX==1&Exclude==0));
scatter(0.1*randn(1,n)+1,TreatmentAHI(PredX==1&Exclude==0),20,[0 0 0],'filled','markerfacealpha',0.3); 
hold('on');
n=length(TreatmentAHI(PredX==0&Exclude==0));
scatter(0.1*randn(1,n)+3,TreatmentAHI(PredX==0&Exclude==0),20,[0 0 0],'filled','markerfacealpha',0.3); 
xlim([0 4]);  
%%
    median(TreatmentAHI(PredX==0&Exclude==0))
    median(TreatmentAHI(PredX==1&Exclude==0))
    
%% Proportions
I = PredT_==1;
A=10;
[sum(DeltaAHIp(I)>50&TreatmentAHI(I)<A) sum(DeltaAHIp(I)>50&TreatmentAHI(I)>=A) sum(DeltaAHIp(I)<50)]
I = PredT_==0;
[sum(DeltaAHIp(I)>50&TreatmentAHI(I)<A) sum(DeltaAHIp(I)>50&TreatmentAHI(I)>=A) sum(DeltaAHIp(I)<50)]

    
    