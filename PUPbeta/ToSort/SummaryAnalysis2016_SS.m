%% Load data
clear
close all

%% set defaults for figures
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesTickDir','out');

%% Load PUPdata

savename = 'J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\Analyzed\FlowAndPes2.mat'

% ShowFigures = 1;
%savename = 'J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\Analyzed\FlowAndPes2.mat'
savename = 'C:\PSG_Data\AirflowDrive_Data\Analyzed\FL_DLM_withDataOut.mat'

% load selectively
load(savename,'SleepData','LG_QualityInfo','DataOut',...
    'BreathDataTable','BreathFLDataTable','LocalSignals');

load(savename,'SpO2data','SleepData','MiscEpochData','Ratings','AnalysisIndex','LGplusinfo','EventsInfo','LG_QualityInfo','DataOut','ArousalDat','fitQual','StoNData','AHIdata','CPAPData')

M=length(LGplusinfo);

%%
clear AHITable
for m=1:M
    AHITable(m,:) = AHIdata{m};
end

clear SpO2Table
for m=1:M
    SpO2Table(m,:) = SpO2data{m};
end

AHInrem = AHITable(:,88);
AHItotal = AHITable(:,80);

AHInrem1 = AHITable(:,112); %allpos
AHInrem2 = AHITable(:,120);
AHInrem3 = AHITable(:,128);
durnrem1 = AHITable(:,112-7); %allpos
durnrem2 = AHITable(:,120-7);
durnrem3 = AHITable(:,128-7);
AHInrem1(durnrem1<5)=NaN; AHInrem2(durnrem2<5)=NaN; AHInrem3(durnrem3<5)=NaN;
deltaAHIp_N1toN2 = (AHInrem1-AHInrem2)./AHInrem1*100; deltaAHIp_N1toN2(AHInrem1<15)=NaN;
deltaAHI_N1toN2 = (AHInrem1-AHInrem2); deltaAHI_N1toN2(AHInrem1<15)=NaN;
deltaAHIp_N2toN3 = (AHInrem2-AHInrem3)./AHInrem2*100; deltaAHIp_N2toN3(AHInrem2<15)=NaN;
deltaAHI_N2toN3 = (AHInrem2-AHInrem3); deltaAHI_N2toN3(AHInrem2<15)=NaN;

Fhypops = (AHITable(:,85)+AHITable(:,87))./AHInrem*100;
NadirNREMSpO2 = SpO2Table(:,10);

%% Import quickly in one open/close from matfiles

[~,~,raw]=xlsread('AnalyzeDataSpreadsheet',1,['G3:G' num2str(3+M-1)]);
for i=1:length(raw), DPW(i,1)=str2double(raw{i}(1:end-4)); end
matfiledirectory = 'J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\'
matvarlist = {'Veupnea','Vactive','Pcrit','Vpassive','ArThresPes','ArthresEdi','Varousal'};
matvarlistfield = {'mean','Vcrit_Feupnea','Pcritvalue','Feupnea','median','Feupnea','Feupnea'};

for m=1:M
    for i=1:length(matvarlist)
        try
            eval([matvarlist{i} '_(m,1)=NaN;']); %default value NaN
            eval(['clear ' matvarlist{i}]); %clear last struct loaded with same name
            load([matfiledirectory raw{m}],matvarlist{i}); %try to load struct
            eval([matvarlist{i} '_(m,1)=' matvarlist{i} '.' matvarlistfield{i} ';']);
        catch me
        end
    end
end

for m=1:M
    try
            VpassiveN_(m,1)=NaN; VpassiveSEM_(m,1)=NaN;
            clear Pcrit
            load([matfiledirectory raw{m}],'Pcrit'); %try to load struct
            VpassiveN_(m,1)=length(Pcrit.include_data);
            VpassiveSEM_(m,1)=Pcrit.VcritSEM/(Veupnea_(m))*100;
    catch me
    end
end
for m=1:M
    try
            VactiveN_(m,1)=NaN; VactiveSEM_(m,1)=NaN;
            clear Vactive
            load([matfiledirectory raw{m}],'Vactive'); %try to load struct
            VactiveN_(m,1)=length(Vactive.Vactive_list);
            VactiveSEM_(m,1)=Vactive.Vcrit_Feupnea_SEM*100;
    catch me
    end
end
for m=1:M
    try
            PVslope_(m,1)=NaN; 
            clear Pcrit
            load([matfiledirectory raw{m}],'Pcrit'); %try to load struct
            PVslope_(m,1)=Pcrit.PVSlope_VE/(Veupnea_(m))*100;
    catch me
    end
end

Vactive_ = 100*Vactive_;
Vpassive_= 100*Vpassive_;
ArthresEdi_=100*ArthresEdi_;
Varousal_=100*Varousal_;

clear filename_
for i=1:length(raw)
    filename_{i,1} = raw{i}(1:end-4);
end

%
[nanmean(VpassiveN_) nanstd(VpassiveN_); ...
nanmean(VactiveN_) nanstd(VactiveN_); ...
nanmean(VpassiveSEM_) nanstd(VpassiveSEM_); ...
nanmean(VactiveSEM_) nanstd(VactiveSEM_); ...
nanmean(PVslope_) nanstd(PVslope_); ...
nanmean(VpassiveSEM_)/nanmean(PVslope_) nanmean(VactiveSEM_)/nanmean(PVslope_); ...
]

%% Import xls data

if 0
    %O2PSGJune2016
    col1='A';
    col2='HJ'
    row1=2;
    row2=47;
    sheet=1;
    %read MAT filenames from xls
    filexls='C:\Users\szs88\Dropbox (Personal)\LG Surgery Database\MASTER SPREADSHEET.xlsx';
    
    loadnewvariables={'BaselineAHI','BMI','Age','Sex','TreatmentAHI'};
    colsnewvariables=[ 67            44    41    40   117 ];
    
    %rangerowsread = [max(rowsnewvariables)];
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
    
    DeltaAHIp = (BaselineAHI-TreatmentAHI)./BaselineAHI*100;
end
%% Exclude patients
Exclude = zeros(M,1);
%Exclude(AHInrem<5)=1
Exclude(AHItotal<5)=1;

%% Analyse PUP loop gain and arthres data

clear LG1 LG2 LGn temp temp2 FVAL VRA1 VRA2 ArThres MeanEx

% events, longestwake, position, FremMAX
minNevents = 0;
maxwakethres = 330; %change to 330 to obtain more data, but then relies on specialized LH scoring.
maxFREM = 0;

%position

%
usefirstX=0;
uselastX=0;
usefirstpercent=0;
%usefirstpercent=0;
% 10.0000    0.4756   -0.9000    1.0000    0.6000    0.0873

usebest=0; Nbest=10;
usemediannotmeanLG1=1;
rem_subjects=0; %ie

removethese=[Exclude];%find((BaselineAHI<20)|(BMI>37.5)) %43 cm / 17 inches is commonly used.

clear  LG1 LG2 MeanEx LGn Tn VRA2 VRA1 ArThres LG3min LG6min LG90s LG1_N


for i=1:M
    try
        tempLG1=LGplusinfo{i}(:,7);
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
            LG1_N(i,1)=NaN;
            
            continue
        end
        
         LGdataallzeros=(0==sum(LGplusinfo{i}'))';
         
        tempLG1=LGplusinfo{i}(:,7); tempLGn=LGplusinfo{i}(:,5);  tempLG2=LGplusinfo{i}(:,8); tempTn=LGplusinfo{i}(:,6);
        tempLG0=LGplusinfo{i}(:,2); temptau=LGplusinfo{i}(:,3);  tempdelay=LGplusinfo{i}(:,9);
        tempVRA1=LGplusinfo{i}(:,10); tempVRA2=LGplusinfo{i}(:,11); tempArThres=LGplusinfo{i}(:,12);
        
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
        %temp3=fitQual{i};
        clear OneMinusRsq
        for ii=1:size(fitQual{i},2)
            temp3pt1=fitQual{i}{ii};
            if length(temp3pt1)>1
                OneMinusRsq(ii,:)=1-temp3pt1(2); %Element #2 of FitQual
            else
                OneMinusRsq(ii,:)=NaN;
            end
        end
        % Rsq=temp3(:,1);
        FVAL=LGplusinfo{i}(:,13);
        
        MeanE=LG_QualityInfo{i}(:,3); %is actual scored events. LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position EXITFLAG];
        SStot=FVAL./OneMinusRsq;
        
        N_events=LG_QualityInfo{i}(:,2); tempN_arousals=LG_QualityInfo{i}(:,1);
        
        if 1
            %criteria=N_events>=x(4)&FVAL<x(1)&MeanE>x(2)&tempN_arousals>=x(5)&tempLG1>x(3)&MeanE<x(7)&tempLG1<x(8);
            criteria=templongestwakedata<=maxwakethres&N_events>=minNevents&tempFREM<=maxFREM&~LGdataallzeros; %(temppositiondata==0|temppositiondata==x(4))&
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
        
        %meanSE=mean(LG1_SE);
        F_incl(i,1)=LG1_N(i,1)./LG1_N_nocriteria(i,1);
        %meanFincl=mean(F_incl);
        %minLGN=min(LG1_N);
        
        %minmeanmedianNmeasures=[min(no_of_LGmeasures) mean(no_of_LGmeasures) median(no_of_LGmeasures)]
        % LGn_=tempLGn(criteria);
        % VRA1_=tempVRA1(criteria);
        % VRA2_=tempVRA2(criteria);
        MeanEx(i,1)=prctile(MeanE(~isnan(MeanE)&criteria),50);
        % LG1(i,1)=prctile(tempLG1(criteria),50)
        LG2(i,1)=prctile(tempLG2(~isnan(tempLG2)&criteria),50);
        LGn(i,1)=prctile(tempLGn(~isnan(tempLGn)&criteria),50);
        Tn(i,1)=prctile(tempTn(~isnan(tempTn)&criteria),50);
        VRA1(i,1)=prctile(tempVRA1(~isnan(tempVRA1)&criteria),50);
        VRA2(i,1)=prctile(tempVRA2(~isnan(tempVRA1)&criteria),50);
        ArThres(i,1)=prctile(tempArThres(~isnan(tempArThres)&criteria),50);
        ArThresCOV(i,1)=nanstd(tempArThres(~isnan(tempArThres)&criteria))/ArThres(i,1);
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
        if 1
            figure(101);
            %set(gcf,'color',[1 1 1]);
            subplot(ceil((M-1)/10),11,i);
            plot(N_events(~isnan(tempLGn)&criteria),tempLGn(~isnan(tempLGn)&criteria),'.',0:20,0.7:0.7,':k');
            xlim([0 20]);
            ylim([0 2]);
            box('off');
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
        ['error:' num2str(i)]  
        LG1(i,1)=NaN;
        LG2(i,1)=NaN;
        MeanEx(i,1)=NaN;
        LGn(i,1)=NaN;
        Tn(i,1)=NaN;
        VRA2(i,1)=NaN;
        VRA1(i,1)=NaN;
        ArThres(i,1)=NaN;
        ArThresCOV(i,1)=NaN;
        LG3min(i,1)=NaN;
        LG6min(i,1)=NaN;
        LG90s(i,1)=NaN;
        LG1_N(i,1)=NaN;
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
    ArThresCOV(removethese)=NaN;
    LG3min(removethese)=NaN;
    LG6min(removethese)=NaN;
    LG90s(removethese)=NaN;
    LG1_N(removethese)=NaN;
end
ArThres(ArThres==0)=NaN;
ArThresPSG1=ArThres*100;

[nanmean(ArThresCOV) nanstd(ArThresCOV)]*100;
%reported data with minwake threshold=30 s to best reflect future clinical use.

%% UA phenotype using model drive
unnormalize = 0;

fontsize_=8;

Vpassive=NaN*zeros(M,1);
Vactive=NaN*zeros(M,1);
VpassiveCI=NaN*zeros(M,2);
VactiveCI=NaN*zeros(M,2);
medianV=NaN*zeros(M,1);
arthresPSG=NaN*zeros(M,1);
NwindowsPSG=NaN*zeros(M,1);
arthresPSG2=NaN*zeros(M,1);
minNwindows = 1;
maxFREM=0;
minNevents=0;
maxwakethres=360;
Nciles=10;


artN1PSG=NaN*zeros(M,1);
artN2PSG=NaN*zeros(M,1);
artN3PSG=NaN*zeros(M,1);
artN1PSG_Nwin=NaN*zeros(M,1);
artN2PSG_Nwin=NaN*zeros(M,1);
artN3PSG_Nwin=NaN*zeros(M,1);

for n=1:M
    if Exclude(n)==1||LG1_N(n)<=minNwindows
        continue
    end
    figure(2);
    
    %title(filename_{n});
    try
        N_events = LG_QualityInfo{n}(:,2);
        %Pos = LG_QualityInfo{n}(:,5);
        FREM = SleepData{n}(:,6);
        Fwake = SleepData{n}(:,1);
        minwake = SleepData{n}(:,7);
        %criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=maxwakethres;
        criteria = N_events>=minNevents & minwake<maxwakethres & FREM<=maxFREM;
        NwindowsPSG(n) = sum(criteria);
        if sum(criteria)==0
            continue
        end
        subplot(ceil(sqrt(M)),ceil(M/ceil(sqrt(M))),n);
        
        plot(1,1,'marker','none'); hold('on');
        title(filename_{n});
        hold('on');
        [t x y a win nota veup hyp e] = VEVdriveArray(DataOut,n,criteria);
        
        x_=1;
        if unnormalize
            y = y.*veup;
            x = x.*veup;
            ymean = mean(y); %re-normalize with mean over whole night
            y = y/ymean;
            x = x/ymean;
        end
        y=y*100;
        
        
        %arthres(n) = ArThresFromPSG(1,x,win);
        %increasingdriveonly, deleteifbeloweupnea(0),
        %setaseupneaifbeloweupnea(1), usemediannotROC(1), falsepositiverate(0.2)
        if 0
            settings=[1 0 1 0 1 0 1 0.2]; 
            %ignorefirstXbreathsaftersleeponset, swapbreathsforhigherdrive, Nbreathsattributabletoarousal, ...
            %increasingdriveonly, deleteifbeloweupnea, setaseupneaifbeloweupnea
            %[arthresPSG(n)] = ArThresNew(a',x',x_,win,settings);
            [arthresPSG(n)] = ArThresNew(a',[x(2:end);NaN]',x_,win,settings);
        elseif 0
            %arthres(n) = ArThresFromPSG(1,x,win);
            arthresPSG(n) = ArThresPSG1(n)/100; % taken from window summary data, median
        elseif 0
            ArOnset = [NaN;diff(a)];    
            arthresPSG(n)= median(x(ArOnset==1));
        else 1 %window by window average as per LGfromFlow (near exact)
            temparthres=NaN*ones(max(win),1);
            temphyp=NaN*ones(max(win),1);
            ArOnset = [NaN;diff(a)];  
            %elast=[NaN;e(1:end-1)];  
            for i=1:max(win)     
                temparthres(i)=nanmean(x(ArOnset==1&win==i)); %&elast==0
                temphyp(i)=mode(hyp(win==i));
            end
            if 0 %sleep state correction factor
                temparthres=temparthres.*(1+(1-temphyp)*0.18);
            end
            arthresPSG(n)=nanmedian(temparthres);

                artN1PSG(n)=nanmedian(temparthres(temphyp==2)); artN1PSG_Nwin(n)=sum(~isnan(temparthres(temphyp==2)));
                artN2PSG(n)=nanmedian(temparthres(temphyp==1)); artN2PSG_Nwin(n)=sum(~isnan(temparthres(temphyp==1)));
                artN3PSG(n)=nanmedian(temparthres(temphyp==0)); artN3PSG_Nwin(n)=sum(~isnan(temparthres(temphyp==0)));
                

%                 [artN1(n),artN1_N(n)] = ArThresNew(a',tempx',x_,win,settings2);
%                 tempx=x; tempx(hyp~=1)=NaN;
%                 [artN2(n),artN2_N(n)] = ArThresNew(a',tempx',x_,win,settings2);
%                 tempx=x; tempx(hyp~=0)=NaN;
%                 [artN3(n),artN3_N(n)] = ArThresNew(a',tempx',x_,win,settings2);
                
        end
        if 0
            figure(34)
            plot(ArThresPSG1/100,arthresPSG,'.');
        end
        if 1 %set minimum for arthres PER SE
            if arthresPSG(n)<1.00
                arthresPSG(n)=1.00;
            end
        end
        arthresPSG2(n) = arthresPSG(n);
        if arthresPSG2(n)<1.05
            arthresPSG2(n)=1.05;
        end

        
        %arthres(n) = 1.6;
        if 0
            x = x(nota==1);
            y = y(nota==1);
        else
            x = x(nota==1&hyp~=4);
            y = y(nota==1&hyp~=4);
        end
        
        
        medianV(n) = median(y);
        
        ploton=1;
        
        plot([1 1],[0 150],'--','color',[0.7 0.7 0.7]);
        hold('on');
        plot([0 3*1],[100 100],'--','color',[0.7 0.7 0.7]);
        [Vpassive(n),Vactive(n),~,~,VpassiveCI(n,:),VactiveCI(n,:)]=VEVdriveAnalysis(x,y,arthresPSG2(n),ploton,x_,Nciles);
        set(gca,'fontsize',fontsize_)
        ylim([0 max([prctile(x,95) 1])]);
        ylim([0 max([prctile(y,95) 1])]);
        ylabel('Ventilation, %eupnea'); xlabel('Vdrive model, %eupnea');
        hold('off');
        
        figure(23)
        subplot(ceil(sqrt(M)),ceil(M/ceil(sqrt(M))),n);
        plot(1,1,'marker','none'); hold('on');
        plot([1 1],[0 150],'--','color',[0.7 0.7 0.7]);
        hold('on');
        plot([0 3*1],[100 100],'--','color',[0.7 0.7 0.7]);
        plot(x/x_*100,y,'k.','markersize',2); box('off');
        
        %title(['ahi:',num2str(round(BaselineAHI(n)),'%u') '->' num2str(round(TreatmentAHI(n)),'%u') ',LGn:' num2str(LGn(n),2)],'fontsize',8);
    catch me
        Vpassive(n)=NaN;
        Vactive(n)=NaN;
        medianV(n)=NaN;
        arthres(n)=NaN;
    end
end

Vcomp = Vactive-Vpassive;
VpassivePSG=Vpassive;
VactivePSG=Vactive;
VcompPSG=Vcomp;
VactivePSGCI=VactiveCI; VpassivePSGCI=VpassiveCI;
clear Vactive Vpassive Vcomp VpassiveCI VactiveCI

% 95%CIs
temp=abs(VpassivePSGCI-[VpassivePSG VpassivePSG]);
temp=mean(temp')';
[nanmean(temp),nanstd(temp)]/1.96

temp2=abs(VactivePSGCI-[VactivePSG VactivePSG]);
temp2=mean(temp2')';
[nanmean(temp2),nanstd(temp2)]/1.96

%% UA phenotype info using Edi/Pes
addpath('E:\Work\MatlabFunctionDatabase');
supinecodes = [0 2 -5]; %Profusion:[1],Spike:[0 2]
maxwakethres = 360;
%minwakearthres = 60;
minNevents = 0;
maxFREM=0;
drivecol=19;   %%[16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak 24:FlowEdi, 23:FlowPes]
dividedrivebyttot=0;
normalizeyusingveupwindow=1;
normalizeyusingknownVeupnea=0;
    Gequals1=0;
    
    Fvol=0; Fpcw=0;

%close all
try close(3); catch me, end
try close(4); catch me, end
try close(5); catch me, end

x_n=NaN*zeros(M,1);
y_n=NaN*zeros(M,1);

x_wake=NaN*zeros(M,1);
y_wake=NaN*zeros(M,1);

Vpassive=NaN*zeros(M,1);
Vactive=NaN*zeros(M,1);
VpassiveCI=NaN*zeros(M,2);
VactiveCI=NaN*zeros(M,2);
medianV=NaN*zeros(M,1);
arthres=NaN*zeros(M,1);
arthres2=NaN*zeros(M,1);
Nwindows=NaN*zeros(M,1);
arthres_N=NaN*zeros(M,1);
arthres2_N=NaN*zeros(M,1);
x_n=NaN*zeros(M,1);
minNwindows = 3;
gvar=NaN*zeros(M,1);
G=NaN*zeros(M,1);
G_N=NaN*zeros(M,1);
threshold=NaN*zeros(M,1);
rwake=NaN*zeros(M,1);
GGpeakwake=NaN*zeros(M,1);

artN1=NaN*zeros(M,1);
artN2=NaN*zeros(M,1);
artN3=NaN*zeros(M,1);
artN1_N=NaN*zeros(M,1);
artN2_N=NaN*zeros(M,1);
artN3_N=NaN*zeros(M,1);

%get sleep eupnea first
maxwakethres2=maxwakethres; 
for n=1:M
    if Exclude(n)==1||LG1_N(n)<=minNwindows
        continue
    end
    try
        N_events = LG_QualityInfo{n}(:,2);
        Pos = LG_QualityInfo{n}(:,5);
        Fwake = SleepData{n}(:,1);
        FREM = SleepData{n}(:,6);
        minwake = SleepData{n}(:,7);
        supine = 0*Pos;
        for j=1:length(supinecodes)
            supine = supine + 1*(Pos==supinecodes(j));
        end
        criteria = N_events>=minNevents & minwake<maxwakethres2 & FREM<=maxFREM ; %& (~isnan(Pos)&supine==1)
        Nwindows(n)=sum(criteria);
        if sum(criteria)==0
            continue
        end
        [t,x,y,a,win,nota,veup,hyp,pcw,t2,e] = VEPesArray(DataOut,n,criteria,drivecol);
        
        if normalizeyusingknownVeupnea
            y_ = Veupnea_(n)/60;
            y = y./y_*100;
        else
            if normalizeyusingveupwindow
                y = y./veup*100; %note Y is unnormalized VE from LGfromflow
                y_ = nanmedian(veup);
            else
                y_ = nanmedian(veup);
                y = y./y_*100;
            end
        end
        y_n(n)=y_;
    catch me
    end
end

%get wake data next (27,23,17) 
for n=1:M
    %ginput(1)
    if Exclude(n)==1||LG1_N(n)<=minNwindows
        continue
    end
    try
        N_events = LG_QualityInfo{n}(:,2);
        Pos = LG_QualityInfo{n}(:,5);
        Fwake = SleepData{n}(:,1);
        FREM = SleepData{n}(:,6);
        minwake = SleepData{n}(:,7);
        %criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=maxwakethres;
        %             supine = 0*Pos;
        %             for j=1:length(supinecodes)
        %                 supine = supine + 1*(Pos==supinecodes(j));
        %             end
        criteria = minwake>=0; %(~isnan(Pos)&supine==1)
        sum(criteria);
        if sum(criteria)==0
            continue
        end

        %[t,x,y,a,win,nota,veup,hyp,pcw,t2] = VEPesArray(DataOut,n,criteria,drivecol); %  original
        [t,x,y,a,win,nota,veup,hyp,pcw,t2,e,FL] = VEPesBreathFLArray(DataOut,pt,criteria,drivecol,BreathFLDataTable); % DLM edit
  
        if dividedrivebyttot
            x=x./t2;
        end
        [ad]=howfarawayfromsleep(a,win); %replace a with (a==1)&(hyp==4) if arousals not scored in wake...
                
        %find breaths during wakefulness and arousals
        if 1
            minNwakebreaths = 50;
            hh=hist(ad,[1:11]); hh(end)=NaN; th=find(hh>minNwakebreaths,1,'last');
            threshold(n) = min([4 th]);
        else
            threshold(n) = 2;
        end
        a1 = ad>threshold(n);
        G(n) = nanmedian(y(a1)./x(a1));
        
        if (Fvol>0||Fpcw>0) && drivecol==16
            x = x + Fvol*(y/G(n)) + Fpcw*pcw; 
            G(n) = nanmedian(y(a1)./x(a1)); %overwrite
        end
        
        if 1
            figure()
            plot(t(ad>=2),y(ad>=2)./x(ad>=2),'r.');
            hold('on');
            plot(t(a1),y(a1)./x(a1),'.');
            plot(t,nanmedian(y(ad>=2)./x(ad>=2))+0*t);
            temp = t(ad>=2);
            data = y(ad>=2)./x(ad>=2);
            [temp,ia,ic] = unique(temp);
            data = data(ia);
            data(data>prctile(data,99))=NaN;
            data(data<prctile(data,1))=NaN;
            temp(isnan(data))=[]; data(isnan(data))=[];
            maxt = 1800;
            Gmedian=NaN*temp;
            for j=1:length(temp)
                temp2 = abs(temp(j)-temp);
                weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
                %Gmedian(j)=nanmean(data(weights>0).*weights(weights>0));
                %Gmedian(j)=wtmedian(data(weights>0),weights(weights>0));
                Gmedian(j)=nansum(data(weights>0).*weights(weights>0));
                %Gmedian(j)=weightedMedian(data(weights>0),weights(weights>0));
            end
            Gtime{n}=sortrows([temp,Gmedian]);
            
            plot(Gtime{n}(:,1),Gtime{n}(:,2),'g-');
        
            
        end
        
        x_wake(n)=median(x(a1));
        y_wake(n)=median(y(a1));
        G_N(n) = length(y(a1)./x(a1));
        g_p25 = prctile(y(a1)./x(a1),25);
        g_p75 = prctile(y(a1)./x(a1),75);
        gvar(n) = 0.5*(g_p75-g_p25)/G(n);
        
        if 0
        figure(31)
        subplot(ceil(sqrt(M)),floor(sqrt(M)),n)
        plot(x(a1),y(a1),'.')
        corrtype='Spearman';
        [rwake(n)]=corr(x(~isnan(x)&~isnan(y)&a1),y(~isnan(x)&~isnan(y)&a1),'type',corrtype);
        end
    catch me
    end
end

[nanmean(x_wake) nanstd(x_wake); 60*nanmean(y_wake) 60*nanstd(y_wake)]

if Gequals1
    Gbackup=G;
    G=0*G+1;
end

try clf(3), catch me, end

%GET BREATH-BY-BREATH SLEEP DATA
for n=1:M
    figure(3);
    if Exclude(n)==1||LG1_N(n)<=minNwindows
        continue
    end
    try
        N_events = LG_QualityInfo{n}(:,2);
        Pos = LG_QualityInfo{n}(:,5);
        FREM = SleepData{n}(:,6);
        Fwake = SleepData{n}(:,1);
        minwake = SleepData{n}(:,7);
        supine = 0*Pos;
        for j=1:length(supinecodes)
            supine = supine + 1*(Pos==supinecodes(j));
        end
        criteria = N_events>=minNevents & minwake<maxwakethres & FREM<=maxFREM ; %& (~isnan(Pos)&supine==1)
        Nwindows(n)=sum(criteria);
        if sum(criteria)==0
            continue
        end
        subplot(ceil(sqrt(M)),ceil(M/ceil(sqrt(M))),n);
        plot(1,1,'marker','none'); hold('on');
        title(filename_{n});
        hold('on');
        
        [t,x,y,a,win,nota,veup,hyp,pcw,t2,e] = VEPesArray(DataOut,n,criteria,drivecol);
        if dividedrivebyttot
            x=x./t2;
        end
        [ad]=howfarawayfromsleep(a,win); %replace a with (a==1)&(hyp==4) if arousals not scored in wake...
        
%         if 1 && drivecol==16
%             x = x+2*pcw; %equivalent ccw is 39 or 7 based on mixed model analysis of all breaths dPes vs dEdi
%             %xeup = Veupnea/nanmedian(y(a==1)'./x(a==1));
%         end
        
        %normalize data
        g = G(n);
        
        if normalizeyusingknownVeupnea
            y_ = Veupnea_(n)/60;
            y = y./y_*100;
        else
            if normalizeyusingveupwindow
                y = y./veup*100; %note Y is unnormalized VE from LGfromflow
                y_ = nanmedian(veup);
            else
                y_ = nanmedian(veup);
                y = y./y_*100;
            end
        end
        
        if (Fvol>0||Fpcw>0) && drivecol==16
            x = x + x_*Fvol*(y/y_/100) + Fpcw*pcw; 
        end
        
        if 0
            x_ = y_/g; %xeup
            x_n(n)=x_;
            y_n(n)=y_;
        else %use y_ from first section
            x_ = y_n(n)/g; %xeup
            x_n(n)=x_;  
        end
        
        if 0
            Gtemp=interp1(Gtime{n}(:,1),Gtime{n}(:,2),t);
            if isnan(Gtemp(1))
                I=find(~isnan(Gtemp),1,'first');
                Gtemp(1:I-1)=Gtemp(I);
            end
            if isnan(Gtemp(end))
                I=find(~isnan(flipud(Gtemp)),1,'first');
                Gtemp((end-I+2):end)=Gtemp(end-I+1);
            end
            if 0
            g=nanmedian(Gtime{n}(:,2));
            x=x./(Gtemp/g);
            x_=y_n(n)/g;
            x_n(n)=x_;  
            else
            x=x./(Gtemp/g);   
            end
        end
        %Arousal threshold
        
        %First find arthres for Vactive/Vcomp:
        settings2=[3 0 1 0 1 0 1 0.2]; 
            %1:ignorefirstXbreathsaftersleeponset [4]
            %2:swapbreathsforhigherdrive [N]
            %3:Nbreathsattributabletoarousal [N]
            %4:increasingdriveonly [N]
            %5:deleteifbeloweupnea [Y]
            %6:setaseupneaifbeloweupnea [N]
            
                tempx=x; tempx(hyp~=2)=NaN;
                [artN1(n),artN1_N(n)] = ArThresNew(a',tempx',x_,win,settings2);
                tempx=x; tempx(hyp~=1)=NaN;
                [artN2(n),artN2_N(n)] = ArThresNew(a',tempx',x_,win,settings2);
                tempx=x; tempx(hyp~=0)=NaN;
                [artN3(n),artN3_N(n)] = ArThresNew(a',tempx',x_,win,settings2);
        
        tempx=x;
        if 0 %drive data not in events are deleted
             tempx(e==1)=NaN;
        end
        [arthres2(n),arthres2_N(n)] = ArThresNew(a',tempx',x_*1.05,win,settings2);
        
        settings2=[3 0 1 0 1 0 1 0.2]; 
        %settings2=[4 0 1 0 1 0 1 0.2]; %5 better for Ar, 4 better for Vcomp
        [arthres(n),arthres_N(n)] = ArThresNew(a',tempx',x_*1.05,win,settings2); %,'minwake',minwake,60
        
        
        
        %For arousal threshold as per modelpsg method
        
        if 0 %turn this off
            temparthres=NaN*ones(max(win),1);
            temphyp=NaN*ones(max(win),1);
            ArOnset = [NaN;diff(a)];  
            %elast=[NaN;e(1:end-1)];  
            for i=1:max(win)     
                temparthres(i)=nanmean(x(ArOnset==1&win==i)); %&elast==0
                temphyp(i)=mode(hyp(win==i));
            end
            arthres(n)=nanmedian(temparthres);
            arthres2(n)=arthres(n);
        end
        
        %For arousal threshold (Vactive/Vcomp)
        if 1 %turn this back on
            if arthres2(n)<(1.05*x_),arthres2(n)=1.05*x_; end
            if arthres(n)<(1.00*x_),arthres(n)=1.00*x_; end
        end
        
        %GG analysis
        try
            [~,~,gg,ggtonic] = GGPesArray(DataOut,n,criteria,drivecol);
        catch me
            gg=NaN*y;
            ggtonic=NaN*y;
        end
        if 1
            a1=(ad>2)&(a==1);
            %lower=0.75; upper=1.25;
            %tempggdata=gg((a1)&(x>lower*x_)&(x<upper*x_)&(y>lower*100)&(y<upper*100));
            tempggdata=gg(a1);
            GGpeakwake(n)=nanmedian(tempggdata);
            GGwakeN(n)=length(tempggdata);
        end
        
        %Remove wakefulness breaths for further analysis
        if 1 %scoredarousalsinwake
            x = x(nota==1);
            y = y(nota==1);
            gg = gg(nota==1);
            ggtonic = ggtonic(nota==1);
        else
            x = x(nota==1&hyp~=4);
            y = y(nota==1&hyp~=4);
            gg = gg(nota==1&hyp~=4);
            ggtonic = ggtonic(nota==1&hyp~=4);
        end
        
        medianV(n) = nanmedian(y);
        
        ploton=1;
        figure(3);
        subplot(ceil(sqrt(M)),ceil(M/ceil(sqrt(M))),n);
        hold('on');
        %plot gridlines
        %plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
        %plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
        plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
        plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
        
        [Vpassive(n),Vactive(n),VEdeciles,Vdrivedeciles,VpassiveCI(n,:),VactiveCI(n,:)]=VEVdriveAnalysis(x(~isnan(y)&~isnan(x))/x_*100,y(~isnan(y)&~isnan(x)),arthres2(n)/x_*100,ploton,x_/x_*100,Nciles);
        set(gca,'fontsize',fontsize_)
        xlim([0 max([prctile(x/x_*100,95) 200])]);
        ylim([0 max([prctile(y,95) 150])]);
        if n==1
            ylabel('Ventilation'); 
        end
        if n==M
            xlabel('Vdrive');
        end
        hold('off');
        
        figure(33)
        subplot(ceil(sqrt(M)),ceil(M/ceil(sqrt(M))),n);
        plot(1,1,'marker','none'); hold('on');
        plot([1 1],[0 150],'--','color',[0.7 0.7 0.7]);
        hold('on');
        plot([0 3*1],[100 100],'--','color',[0.7 0.7 0.7]);
        h=plot(x/x_*100,y,'k.','markersize',2); box('off');
        
        try %GG analysis
            figure(4);
            subplot(ceil(M/ceil(sqrt(M))),ceil(sqrt(M)),n);
            hold('on');
            %plot gridlines
            plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
            
            %plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
            if 0
                [GGpassive(n),GGactive(n),GGdeciles,Vdrivedeciles]=VEVdriveAnalysis(x(~isnan(gg)&~isnan(x)),gg(~isnan(gg)&~isnan(x)),arthres(n),ploton,x_,Nciles);
                hold('off');
                set(gca,'fontsize',fontsize_)
                xlim([0 max([prctile(x,95) 0.001])]);
                ylim([0 max([prctile(gg,95) 0.001])]/GGpeakwake(n)*100);
                ylabel('GGpeak, %max'); xlabel('Pes, cmH2O');
            else
                [GGpassive(n),GGactive(n),GGdeciles,Vdrivedeciles]=VEVdriveAnalysis(x(~isnan(gg)&~isnan(x)),gg(~isnan(gg)&~isnan(x))/GGpeakwake(n)*100,arthres(n),ploton,x_,Nciles);
                hold('off');
                set(gca,'fontsize',fontsize_)
                xlim([0 max([prctile(x,95) 0.001])]);
                ylim([0 max([prctile(gg,95) 0.001])]/GGpeakwake(n)*100);
                ylabel('GGpeak, %wake'); xlabel('Pes, cmH2O');
            end
            
            figure(5);
            subplot(ceil(M/ceil(sqrt(M))),ceil(sqrt(M)),n);
            hold('on');
            %plot gridlines
            plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
            
            %plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
            [GGtonicpassive(n),GGtonicactive(n),GGtonicdeciles,Vdrivedeciles]=VEVdriveAnalysis(x(~isnan(ggtonic)&~isnan(x)),ggtonic(~isnan(gg)&~isnan(x)),arthres(n),ploton,x_,Nciles);
            hold('off');
            set(gca,'fontsize',fontsize_)
            xlim([0 max([prctile(x,95) 0.01])]);
            ylim([0 max([prctile(gg,95) 0.01])]);
            ylabel('GGtonic, %max'); xlabel('Pes, cmH2O');
        catch me
            %disp('no GG analysis');
        end
        %title(['ahi:',num2str(round(BaselineAHI(n)),'%u') '->' num2str(round(TreatmentAHI(n)),'%u') ',LGn:' num2str(LGn(n),2)],'fontsize',8);
        
    catch me
        ['error:' num2str(n)] 
        Vpassive(n)=NaN;
        Vactive(n)=NaN;
        medianV(n)=NaN;
        arthres(n)=NaN;
    end
    
end
Vcomp = Vactive-Vpassive;

    %Vcomp(Vpassive>90&Vactive>90&Vcomp<10&Vcomp>-10)=NaN
    %VcompPSG(VpassivePSG>95&VactivePSG>95&VcompPSG<10&VcompPSG>-10)=NaN
%end

% Study-specific Comparisons

%Vcomp(Vactive==0|Vpassive==0)=NaN;
Vcomp_ = Vactive_-Vpassive_;
if 1
    Vcomp_(Vactive_<=-100|Vpassive_<=-100)=NaN;
    Vactive_(Vactive_<-100)=-100;
    Vpassive_(Vpassive_<-150)=-150;
end

% 95%CIs
temp=abs(VpassiveCI-[Vpassive Vpassive]);
temp=mean(temp')';
[nanmean(temp),nanstd(temp)]/1.96

temp2=abs(VactiveCI-[Vactive Vactive]);
temp2=mean(temp2')';
[nanmean(temp2),nanstd(temp2)]/1.96

[nanmean(Nwindows),nanstd(Nwindows)]

%% Plot comparisons A
plotwithnumbers=0;

corrtype = 'Pearson';% 'Pearson' 'Spearman' %Pearson generally stronger

figure(10);
subplot(1,3,1); 
x=Vpassive; y=VpassivePSG;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 0
    plotregressionwithSEM(x,y); hold('on'); 
   [r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
   %[r,p,pftest,xmodelpassive,xmodelpassive_sem]=plotregressionwithSEM_powerlaw(x,y);
   [r,p,pftest,xmodelpassive,xmodelpassive_sem]=plotregressionwithSEM_sqrt(x,y);
end
title({'Vpassive';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylabel('psg model');
xlabel('psg actual edi/pes');
ylim([-1 107]);
xlim([-1 100]);

subplot(1,3,2); 
x=Vactive; y=VactivePSG;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p,slope,int,sem]=plotregressionwithSEM(x,y); hold('on');      
hold('on'); %plot([0 1],[0 1],'-','color',[0.5 0.5 0.5]);
 %[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vactive';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylim([-3 115]);
xlim([-3 115]);

subplot(1,3,3); 
x=Vcomp; y=VcompPSG; 
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 1
[r,p,slope,int,sem]=plotregressionwithSEM(x,y); hold('on');      
else
[r,p,slope,p2]=plotregressionwithSEM_linearsegment(x,y); hold('on');
end
hold('on'); %plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vcomp';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
xlim([-70 70]);
ylim([-60 40]);

%% Plot Arthres
figure(11)
subplot(1,2,1); 
x=arthres./x_n*100; y=arthresPSG*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 1
[r,p,arthresslope]=plotregressionwithSEM(x,y); hold('on'); 
% [r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
[r,p,slope,p2]=plotregressionwithSEM_linearsegment(x,y); hold('on');
end
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
% xlim([97 370]);
 ylim([97 199]);

%thres = prctile(x,100/3);
clear data
thres = prctile(x,100/2);
[thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(x<thres,y)
data(1,:) = [thres,thresX,AUC,SEM,p,sensspec]

thres = prctile(x,100/3);
[thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(x<thres,y)
data(2,:) = [thres,thresX,AUC,SEM,p,sensspec]

thres = prctile(x,100/4);
[thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(x<thres,y)
data(3,:) = [thres,thresX,AUC,SEM,p,sensspec]

thres = prctile(x,100/5);
[thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(x<thres,y)
data(4,:) = [thres,thresX,AUC,SEM,p,sensspec]

performance105 = PredictiveValue(x<prctile(x,100/4),y<105,x)
performance110 = PredictiveValue(x<prctile(x,100/2),y<110,x)
performance115 = PredictiveValue(x<prctile(x,100/2),y<115,x)
performance120 = PredictiveValue(x<prctile(x,100/2),y<120,x)

 if 0
subplot(1,2,2); 
x=arthres; y=arthresPSG*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 1
plotregressionwithSEM(x,y); hold('on'); 
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
[r,p,slope,p2]=plotregressionwithSEM_linearsegment(x,y); hold('on');
end
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
%xlim([97 370]);
%ylim([97 190]);
 end
[nanmean(arthres_N) nanstd(arthres_N)]



if drivecol==16
    PesarthresFeupnea = x;
elseif drivecol==19
    EdiarthresFeupnea = x;
end
if exist('EdiarthresFeupnea')&&exist('PesarthresFeupnea')
    figure(61)
    [r,p,slope]=plotregressionwithSEM(EdiarthresFeupnea,PesarthresFeupnea); 
end
    
%% Plot C
figure(12)
subplot(1,3,1); Vpassive
x=Vpassive_; y=Vpassive;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 0
plotregressionwithSEM(x,y); hold('on'); 
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
[r,p,slope,p2,~,~,~,~,model,sem]=plotregressionwithSEM_linearsegment(x,y); hold('on');
end
title({'Vpassive';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylabel('psg actual edi/pes');
xlabel('cpap drops');
ylim([-5 105]);
xlim([-150 100]);
set(gca,'ytick',0:20:100)
set(gca,'Position',[0.209829867674858 0.11 0.196597353497164 0.790793653490052]);

subplot(1,3,2); 
x=Vactive_; y=Vactive; 
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 0
plotregressionwithSEM(x,y); hold('on'); 
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
lowerlimitonupperbreakpoint=60; %needed to match sigmaplot outcome
[r,p,slope,p2,~,~,~,~,model,sem]=plotregressionwithSEM_linearsegment(x,y,'lowerlimitonupperbreakpoint',lowerlimitonupperbreakpoint); hold('on');
end
title({'Vactive';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
set(gca,'ytick',0:20:100)
ylim([-5 105]);
xlim([-100 210]);
set(gca,'Position',[0.456836798991808 0.11 0.227473219911783 0.790793653490052]);

subplot(1,3,3); 
x=Vcomp_; y=Vcomp;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p,slope,int,sem]=plotregressionwithSEM(x,y); hold('on'); 
title({'Vcomp';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
% subplot(3,4,4); plot(arthresgs,arthres./x_n*100,'.','markersize',14); box('off'); hold('on'); plot([1 3],[1 3],'-','color',[0.5 0.5 0.5]);
% x=arthresgs; y=arthres./x_n*100; [r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
% title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
% if plotwithnumbers, plotnumbers(x,y,DPW); end
ylim([-70 70]);
xlim([-80 170]);
set(gca,'Position',[0.730938878386894 0.11 0.174061121613107 0.790793653490052]);

%% Plot D
figure(13)
subplot(1,3,1); 
x=Vpassive_; y=VpassivePSG;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 0
plotregressionwithSEM(x,y); hold('on'); 
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
[r,p,VpassiveSlope,p2,~,~,~,~,xmodel,sem]=plotregressionwithSEM_linearsegment(x,y); hold('on');
end
hold('on'); %plot([0 1],[0 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vpassive';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylabel('psg model');
xlabel('cpap drops');
xlim([-150 100]);
ylim([-5 115]);
set(gca,'Position',[0.209829867674858 0.11 0.196597353497164 0.790793653490052]);
set(gca,'ytick',0:20:100)

subplot(1,3,2); 
x=Vactive_; y=VactivePSG;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
if 0
plotregressionwithSEM(x,y); hold('on'); 
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
else
[r,p,VactiveSlope,p2,~,~,~,~,xmodel,sem]=plotregressionwithSEM_linearsegment(x,y); hold('on');
end
hold('on'); plot([0 1],[0 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vactive';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
xlim([-100 210]);
ylim([-5 115]);
set(gca,'Position',[0.456836798991808 0.11 0.227473219911783 0.790793653490052]);


subplot(1,3,3); 
x=Vcomp_; y=VcompPSG; 
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p,slope,int,sem]=plotregressionwithSEM(x,y);
hold('on'); %plot([0 1],[0 1],'-','color',[0.5 0.5 0.5]);
title({'Vcomp';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylim([-53 35]);
xlim([-80 170]);
set(gca,'Position',[0.730938878386894 0.11 0.174061121613107 0.790793653490052]);

%% save plots
savefigdirectory='C:\Users\szs88\Dropbox (Personal)\Phenotyping collapsbility and responsiveness from PSG\';
saveas(10,[savefigdirectory '10.fig']);
saveas(11,[savefigdirectory '11.fig']);
saveas(12,[savefigdirectory '12.fig']);
saveas(13,[savefigdirectory '13.fig']);

%%
return
%%
AAtemp=[Vpassive_ Vpassive VpassivePSG Vactive_ Vactive VactivePSG];
%% Veupnea 
figure(80)
plot(Veupnea_,y_n*60,'.','markersize',15); box('off');
xlabel('Veupnea, CPAP');
ylabel('mean ventilation'); %assumes method used above...
plotnumbers(Veupnea_,y_n*60,DPW);
[r,p]=plotregressionwithSEM(Veupnea_,y_n*60); hold('on'); 

nanstd((y_n*60-Veupnea_)./Veupnea_*100)
temp=Veupnea_-y_n*60;
temp(14)=NaN;

nanmean(temp./Veupnea_*100)
nanstd((temp)./Veupnea_*100);



%% GGcomp vs Vcomp
figure(11)

corrtype = 'Pearson';

GGactive(GGactive==0)=NaN;
GGpassive(GGpassive==0)=NaN;
GGcomp = GGactive-GGpassive;
GGcompp = (GGactive-GGpassive)./GGactive*100;

xdata = GGcompp; %%%%%%%%%%%%%%

x=xdata(:); y=Vcomp_(:);
subplot(1,3,1); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vcomp';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylabel('CPAP drops');
xlabel('GGcomp');
x=xdata(:); y=Vcomp(:);
subplot(1,3,2); plot(x,y,'.','markersize',14); title('Vcomp'); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vcomp';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylabel('psg edi/pes');
xlabel('GGcomp');
x=xdata(:); y=VcompPSG(:);
subplot(1,3,3); plot(x,y,'.','markersize',14); title('Vcomp'); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vcomp';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end
ylabel('psg model');

xlabel('GGcomp');

%% mixed model analysis of Pes = (xEdi +jVol) | subject

pN=[];
eN=[];
vN=[];
sN=[];
%[16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak]
minNevents=1
for n=1:M
    try
        N_events = LG_QualityInfo{n}(:,2);
        %Pos = LG_QualityInfo{n}(:,5);
        FREM = SleepData{n}(:,6);
        minwake = SleepData{n}(:,7);
        %criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=maxwakethres;
        criteria = N_events>=minNevents & minwake<=maxwakethres & FREM<=maxFREM;
        
        drivecol=16;
        [t,x,y,a,win,nota,veup,hyp,pcw] = VEPesArray(DataOut,n,criteria,drivecol);
        drivecol=19;
        [~,edi] = VEPesArray(DataOut,n,criteria,drivecol);
        vol = pcw/5;
        edihigh = nanmedian(edi); vlow = prctile(vol(edi>edihigh),0.05);
        ratiouvpercmH2O(n) = nanmedian(edi(vol<vlow&edi>edihigh)./x(vol<vlow&edi>edihigh));
        subject=0*x+n;
        
        pN=[pN;x];
        sN=[sN;subject];
        vN=[vN;vol];
        eN=[eN;edi/ratiouvpercmH2O(n)*2];
    catch me
    end
end

%%
eN_ = eN.^1 - 0*vN;
pN_ = pN.^1 + vN*0; % + vN*30/0.77
vN_ = vN.^1;
Tbl=table(sN,eN_,vN_,pN_, ...
    'VariableNames',{'Subjects','Edi','Vol','Pes'});

Tbl.Subjects = nominal(Tbl.Subjects);
%         Tbl.Edi = ordinal(Tbl.Edi);
%         Tbl.Vol = ordinal(Tbl.Vol);
%         Tbl.Pes = ordinal(Tbl.Pes);
%lme = fitlme(Tbl,['Pes ~ 1 + Edi + Vol + (Edi | Subjects) '])

%lme = fitlme(Tbl,['Edi ~ 1 + Pes + Vol + (Vol | Subjects)'])
lme = fitlme(Tbl,['Edi ~ 1 + Pes + Vol + (Pes | Subjects) + (Vol | Subjects) '])
%lme = fitlme(Tbl,['Vol ~ 1 + Edi + Pes + (Edi| Subjects)'])

lme.Rsquared

%     figure(1)
%     plot(pN_,eN_,'.','markersize',2);



%% nonlinear regression for Vpassive. Sometimes not as good as Sigmaplot nonlinear regression...
addpath('E:\Work\Projects SS\Phenotyping and Edi');
global fixedslope upperlimit
fixedslope=NaN;
upperlimit.on=1;
figure(99)
[yintsem,yint,xintsem,xint,slope,rsquared]=plotregressionwithSEM_linearsegment(Vpassive_',Vpassive');

tempx=Vpassive_(~isnan(Vpassive_)&~isnan(Vpassive))
tempy=Vpassive(~isnan(Vpassive_)&~isnan(Vpassive))
ypred = tempx*slope+yint;
hold('on');
plot(tempx,ypred,'g.')
Fres = sum((tempy-ypred).^2)
Ftot = sum((tempy-mean(tempy)).^2)
rsq = 1-Fres/Ftot;

[r,p]=corr(Vpassive_(~isnan(Vpassive_)&~isnan(Vpassive)),Vpassive(~isnan(Vpassive_)&~isnan(Vpassive)))

[p,S,mu] = polyfit(Vpassive_(~isnan(Vpassive_)&~isnan(Vpassive)),Vpassive(~isnan(Vpassive_)&~isnan(Vpassive)),1)
xtemp=[min(Vpassive_) max(Vpassive_)]
[ytemp,deltatemp] = polyval(p,xtemp,S,mu);

ypred2= polyval(p,tempx,S,mu);
Fres2 = sum((tempy-ypred2).^2);
rsq2 = 1-Fres2/Ftot;
hold('on');
plot(xtemp,ytemp,'k:')
r^2

%% Does Vpassive just correlate with Vpassive_ because Vactive and Vpassive are related?

figure(55)


x=Vpassive_(:); y=Vpassive(:);
subplot(1,6,1); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vpassive, Edi vs CPAP';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end

x=Vactive(:); y=Vactive_(:);
subplot(1,6,2); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'Vactive, Edi vs CPAP';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end

x=Vpassive(:); y=Vactive(:);
subplot(1,6,3); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'VactiveEdi vs VpassiveEdi';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end

x=Vpassive_(:); y=Vactive_(:);
subplot(1,6,4); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'VactiveCPAP vs VpassiveCPAP';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end

x=Vpassive_(:); y=Vactive(:);
subplot(1,6,5); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'VactiveEdi vs VpassiveCPAP';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end

x=Vpassive(:); y=Vactive_(:);
subplot(1,6,6); plot(x,y,'.','markersize',14); box('off'); hold('on'); plot([-0.5 1],[-0.5 1],'-','color',[0.5 0.5 0.5]);
[r,p]=corr(x(~isnan(x)&~isnan(y)),y(~isnan(x)&~isnan(y)),'type',corrtype);
title({'VactiveCPAP vs VpassiveEdi';['r=',num2str(r,2),', p=',num2str(p,2)]});
if plotwithnumbers, plotnumbers(x,y,DPW); end


%% UA phenotype using model drive in REM (under construction)
unnormalize = 0;

fontsize_=8;

VpassivePSG_REM=NaN*zeros(M,1);
VactivePSG_REM=NaN*zeros(M,1);
medianVPSG_REM=NaN*zeros(M,1);
arthresPSG_REM=NaN*zeros(M,1);
arthresPSG2_REM=NaN*zeros(M,1);
minNwindows = 1;
minFREM=0.33;
minNevents=0;
maxwakethres=360;
Nciles=10;

for n=1:M
    if Exclude(n)==1||LG1_N(n)<=minNwindows
        continue
    end
    figure(2);
    
    %title(filename_{n});
    try
        N_events = LG_QualityInfo{n}(:,2);
        %Pos = LG_QualityInfo{n}(:,5);
        FREM = SleepData{n}(:,6);
        minwake = SleepData{n}(:,7);
        %criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=maxwakethres;
        criteria = N_events>=minNevents & minwake<=maxwakethres & FREM>minFREM;
        sum(criteria);
        if sum(criteria)==0
            continue
        end
        subplot(ceil(sqrt(M)),ceil(M/ceil(sqrt(M))),n);
        
        plot(1,1,'marker','none'); hold('on');
        title(filename_{n});
        hold('on');
        [t x y a win nota veup hyp] = VEVdriveArray(DataOut,n,criteria);
        
        
        
        
        
        x_=1;
        if unnormalize
            y = y.*veup;
            x = x.*veup;
            ymean = mean(y); %re-normalize with mean over whole night
            y = y/ymean;
            x = x/ymean;
        end
        y=y*100;
        
        
        %arthres(n) = ArThresFromPSG(1,x,win);
        %increasingdriveonly, deleteifbeloweupnea(0),
        %setaseupneaifbeloweupnea(1), usemediannotROC(1), falsepositiverate(0.2)
        temparthres=NaN*ones(max(win),1);
        ArOnset = [NaN;diff(a)];  
        for i=1:max(win)     
            temparthres(i)=nanmean(x(ArOnset==1&win==i));
        end
        arthresPSG(n)=nanmedian(temparthres);
        
        if 1 %set minimum for arthres PER SE
            if arthresPSG(n)<1.00
                arthresPSG(n)=1.00;
            end
        end
        arthresPSG2(n) = arthresPSG(n);
        if arthresPSG2(n)<1.05
            arthresPSG2(n)=1.05;
        end

        
        x = x(nota==1&hyp==3);
        y = y(nota==1&hyp==3);
        
        
        medianV(n) = median(y);
        
        ploton=1;
        
        plot([1 1],[0 150],'--','color',[0.7 0.7 0.7]);
        hold('on');
        plot([0 3*1],[100 100],'--','color',[0.7 0.7 0.7]);
        [Vpassive_REM(n),VactivePSG_REM(n)]=VEVdriveAnalysis(x,y,arthresPSG2(n),ploton,1,Nciles);
        set(gca,'fontsize',fontsize_)
        ylim([0 max([prctile(x,95) 1])]);
        ylim([0 max([prctile(y,95) 1])]);
        ylabel('Ventilation, %eupnea'); xlabel('Vdrive model, %eupnea');
        hold('off');
        
        
        %title(['ahi:',num2str(round(BaselineAHI(n)),'%u') '->' num2str(round(TreatmentAHI(n)),'%u') ',LGn:' num2str(LGn(n),2)],'fontsize',8);
    catch me
        VpassivePSG_REM(n)=NaN;
        VactivePSG_REM(n)=NaN;
        medianVPSG_REM(n)=NaN;
        arthresPSG_REM(n)=NaN;
    end
end

VcompPSG_REM = VactivePSG_REM-VpassivePSG_REM;

%% Load arthres Phenotype PTSS2013
addpath('E:\Work\MatlabFunctionDatabase');
corrtype='Pearson'
if 0
dirandfile='E:\Work\Projects SS\LG from Polysomnography\PT+SS\Writing\phenotype_parameters_20131005.xlsx'
sheet='20130205'
art_CPAPrange = 'BR24:BR54';
art_PSGrange = 'BZ24:BZ54';
art_CPAP=xlsread(dirandfile,sheet,art_CPAPrange);
art_PSG=xlsread(dirandfile,sheet,art_PSGrange);
art_PSG = (art_PSG+1)*100;
art_CPAP = (1+art_CPAP)*100
elseif 1
dirandfile='E:\Work\Projects SS\LG from Polysomnography\PT+SS\Writing\ERJ\Pheno_O2_ACZ_LGresults_20140322_flippedflowfix_old_SS2017.xlsx'
sheet='Pheno PSGSummaryData - 100pc'
art_CPAPrange = 'AC40:BG40';
art_PSGrange = 'AC35:BG35';
art_CPAP=xlsread(dirandfile,sheet,art_CPAPrange);
art_PSG=xlsread(dirandfile,sheet,art_PSGrange);
art_PSG = (art_PSG)*100;
art_CPAP = (1+art_CPAP)*100
elseif 0
dirandfile='E:\Work\Projects SS\LG from Polysomnography\PT+SS\Writing\ERJ\Pheno_O2_ACZ_LGresults_20140322_flippedflowfix_old_SS2017.xlsx'
sheet='Pheno PSGSummaryData - 50pc'
art_CPAPrange = 'AC57:BG57';
art_PSGrange = 'AC52:BG52';
art_CPAP=xlsread(dirandfile,sheet,art_CPAPrange);
art_PSG=xlsread(dirandfile,sheet,art_PSGrange);
art_PSG = (art_PSG)*100;
art_CPAP = (1+art_CPAP)*100
end

figure(19)
subplot(1,1,1); 
x=art_CPAP; y=art_PSG;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
% if plotwithnumbers, plotnumbers(x,y,DPW); end
xlim([97 300]);
ylim([97 230]);

set(gca,'fontname','arial narrow','fontsize',12)


%% MORE

%% Arthres Edi/Pes versus sleep state

include = 1==ones(1,length(artN3))' %~isnan(artN3)&artN3~=0
Nincl = sum(include)

arstates = [artN1(include)./x_n(include) artN2(include)./x_n(include) artN3(include)./x_n(include)]*100;

[artN1_N(include)' artN2_N(include)' artN3_N(include)']

nanmean(arstates)

arstates2 = [artN1(include)./x_n(include);artN2(include)./x_n(include);artN3(include)./x_n(include)]*100;
subjects = [1:Nincl 1:Nincl 1:Nincl]';
states = [1+0*(1:Nincl) 2+0*(1:Nincl) 3+0*(1:Nincl)]';

%stats, mixed model:
Tbl=table(subjects,states-2,arstates2, ...
    'VariableNames',{'Subjects','State','ARthres'});
Tbl.Subjects = nominal(Tbl.Subjects);
lme = fitlme(Tbl,['ARthres ~ 1 + State + (1 | Subjects)'])
lme.Rsquared;

%% ArPSG versus sleep state [trend for reduction]

include=1==ones(1,length(artN3PSG))' %~isnan(artN3PSG)&artN3PSG~=0
Nincl = sum(include)

arstates = [artN1PSG(include) artN2PSG(include) artN3PSG(include)]*100;

[artN1PSG_Nwin(include) artN2PSG_Nwin(include) artN3PSG_Nwin(include)]

nanmean(arstates)
arstates2 = [artN1PSG(include);artN2PSG(include);artN3PSG(include)]*100;
subjects = [1:Nincl 1:Nincl 1:Nincl]';
states = [1+0*(1:Nincl) 2+0*(1:Nincl) 3+0*(1:Nincl)]';

%stats, mixed model:
Tbl=table(subjects,states,arstates2, ...
    'VariableNames',{'Subjects','State','ARthres'});
Tbl.Subjects = nominal(Tbl.Subjects);
lmePSG = fitlme(Tbl,['ARthres ~ 1 + State + (1 | Subjects)'])
lmePSG.Rsquared;

%%


AHInrem
Fhypops
NadirNREMSpO2


%% vs AHI, SpO2, Fhypops [arthres arthresPSG]

figure(111)
subplot(3,3,1); 
x=AHInrem; y=arthresPSG*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
subplot(3,3,1+3); 
x=AHInrem; y=arthres./x_n*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
subplot(3,3,1+6); 
x=AHInrem; y=arthres;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(3,3,2); 
x=NadirNREMSpO2; y=arthresPSG*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
subplot(3,3,2+3); 
x=NadirNREMSpO2; y=arthres./x_n*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
subplot(3,3,2+6); 
x=NadirNREMSpO2; y=arthres;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(3,3,3); 
x=Fhypops; y=arthresPSG*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
subplot(3,3,3+3); 
x=Fhypops; y=arthres./x_n*100;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});
subplot(3,3,3+6); 
x=Fhypops; y=arthres;
removenan=isnan(x)|isnan(y); x(removenan)=[]; y(removenan)=[]; 
[r,p]=plotregressionwithSEM(x,y); hold('on'); 
title({'Ar.Thres';['r=',num2str(r,2),', p=',num2str(p,2)]});



aadata = [arthres./x_n*100 arthresPSG*100 AHInrem NadirNREMSpO2 Fhypops]


%%

minAHIbaseline=5;
deltaAHIp_N1toN2 = (AHInrem1-AHInrem2)./AHInrem1*100; deltaAHIp_N1toN2(AHInrem1<minAHIbaseline)=NaN;
deltaAHI_N1toN2 = (AHInrem1-AHInrem2); deltaAHI_N1toN2(AHInrem1<minAHIbaseline)=NaN;
deltaAHIp_N2toN3 = (AHInrem2-AHInrem3)./AHInrem2*100; deltaAHIp_N2toN3(AHInrem2<minAHIbaseline)=NaN;
deltaAHI_N2toN3 = (AHInrem2-AHInrem3); deltaAHI_N2toN3(AHInrem2<minAHIbaseline)=NaN;



figure(111)
subplot(2,2,1)
[r,p]=plotregressionwithSEM(artN1PSG*100,deltaAHIp_N1toN2);
title({'DeltaAHI N2-N1 % vs ArThresPSGN1';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(2,2,2)
[r,p]=plotregressionwithSEM(artN1PSG*100,deltaAHI_N1toN2);
title({'DeltaAHI N2-N1 vs ArThresPSGN1';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(2,2,3)
[r,p]=plotregressionwithSEM(artN1./x_n*100,deltaAHIp_N1toN2);
title({'DeltaAHI N2-N1 % vs ArThresN1';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(2,2,4)
[r,p]=plotregressionwithSEM(artN1./x_n*100,deltaAHI_N1toN2);
title({'DeltaAHI N2-N1 vs ArThresN1';['r=',num2str(r,2),', p=',num2str(p,2)]});


figure(112)
subplot(2,1,1)
[r,p]=plotregressionwithSEM(AHInrem1,deltaAHIp_N1toN2);
title({'DeltaAHI N2-N1 % vs ArThresPSGN1';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(2,1,2)
[r,p]=plotregressionwithSEM(AHInrem1,deltaAHI_N1toN2);
title({'DeltaAHI N2-N1 vs ArThresPSGN1';['r=',num2str(r,2),', p=',num2str(p,2)]});



figure(113)
subplot(2,1,1)
[r,p]=plotregressionwithSEM(LG1,deltaAHIp_N1toN2);
title({'DeltaAHI N2-N1 % vs ArThresPSGN1';['r=',num2str(r,2),', p=',num2str(p,2)]});

subplot(2,1,2)
[r,p]=plotregressionwithSEM(LG1,deltaAHI_N1toN2);
title({'DeltaAHI N2-N1 vs ArThresPSGN1';['r=',num2str(r,2),', p=',num2str(p,2)]});



save 123321.mat