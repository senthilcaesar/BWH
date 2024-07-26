%% Start
clear all

logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p));
normalci = @(p,n) 1.96*((p.*(1-p)./n).^0.5);

sqrttrans = @(x) (1-(1-x).^0.5);
sqrttransi = @(x) (1-(1-x).^2);

global settings
settings.CurrentCodeVersion = 'PUPbeta20190629';
mydir  = mfilename('fullpath'); %only works under "run" not "evaluate line"
idcs   = strfind(mydir,filesep);
settings.workdir = mydir(1:idcs(end-1));
settings.codedir = ['C:\Users\szs88\Dropbox (Personal)\PUPbeta_git\' settings.CurrentCodeVersion '\'];

addpath(genpath('G:\Dropbox (Personal)\PUPbeta_git\PUPbeta20190629\'));
addpath('G:\Dropbox (Personal)\PhenotypeDrive2018\Workspaces');

%% Only if full rerun; do not clear if combining cohorts
%%
plotson=1;
plottraces=1;
dT = 3;
NoisePredThres=0.50;
MinNoiseTime=60;
savedatalargetable=1;
savedatalargetablebest=1;
useSpO2off=1; %zero will ignore SpO2 in final run; turning off [zet to zero) makes little difference (tested in Phenotype, higher median acc, lower mean acc 0.5%)

%default
%%
for Exp=4
    %Exp=4
    switch Exp
        case 1
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\PhenotypeDrive2018\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Personal)\PhenotypeDrive2018\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1)
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [patients{i,2} patients{i,1}]
            end
            settings.savename='PhenoDrive2019';
            noscoredarinwake=0; %0 is have scoring of arousals in wake

        case 2
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Partners HealthCare)\MAD-OX\Traits\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Partners HealthCare)\MAD-OX\Traits\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            idcs   = strfind(settings.AMasterdir,filesep);
            settings.workdir = settings.AMasterdir(1:idcs(end-1));
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1)
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [settings.workdir 'Converted\' patients{i,1}]
            end
            settings.savename='MADOX2019';
            noscoredarinwake=0; %0 is have scoring of arousals in wake

        case 3
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\SaraODB OAT\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Personal)\SaraODB OAT\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            idcs   = strfind(settings.AMasterdir,filesep);
            settings.workdir = settings.AMasterdir(1:idcs(end-1));
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1)
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [settings.workdir 'Converted\' patients{i,1}]
            end
            settings.savename='SaraOBDOAT2019';
            noscoredarinwake=0; %0 is have scoring of arousals in wake

        case 4 %C3-A2, C4-A1, O1-A2, Cz-A1
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            idcs   = strfind(settings.AMasterdir,filesep);
            settings.workdir = settings.AMasterdir(1:idcs(end-1));
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1);
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [settings.workdir 'Converted\' patients{i,1}];
            end
            MultiScorer = xlsread(AMasterSpreadsheet,1,'R4:R10003');
            filedir(MultiScorer>0)=[];
            Npatients = size(filedir,1)
            settings.savename='RICCADSA2019';
            noscoredarinwake=0; %0 is have scoring of arousals in wake
        case 5
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\MESA_Converted_Partial\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Personal)\MESA_Converted_Partial\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            idcs   = strfind(settings.AMasterdir,filesep);
            settings.workdir = settings.AMasterdir(1:idcs(end-1));
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1)
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [settings.workdir patients{i,1}];
            end
            filedir
            settings.savename='MESA2019';
            noscoredarinwake=1;
        case 6
            settings.AMasterdir = 'J:\PEOPLE\FACULTY\SANDS\O2PSG\TraitsAnalysis\SSO2Rerun2019\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            idcs   = strfind(settings.AMasterdir,filesep);
            settings.workdir = settings.AMasterdir(1:idcs(end-1));
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1)
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [settings.workdir 'Converted\' patients{i,1}];
            end
            filedir
            settings.savename='SSO2PSGreun2019';
            noscoredarinwake=0;
        case 9 %RICCADSAMultiscoring
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\PUPstart\';
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            idcs   = strfind(settings.AMasterdir,filesep);
            settings.workdir = settings.AMasterdir(1:idcs(end-1));
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1);
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [settings.workdir 'Converted\' patients{i,1}];
            end
            MultiScorer = xlsread(AMasterSpreadsheet,1,'R4:R10003');
            Iexcl = (MultiScorer==0|MultiScorer==99);
            filedir(Iexcl)=[];
            MultiScorerID = xlsread(AMasterSpreadsheet,1,'q4:q10003');
            MultiScorerID(Iexcl)=[];
            settings.savename='RICCADSAMultiscoring';
            noscoredarinwake=0; %0 is have scoring of arousals in wake
            Npatients = size(filedir,1)
    end
    
    %%
    addpath(genpath(settings.workdir));
    addpath(genpath(settings.codedir));
    addpath(genpath(pwd));
    
    %%
    
    %load('workspaceXYZ.mat', 'RefTable', 'mdlA');
    %save('RefTable','RefTable');
    %save('mdlA','mdlA')
    load('mdlA');
    load('RefTable');
    load('mdlAcc');
    load('mdlNoise');
    load('DecilesLookupRef','DecilesLookupRef')
    
%     save RefTable RefTable
%     save mdlA mdlA
    
    %% Loop to make giant data table all patients all breaths, nonoverlapping
    
    
    
   
    %15,17 is missing Pbetalogfilt1
    %47 is missing SpO2, fix this
    %% 
    clear FullTable
    %% 
    
    plotfig1=1;
    %%
    %currently using time when SpO2 is on
    for n=1:Npatients% (970 no pbetafilt1,982 data matrix X)
        try
        %%
        %pause
        disp(['start: ' num2str(n)]);
        loadpath=[filedir{n}];
        
        try
            load(loadpath);
        catch me
            disp(me.message);
            'no data, skipped'
            continue
        end
        
        dt = DataEventHypnog_Mat(2,1) - DataEventHypnog_Mat(1,1);
        Fs = 1./dt;
        %%
        dT = dt;
        
        
        
        %%
        dsfactor = dT/dt;
        DataEventHypnog_Mat_ds = DataEventHypnog_Mat(1:round(dsfactor):end,:);
        DataEventHypnog_Mat_ds = array2table(DataEventHypnog_Mat_ds);
        DataEventHypnog_Mat_ds.Properties.VariableNames = ChannelsList;
        
        
        DataEventHypnog_Mat_ds.Epochs(DataEventHypnog_Mat_ds.Epochs<0|DataEventHypnog_Mat_ds.Epochs>4)=NaN;
        %
        
        
        DataEventHypnog_Mat_ds.WakeNoAR = (DataEventHypnog_Mat_ds.Epochs==4)*1;
            %DataEventHypnog_Mat_ds.WakeNoAR(isnan(DataEventHypnog_Mat_ds.Epochs))=NaN;
        DataEventHypnog_Mat_ds.ExcludeAR = zeros(length(DataEventHypnog_Mat_ds.WakeNoAR),1);
        if ~noscoredarinwake
        DataEventHypnog_Mat_ds.ExcludeAR(DataEventHypnog_Mat_ds.Epochs==4&DataEventHypnog_Mat_ds.EventsAr<0.5)=1;
        end
        DataEventHypnog_Mat_ds.ExcludeAR(DataEventHypnog_Mat_ds.Epochs~=4&DataEventHypnog_Mat_ds.EventsAr>0.5)=1;
        
        DataEventHypnog_Mat_ds.ExcludeAR = logical(DataEventHypnog_Mat_ds.ExcludeAR);
        
        if ~ismember(['Pbetalogfilt1'],DataEventHypnog_Mat_ds.Properties.VariableNames)
            disp('missing Pbetalogfilt1');
            continue
        end
        
        try
            SpO2off = 1*(DataEventHypnog_Mat_ds.SpO2<50);
            SpO2off = RemoveShortSegments(SpO2off,60,dT);
            
            SpO2starti = find(DataEventHypnog_Mat_ds.SpO2>50,1,'first');
            SpO2endi = find(DataEventHypnog_Mat_ds.SpO2>50,1,'last');
            SpO2ever = 1*SpO2off; SpO2ever(SpO2starti:SpO2endi)=0; SpO2ever=logical(SpO2ever);
        catch me
            disp('missing SpO2 signal, skipping');
            continue
        end
        
        if sum(SpO2ever==1)>0.90*length(SpO2ever)
            disp('missing SpO2 signal, skipping');
            continue
        end
        
        
        %% 
        %DataEventHypnog_Mat_ds
        
        filter_HFcutoff_butter1 = 4; %4,0.3
        filter_LFcutoff_butter1 = 1;
        filter_order0 = 2;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dT/2));
        LOCfilt = nanfilter(B_butter0,A_butter0,DataEventHypnog_Mat_ds.LOC,1);
        ROCfilt = nanfilter(B_butter0,A_butter0,DataEventHypnog_Mat_ds.ROC,1);
        
        if 0
        LOCfilt(DataEventHypnog_Mat_ds.EventsAr==1)=NaN;
        ROCfilt(DataEventHypnog_Mat_ds.EventsAr==1)=NaN;
        end

        
        figure(1); clf(1); 
        ax(1)=subplot(4,1,1);
        plot(DataEventHypnog_Mat_ds.Time,[DataEventHypnog_Mat_ds.LOC LOCfilt]);
        ax(2)=subplot(4,1,2);
        plot(DataEventHypnog_Mat_ds.Time,[DataEventHypnog_Mat_ds.ROC ROCfilt]);
        ax(3)=subplot(4,1,3);
         plot(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.Epochs);
         hold on
         plot(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.EventsAr + 4.5);
        
        
        widthrms = 6;
        
        
        BufferI = buffer(1:length(DataEventHypnog_Mat_ds.Time),round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
        t = DataEventHypnog_Mat_ds.Time(BufferI);% buffer(DataEventHypnog_Mat_ds.Time,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
        t_ = t(1,:)' + widthrms/2;
        
        LOCbuffer = LOCfilt(BufferI); %buffer(LOCfilt,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay'); 
        LOCrms = (nanmean((LOCbuffer.^2)).^0.5)';
        
        ROCbuffer = ROCfilt(BufferI);% buffer(ROCfilt,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
        ROCrms = (nanmean((ROCbuffer.^2)).^0.5)'; %rms of yR
        
        if 0
        LOCrms=LOCrms.^2;
        ROCrms=ROCrms.^2;
        end
        
        %y2 = movmean(LOCfilt,[round(widthrms*1/dT/2) round(widthrms*1/dT/2)],'omitnan');  

        % later replace this with EEG signal

        %ar = interp1(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.EventsAr,t_,'nearest'); %far faster than  
        ar = DataEventHypnog_Mat_ds.EventsAr(BufferI); 
            ar = nanmean(ar)';

        
        
        clear rval
        rval = NaN*zeros(length(ROCrms),1);
        for i=1:length(ROCrms) %%%%%%%%% approx 5 sec
             rval(i)=corr(LOCbuffer(:,i),ROCbuffer(:,i));  
        end
        
        if 0
        figure(8)
        plot(-1:0.01:1,logitinverse(((-1:0.01:1)-0.5)*10))
        end
        
        rval2 = logitinverse((-rval-0.5)*10);
        rmssumbyr=(LOCrms + ROCrms).*rval2;
        rmssum=(LOCrms + ROCrms);
        rmssumbyrbyar=(LOCrms + ROCrms).*rval2.*(1-0.9999*ar);
        rmssumbyrnanar=rmssumbyr;
            rmssumbyrnanar(ar>0.5)=NaN;
        
            if 0
        rmssumbyrnanar=rmssumbyrnanar.^2;
        rmssumbyrnanar=rmssumbyrnanar.^2;
        rmssum=rmssumbyrnanar.^2;
            end
        
        refsig = log10(interp1(t_,rmssum,DataEventHypnog_Mat_ds.Time,'nearest'));
        
        
    
    if plotfig1
        figure(1)
        ax(1)=subplot(4,1,1);
        hold on
        plot(t_,LOCrms*2);
        ax(2)=subplot(4,1,2);
        hold on
        plot(t_,ROCrms*2);
        linkaxes(ax,'x')
    end
    
        widthmta = 90;
        
        BufferI2 = buffer(1:length(t_),round(widthmta*1/(t_(2)-t_(1))),round(widthmta*1/(t_(2)-t_(1)))-1,'nodelay');
        ymtaAr90 = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrnanar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
        ymtaAr90b = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrbyar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
        
        refsig = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssum(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
        %refsig = ymtaAr90; not as good
%         
        widthmta = 30;        
        BufferI2 = buffer(1:length(t_),round(widthmta*1/(t_(2)-t_(1))),round(widthmta*1/(t_(2)-t_(1)))-1,'nodelay');
        ymtaAr30 = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrnanar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));

        widthmta = 180;        
        BufferI2 = buffer(1:length(t_),round(widthmta*1/(t_(2)-t_(1))),round(widthmta*1/(t_(2)-t_(1)))-1,'nodelay');
        ymtaAr180 = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrnanar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
        ymtaAr180b = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrbyar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
        
        widthmta = 120;        
        BufferI2 = buffer(1:length(t_),round(widthmta*1/(t_(2)-t_(1))),round(widthmta*1/(t_(2)-t_(1)))-1,'nodelay');
        ymtaAr120 = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrnanar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
        ymtaAr120b = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrbyar(BufferI2)),DataEventHypnog_Mat_ds.Time,'nearest'));
%         t = buffer(DataEventHypnog_Mat_ds.Time,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
%            t(I,end)=NaN;
%            
%         t_ = (nanmean(t))';
        
           Exclude = SpO2off | DataEventHypnog_Mat_ds.Epochs>4 | isnan(DataEventHypnog_Mat_ds.Epochs);   
           ExcludeREM = SpO2off | DataEventHypnog_Mat_ds.EventsAr==1 | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>3 | isnan(DataEventHypnog_Mat_ds.Epochs);
        
        REMvsNREM = (DataEventHypnog_Mat_ds.Epochs==3)*1;       
        REMvsNREM(ExcludeREM)=NaN;
        balance = nansum(REMvsNREM)/sum(~isnan(REMvsNREM));
        
        
        
        p50 = prctile(refsig(~Exclude),50) * ones(length(REMvsNREM),1);
        

        
        EOGrmsConj_p50_90s = ymtaAr90 - p50;
        EOGrmsConj_p50_30s = ymtaAr30 - p50;
        EOGrmsConj_p50_180s = ymtaAr180 - p50;
        EOGrmsConj_p50_120s = ymtaAr120 - p50;
        
        EOGrmsConj_p50_180sb = ymtaAr180b - p50;
        EOGrmsConj_p50_120sb = ymtaAr120b - p50;
        EOGrmsConj_p50_90sb = ymtaAr90b - p50;

        p5 = prctile(refsig(~Exclude),5) * ones(length(REMvsNREM),1);
            p5_p50  = (p5 - p50);
        p95 = prctile(refsig(~Exclude),95) * ones(length(REMvsNREM),1);
            p95_p50 = (p95 - p50);
        p25 = prctile(refsig(~Exclude),25) * ones(length(REMvsNREM),1);
            p25_p50 = (p25 - p50);
        p75 = prctile(refsig(~Exclude),75) * ones(length(REMvsNREM),1);
            p75_p50 = (p75 - p50);
        
        Subj = ones(length(REMvsNREM),1)*n+Exp*10000;
        
        Time = DataEventHypnog_Mat_ds.Time;
        
        
        Tbl = table(Subj,Time,REMvsNREM,ExcludeREM,p50,p5_p50,p95_p50,p25_p50,p75_p50,EOGrmsConj_p50_90s,EOGrmsConj_p50_30s,EOGrmsConj_p50_120s,EOGrmsConj_p50_180s,EOGrmsConj_p50_90sb,EOGrmsConj_p50_120sb,EOGrmsConj_p50_180sb);
        Tbl = Tbl(1:round(3/dt):end,:);
        
        if 1
        weightsREM = zeros(length(REMvsNREM),1);
        weightsREM(REMvsNREM==1)=1-balance;
        weightsREM(REMvsNREM==0)=balance;
        weightsREM = weightsREM(1:round(3/dt):end);
        
            mdlREMa = fitglm(Tbl,'REMvsNREM ~ EOGrmsConj_p50_180sb','Distribution','binomial','weights',weightsREM,'Exclude',isnan(Tbl.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
            mdlREMa.Rsquared.Ordinary
            predout = predict(mdlREMa,Tbl);
            try
                predout1 = predict(mdlREM,Tbl);
                if plotfig1
                    subplot(4,1,3); hold('on'); plot(Tbl.Time,predout1,'linewidth',1.5)
                end
            catch me 
                'failed model mdlREM'
            end
%             subplot(4,1,3); hold('on'); plot(Tbl.Time,Tbl.rmta_)
            if plotfig1
                ax(4)=subplot(4,1,4); 
                hold on
                %plot(tmta_,log10(ymta_),'linewidth',2);
                if exist('mdlREM')
                    thres1 = p50 + -(mdlREM.Coefficients.Estimate(1) + p5_p50*mdlREM.Coefficients.Estimate(2)  + p95_p50*mdlREM.Coefficients.Estimate(3))/mdlREM.Coefficients.Estimate(4);
                else
                    thres1 = NaN*p50;
                end
                plot(t_,rmssumbyrnanar.^0.5); 
                plot(DataEventHypnog_Mat_ds.Time,(10.^[thres1 ymtaAr180b]).^0.5,'linewidth',1);
                
                 
%                 hold on
%                 plot(t_,[rmssum rmssumbyr rmssumbyrbyar],'linewidth',1);
              linkaxes(ax,'x');

            end
        end
        
        if 1
            if exist('FullTable')
                FullTable = [FullTable;Tbl];
                %FullTable = outerjoin(FullTable,Tbl,'MergeKeys',true);
                %BreathDataFullTable = [BreathDataFullTable;BreathDataTable2];
            else
                FullTable = [Tbl]; %start new table
            end
        end
        
        catch me1
           me1.message
        end
    end
end

%%
I = isnan(FullTable.REMvsNREM);
FullTable(I,:)=[];
%I = (FullTable.Subj==40004);
%FullTable(I,:)=[];


balance = nansum(FullTable.REMvsNREM)/sum(~isnan(FullTable.REMvsNREM));
balanceadjust = 3;        %1.7
        weightsREM = zeros(length(FullTable.REMvsNREM),1);
        weightsREM(FullTable.REMvsNREM==1)=1-balance*balanceadjust;
        weightsREM(FullTable.REMvsNREM==0)=balance*balanceadjust;
        sum(weightsREM(FullTable.REMvsNREM==1))
        sum(weightsREM(FullTable.REMvsNREM==0))
        
weightsREM = ones(length(FullTable.REMvsNREM),1);
        
%p5_p50,p10_p50,p25_p50,p33_p50,p67_p50,p75_p50,p90_p50,p95_p50,EOGrmsConj_p50
% + p33_p50 + p67_p50 + p90_p50
mdlREM = fitglm(FullTable,'REMvsNREM ~ EOGrmsConj_p50_180sb + p25_p50 + p75_p50 + p95_p50 + p5_p50','Distribution','binomial','weights',weightsREM,'Exclude',isnan(FullTable.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
        mdlREM.Rsquared.Ordinary
        
        REMpred = predict(mdlREM,FullTable);
 performance = PredictiveValue(1*(FullTable.REMvsNREM>0.5),1*(REMpred>0.5),FullTable.REMvsNREM)
 performance.TP_FP_TN_FN
        [x,y,t,AUC,~] = perfcurve(1*(FullTable.REMvsNREM(~isnan(FullTable.REMvsNREM))>0.5),REMpred(~isnan(FullTable.REMvsNREM)),1); %need to find the threshold value that gives the OPTROCPT!
        [~,I]=max(y+(1-x));
        thresFind=mean(t(I:(I+1)))
        
        
        
        %%
        if 1
            i=1
            uniquesubjlist = unique(FullTable.Subj);
            subj=uniquesubjlist(i);
            figure(2);clf(2);
            plot(1:length(FullTable.REMvsNREM(FullTable.Subj==subj)),FullTable.REMvsNREM(FullTable.Subj==subj));
            hold('on')
            plot(1:length(FullTable.REMvsNREM(FullTable.Subj==subj)),REMpred(FullTable.Subj==subj));
        end

       
%%
constantslist2 = {'p5_p50','p95_p50'};
constant = 1*(sum((string(mdlREM.Coefficients.Properties.RowNames) == string(constantslist2)),2)==1);
constant(1)=NaN;
I1 = find(constant==1);
constantnames = mdlREM.Coefficients.Properties.RowNames(I1);
coeffsconstant = mdlREM.Coefficients.Estimate(I1);
I2 = find(constant==0);
coeffs1 = mdlREM.Coefficients.Estimate(I2);
terms = 0;
temptemp=1/sum(coeffs1)*coeffsconstant;
for i=1:length(I1)
    terms = terms + 1/sum(coeffs1)*coeffsconstant(i)*eval(['FullTable.' constantnames{i}]);
end
FullTable.ref = FullTable.p50 + terms;


%%
save FullTable FullTable -v7.3