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
    %%
    load('mdlREM')
    %%
    clear FullTable
    %%
    
    plotfig1=1;
    %%
    % fix 220, R2=-Inf
    % currently using time when SpO2 is on
    for n=25%1:2:Npatients% (970 no pbetafilt1,982 data matrix X)
        %try
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
            
            dT = DataEventHypnog_Mat_ds.Time(2)-DataEventHypnog_Mat_ds.Time(1);
            
            DataEventHypnog_Mat_ds.Epochs(DataEventHypnog_Mat_ds.Epochs<0|DataEventHypnog_Mat_ds.Epochs>4)=NaN;
            %
            
%             
%             DataEventHypnog_Mat_ds.WakeNoAR = (DataEventHypnog_Mat_ds.Epochs==4)*1;
%             %DataEventHypnog_Mat_ds.WakeNoAR(isnan(DataEventHypnog_Mat_ds.Epochs))=NaN;
%             DataEventHypnog_Mat_ds.ExcludeAR = zeros(length(DataEventHypnog_Mat_ds.WakeNoAR),1);
%             if ~noscoredarinwake
%                 DataEventHypnog_Mat_ds.ExcludeAR(DataEventHypnog_Mat_ds.Epochs==4&DataEventHypnog_Mat_ds.EventsAr<0.5)=1;
%             end
%             
%             DataEventHypnog_Mat_ds.ExcludeAR(DataEventHypnog_Mat_ds.Epochs~=4&DataEventHypnog_Mat_ds.EventsAr>0.5)=1;
%             DataEventHypnog_Mat_ds.ExcludeAR = logical(DataEventHypnog_Mat_ds.ExcludeAR);
%             
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
            
            Exclude = SpO2off==1 | isnan(DataEventHypnog_Mat_ds.Epochs) | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>4;
            
            %for judging methods
            ExcludeREM = SpO2off | DataEventHypnog_Mat_ds.EventsAr==1 | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>3 | isnan(DataEventHypnog_Mat_ds.Epochs);
            
            REMvsNREM = (DataEventHypnog_Mat_ds.Epochs==3)*1;
            REMvsNREM(ExcludeREM)=NaN;
            
            %%
            
            if ~noscoredarinwake
                ExcludeAR = Exclude | DataEventHypnog_Mat_ds.Epochs==4 & DataEventHypnog_Mat_ds.EventsAr==0 | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
            else
                ExcludeAR = Exclude | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
            end
            
            dsfactor = round(3/dT);
            I = (1:dsfactor:length(Exclude))';
            clear WSpredlogit acc
            for j=1:8
                try
                    [~,acc(j)] = WSfromoneEEG(Exclude(I),mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds(I,:),j,3,ExcludeAR(I));
                catch me
                    
                end
            end
            [~,accmaxi] = max(acc);
            
            [WSpredlogit,Acc,Exclude_,NoiseBinary] = WSfromoneEEG(Exclude,mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds,accmaxi,dT,ExcludeAR);
            %Exclude is updated here to include noise, caused trouble e.g.
            %for n=53 RICCADSA
            
           %% 
            tic
            Options = [];
            %ArSignal = DataEventHypnog_Mat_ds.EventsAr;
            %ArSignal = logitinverse(WSpredlogit)>0.25;
            ArSignal = 1*(DataEventHypnog_Mat_ds.EventsAr==1 | logitinverse(WSpredlogit)>0.75);
            
            widthrms=120;
            BufferI = buffer(1:length(WSpredlogit),round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
                    BufferI(:,sum(BufferI==0)>0)=[];       
            t = DataEventHypnog_Mat_ds.Time(BufferI);% buffer(DataEventHypnog_Mat_ds.Time,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
            t_ = t(1,:)' + widthrms/2;
            WSslow = nanmedian(WSpredlogit(BufferI))';
            WSslow = interp1(t_,WSslow,DataEventHypnog_Mat_ds.Time,'linear');
            
            ArSignal = 1*(DataEventHypnog_Mat_ds.EventsAr==1 | logitinverse(WSpredlogit)>0.75) | WSslow>0.75;
            
            [PrREM,Tbl,EOGrmsConj_p50,ymtaAr180b,t_,rmssumbyrbyar,rmssumbyr,rmssum,Thres,LOCfilt,ROCfilt] = REMfromEOG(DataEventHypnog_Mat_ds.LOC,DataEventHypnog_Mat_ds.ROC,DataEventHypnog_Mat_ds.Time,ArSignal,Exclude,Options,mdlREM);
            if 0
                PrREM(NoiseBinary==1)=NaN;
            end
            toc
            
            % Plot
            if 1
                Tbl.Time = DataEventHypnog_Mat_ds.Time;
                Tbl.Subj = ones(length(REMvsNREM),1)*n+Exp*10000;
                Tbl.REMvsNREM = REMvsNREM;
                Tbl.ExcludeREM = ExcludeREM;
                Tbl.WSslowPr = logitinverse(WSslow);
                Tbl.WSslowLogOdds = WSslow;
            end
            
            ThresF = 5;
            TonicREM_ = DataEventHypnog_Mat_ds.Epochs==3 & rmssumbyr<ThresF*(10.^Thres) & rmssum<(ThresF*50*(10.^Thres));
            PhasicREM_ = DataEventHypnog_Mat_ds.Epochs==3 & rmssumbyr>ThresF*(10.^Thres);
            
            TonicREM = RemoveShortSegments(TonicREM_,6,dT,0);
            PhasicREM = 1-RemoveShortSegments(1-PhasicREM_,2,dT,0);
            
            if plotfig1
                figure(1); clf(1); set(gcf,'color',[1 1 1]);
                ax(1)=subplot(4,1,1);
                plot(DataEventHypnog_Mat_ds.Time,[DataEventHypnog_Mat_ds.LOC LOCfilt]);
                set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                ax(2)=subplot(4,1,2);
                plot(DataEventHypnog_Mat_ds.Time,[DataEventHypnog_Mat_ds.ROC ROCfilt]);
                set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                ax(3)=subplot(4,1,3);
                plot(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.Epochs);
                hold on
                plot(DataEventHypnog_Mat_ds.Time,[DataEventHypnog_Mat_ds.EventsAr] + 4.5); 
                plot(DataEventHypnog_Mat_ds.Time,logitinverse(WSpredlogit) + 4.5);
                plot(DataEventHypnog_Mat_ds.Time,logitinverse(WSslow) + 4.5,'color','r');
                plot(DataEventHypnog_Mat_ds.Time,[TonicREM_+5.9 TonicREM+6 PhasicREM_+7.4 PhasicREM+7.5]);
                plot(DataEventHypnog_Mat_ds.Time,[ArSignal*1+8.75]);
                plot(Tbl.Time,PrREM,'linewidth',1.5);
                set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                ax(4)=subplot(4,1,4);
                
                plot(DataEventHypnog_Mat_ds.Time,[rmssum.^0.5 -rmssum.^0.5],'color',[0.6 0.6 0.6]);
                hold on
                plot(DataEventHypnog_Mat_ds.Time,[-rmssumbyr.^0.5],'color',[0.2 0.2 0.2]);
                plot(DataEventHypnog_Mat_ds.Time,[rmssumbyrbyar.^0.5],'color',[0.9 0.1 0.2]);
                plot(DataEventHypnog_Mat_ds.Time,(10.^[Thres ymtaAr180b]).^0.5,'color',[0.4 0.3 0.2]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                set(gca,'box','off','tickdir','out');
                %         figure(8)
                %         plot(-1:0.01:1,logitinverse(((-1:0.01:1)-0.5)*10))
                linkaxes(ax,'x');
            end
            
            %Tbl = table(Subj,Time,REMvsNREM,ExcludeREM,p50,p5_p50,p95_p50,p25_p50,p75_p50,EOGrmsConj_p50_90s,EOGrmsConj_p50_30s,EOGrmsConj_p50_120s,EOGrmsConj_p50_180s,EOGrmsConj_p50_90sb,EOGrmsConj_p50_120sb,EOGrmsConj_p50_180sb);
            Tbl_ = Tbl(1:round(3/dT):end,:);                      
            mdlREMa = fitglm(Tbl_,'REMvsNREM ~ EOGrmsConj_p50','Distribution','binomial','Exclude',isnan(Tbl_.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
            mdlREMa.Rsquared.Ordinary
            
            
            if 0
                if exist('FullTable')
                    FullTable = [FullTable;Tbl_];
                    %FullTable = outerjoin(FullTable,Tbl,'MergeKeys',true);
                    %BreathDataFullTable = [BreathDataFullTable;BreathDataTable2];
                else
                    FullTable = [Tbl_]; %start new table
                end
            else
                'warning: not building fulltable'
            end
            
        %catch me1
        %    me1.message
        %end
    end
end
%%

%%
I = isnan(FullTable.REMvsNREM);
FullTable(I,:)=[];
%I = (FullTable.Subj==40004);
%FullTable(I,:)=[];

mdlREM = fitglm(FullTable,'REMvsNREM ~ EOGrmsConj_p50 + p95_p50 + p5_p50 + p75_p50 + p25_p50','Distribution','binomial','Exclude',isnan(FullTable.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
%mdlREM = fitglm(FullTable,'REMvsNREM ~ EOGrmsConj_p50 + p95_p50','Distribution','binomial','Exclude',isnan(FullTable.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
mdlREM.Rsquared.Ordinary
mdlREM=compact(mdlREM);
save mdlREM mdlREM -v7.3
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
constantslist2 = {'p5_p50','p95_p50','p25_p50','p75_p50'};
%constantslist2 = {'p95_p50'};
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
    terms = terms - 1/sum(coeffs1)*coeffsconstant(i)*eval(['FullTable.' constantnames{i}]);
end
FullTable.ref = FullTable.p50 - terms; %not tested
FullTable.EOGrmsConj_ref = FullTable.EOGrmsConj_p50 - terms; %good

mdlREMref = fitglm(FullTable,'REMvsNREM ~ EOGrmsConj_ref','Distribution','binomial','Exclude',isnan(FullTable.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
mdlREMref.Rsquared.Ordinary
mdlREM.Rsquared.Ordinary

%%
FullTable.EOGrmsConj_refT = (10.^FullTable.EOGrmsConj_ref).^0.1;
mdlREMrefT = fitglm(FullTable,'REMvsNREM ~ EOGrmsConj_refT','Distribution','binomial','Exclude',isnan(FullTable.REMvsNREM)) % +pred1_p50_+pred1_p75_+pred1_p95_
mdlREMrefT.Rsquared.Ordinary


REMpred = predict(mdlREMrefT,FullTable);
performance = PredictiveValue(1*(FullTable.REMvsNREM>0.5),1*(REMpred>0.5),FullTable.REMvsNREM)
performance.TP_FP_TN_FN
[x,y,t,AUC,~] = perfcurve(1*(FullTable.REMvsNREM(~isnan(FullTable.REMvsNREM))>0.5),REMpred(~isnan(FullTable.REMvsNREM)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))

%%
save FullTable FullTable -v7.3