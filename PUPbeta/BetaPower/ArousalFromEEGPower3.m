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

%% Only if full rerun; do not clear if combining cohorts
if 1
clear Tacc BreathDataFullTable BreathDataFullTablebest 
end
%%
if 0
load('workspacetemp.mat')
end
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
for Exp=2
    %Exp=4
    switch Exp
        case 1
            settings.AMasterdir = 'C:\Users\szs88\Dropbox (Personal)\PhenotypeDrive2018\PUPstart\';
            settings.AMasterdir = 'G:\Dropbox (Personal)\PhenotypeDrive2018\PUPstart\';
    
            
            settings.savename='PhenoDrive2019';
            
            AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %
            [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
            Npatients = size(patients,1);
            clear filedir
            for i=1:Npatients
                filedir{i,1} = [patients{i,2} patients{i,1}];
            end
            
            clear filedirC
            parts = strsplit(settings.AMasterdir,'\');
            parent_path = [strjoin(parts(1:end-2),'\') '\'];
            for i=1:Npatients
                filedirC{i,1} = [parent_path 'Analyzed\' settings.savename '_' num2str(i)];
            end
            
            noscoredarinwake=0; %0 is have scoring of arousals in wake
            
            [~,~,raw] = xlsread(AMasterSpreadsheet,1,'AD4:AW10003');

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
    
    
    
    tic
    %15,17 is missing Pbetalogfilt1
    %47 is missing SpO2, fix this
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
        
        clear EEGRef WSpredlogit rareIDs100 noncentralness nonperipheralness Exclude acc
        for j=1:8
            try
                if ~ismember(['Pbetalogfilt' num2str(j)],DataEventHypnog_Mat_ds.Properties.VariableNames)
                    j=j-1;
                    break
                end
                
                
                Tselect=table();
                Tselect.Pbeta = eval(['DataEventHypnog_Mat_ds.Pbetalogfilt' num2str(j)]);
                Tselect.Palpha = eval(['DataEventHypnog_Mat_ds.Palphalogfilt' num2str(j)]);
                Tselect.Ptheta = eval(['DataEventHypnog_Mat_ds.Pthetalogfilt' num2str(j)]);
                Tselect.Pdelta = eval(['DataEventHypnog_Mat_ds.Pdeltalogfilt' num2str(j)]);
                
                Exclude(:,j) = SpO2off==1 | Tselect.Pbeta==-Inf | isnan(DataEventHypnog_Mat_ds.Epochs) | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>4;
                %before: Exclude = SpO2off==1 | isnan(DataEventHypnog_Mat_ds.Epochs);
                
                [WSpredlogit(:,j),Tselect,EEGRef(j,1)] = PredWakeSleep(Tselect,mdlA,RefTable,Exclude(:,j));
                Tselect.Total = Tselect.Pbeta_ref + Tselect.Palpha_ref + Tselect.Ptheta_ref + Tselect.Pdelta_ref;
                Tselect.NoiseBinary = RemoveShortSegments(1*(predict(mdlNoise,Tselect)>NoisePredThres),MinNoiseTime,dT);
                %[EEGRef(j,1),Tselect] = EEGreference(Tselect,RefTable,Exclude);
                %WSpredlogit(:,j) = logit(predict(mdlA,Tselect));
                
                %
                if plottraces
                    figure(j); clf(j);
                    
                    
                    ax1(1)=subplot(3,1,1);
                    plot(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.Epochs);
                    hold('on')
                    plot(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.EventsAr);
                    ax1(2)=subplot(3,1,2);
                    plot(DataEventHypnog_Mat_ds.Time,WSpredlogit(:,j));
                    ax1(3)=subplot(3,1,3);
                    plot(DataEventHypnog_Mat_ds.Time,DataEventHypnog_Mat_ds.SpO2);
                    
                    linkaxes(ax1,'x');
                end
                
                %generated noise model in here
                
                %noise based on j=last EEG signal only
                Tselect.Total = Tselect.Pbeta_ref + Tselect.Palpha_ref + Tselect.Ptheta_ref + Tselect.Pdelta_ref;
                Tselect.WSpredlogit = (WSpredlogit(:,j));
                Tselect.NoiseBinary = RemoveShortSegments(1*(predict(mdlNoise,Tselect)>NoisePredThres),MinNoiseTime,dT);
%                 
%                 NoisePred = predict(mdlNoise,Tselect);
%                 NoiseBinary = 1*(NoisePred>NoisePredThres);
%                 NoiseBinary = RemoveShortSegments(NoiseBinary,MinNoiseTime,dT);
                if plottraces
                    ax1(2)=subplot(3,1,2);
                    hold('on')
                    %plot(DataEventHypnog_Mat_ds.Time,NoisePred*10,'k');
                    plot(DataEventHypnog_Mat_ds.Time,Tselect.NoiseBinary*10,'g');
                    plot(DataEventHypnog_Mat_ds.Time,Tselect.Total,'g');
                end
                
                Exclude(:,j) = Tselect.NoiseBinary==1 | useSpO2off*(SpO2off==1) | Tselect.Pbeta==-Inf | isnan(DataEventHypnog_Mat_ds.Epochs) | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>4;
                [WSpredlogit(:,j),Tselect,EEGRef(j,1)] = PredWakeSleep(Tselect,mdlA,RefTable,Exclude(:,j));
                
                if ~noscoredarinwake
                    ExcludeAR = Exclude(:,j) | DataEventHypnog_Mat_ds.Epochs==4 & DataEventHypnog_Mat_ds.EventsAr==0 | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
                else
                    ExcludeAR = Exclude(:,j) | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
                end
                
                performance = PredictiveValue(1*(DataEventHypnog_Mat_ds.WakeNoAR(~ExcludeAR)>0.5),1*(WSpredlogit(~ExcludeAR,j)>0),DataEventHypnog_Mat_ds.WakeNoAR(~ExcludeAR));
                acc(j,1)=performance.Acc_sem_chance_p(1);
                
                
                %Exclude = SpO2off | Tselect.NoiseBinary>0 | Tselect.Pbeta_ref == -Inf;% | SpO2off==1;
                %Exclude = SpO2ever | isnan(DataEventHypnog_Mat_ds.Epochs) | NoiseBinary>0;% | SpO2off==1;
%                 [EEGRef(j,1),Tselect] = EEGreference(Tselect,RefTable,Exclude);
%                 
                if 0 %unnecessary complexity, not very helpful
                %ID method to assess rarity of power combinations
                IDa = [sum(Tselect.Pdelta_ref>DecilesLookupRef.PdeltaCut,2), ...
                    sum(Tselect.Ptheta_ref>DecilesLookupRef.PthetaCut,2), ...
                    sum(Tselect.Palpha_ref>DecilesLookupRef.PalphaCut,2), ...
                    sum(Tselect.Pbeta_ref>DecilesLookupRef.PbetaCut,2)];
                ID = IDa(:,1)*1000 + IDa(:,2)*100 + IDa(:,3)*10 + IDa(:,4)*1;
                Tselect.IDlikelihood = DecilesLookupRef.N(ID+1)/sum(DecilesLookupRef.N)*10000;
                
                %             figure(300);
                %             hist(Tselect.IDlikelihood,1000)
                rareIDs0p1(j,1)=sum(Tselect.IDlikelihood<0.1)/length(Tselect.IDlikelihood);
                rareIDs1(j,1)=sum(Tselect.IDlikelihood<1)/length(Tselect.IDlikelihood);
                
                end
                
                %WSpredlogit(:,j) = logit(predict(mdlA,Tselect));
                if plotson
                    ax1(2)=subplot(3,1,2);
                    plot(DataEventHypnog_Mat_ds.Time,WSpredlogit(:,j));
                    %plot(DataEventHypnog_Mat_ds.Time,Tselect.Total);
                    %plot(DataEventHypnog_Mat_ds.Time,Tselect.IDlikelihood*5-10);
                end
            catch me
                me.message
            end
        end
        
        
        clear maxcorr FextremeUnderNeg10 FextremeOver10 fitdata h1 meandiffZ tstatdiffWS tstatdiffWSadj bimodalitycoef diffhist pooledSD meandiff mean_ SD_ skewness_ kurtosis_ non5ness
        if plotson
            figure(99); clf(99); set(gcf,'color',[1 1 1]);
        end
        
        if 1
        TAccFeatures = AccFeatures(WSpredlogit,Exclude);
        PredAcc = predict(mdlAcc,TAccFeatures);
            PredAcc(PredAcc<0)=0; PredAcc(PredAcc>1)=1;
        end
        
        for i=1:j
            if plotson
                subplot(1,j,i)
            end
            
            dStep=0.25;
            Centers=-10:dStep:10;
            Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
            
            [h1(i,:),edges] = histcounts(WSpredlogit(~Exclude(:,i),i),Edges);
            h1(i,:) = h1(i,:)/sum(h1(i,:));
            
            %stairs(Centers-dStep/2,h1);
            %stairs(Centers-dStep/2,h2);
            if plotson
                bar(Centers,h1(i,:),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
                hold('on');
            end
            %tempexcl = (~DataEventHypnog_Mat_ds.ExcludeAR & ~Exclude);
            %performance = PredictiveValue(DataEventHypnog_Mat_ds.WakeNoAR(tempexcl)>0.5,WSpredlogit(tempexcl,i)>0,DataEventHypnog_Mat_ds.WakeNoAR(tempexcl));
            
            
            Iw = WSpredlogit(:,i)>0 & ~Exclude(:,i);
            Is = WSpredlogit(:,i)<=0 & ~Exclude(:,i);
            pooledSD(i,1) = (((sum(Iw)-1)*var(WSpredlogit(Iw,i)) + (sum(Is)-1)*var(WSpredlogit(Is,i)))/(sum(Iw)+sum(Is)))^0.5;
            meandiffZ(i,1) = (nanmean(WSpredlogit(Iw,i)) - nanmean(WSpredlogit(Is,i)))/pooledSD(i,1);
            meandiff(i,1) = (nanmean(WSpredlogit(Iw,i)) - nanmean(WSpredlogit(Is,i)));
            
            %%%
            mean_(i,1) = nanmean(WSpredlogit(~Exclude(:,i),i));
            SD_(i,1) = nanstd(WSpredlogit(~Exclude(:,i),i));
            
            skewness_(i,1) = skewness(WSpredlogit(~Exclude(:,i),i));
            kurtosis_(i,1) = kurtosis(WSpredlogit(~Exclude(:,i),i));
            FextremeOver10(i,1) = sum(WSpredlogit(~Exclude(:,i),i)>10)/sum(WSpredlogit(~Exclude(:,i),i)>-Inf);
            FextremeUnderNeg10(i,1) = sum(WSpredlogit(~Exclude(:,i),i)<-10)/sum(WSpredlogit(~Exclude(:,i),i)>-Inf);
            
            try
            ft = fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)');
            ftoptions = fitoptions(ft);
            ftoptions.StartPoint = [0.1 0.1 5 -5 3 3];
            ftoptions.Lower = [0 0 0  -10  0 0];
            ftoptions.Upper = [1 1 10  0 10 10];
            [fitobject,gof] = fit(Centers',h1(i,:)',ft,ftoptions);
            fitdata(i,:) = [fitobject.a1 fitobject.a2 fitobject.b1 fitobject.b2 fitobject.c1 fitobject.c2 gof.rsquare];
            if plotson
                plot(fitobject)
            end
            catch me
               fitdata(i,:) = [NaN NaN NaN NaN NaN NaN NaN];
            end
            
            %%%
            bimodalitycoef(i,1) = (skewness(WSpredlogit(~Exclude(:,i),i))^2 + 1)/kurtosis(WSpredlogit(~Exclude(:,i),i));
            
            %%%
            noncentralness(i,1) = (nanmean((WSpredlogit(~Exclude(:,i),i)).^2))^0.5;
            
            %%%
            temp = WSpredlogit(~Exclude(:,i),i);
            temp(temp<-10)=-10; temp(temp>10)=10;
            temp = min([abs(temp-10) abs(temp+10)]')';
            nonperipheralness(i,1) = (nanmean(temp.^2))^0.5;
            
            temp = WSpredlogit(~Exclude(:,i),i);
            temp = min([abs(temp-5) abs(temp+5)]')';
            non5ness(i,1) = (nanmean(temp.^2))^0.5;
            %         [~,~,~,stats] = ttest2(WSpredlogit(~Exclude,i)>0,WSpredlogit(~Exclude,i)<=0);
            %         tstatdiffWS(i) = stats.tstat;
            %         tstatdiffWSadj(i) = stats.tstat / sqrt(nansum(~Exclude));
            if plotson
                title(num2str(acc(i,1),2),'FontWeight','normal');
            end
            
            
            %%% correlations with others
            if j==1 %only one signal, thus can not compare to others
                maxcorr(i,1)=0;
            else
            clear corrtemp
            temp = WSpredlogit(:,i);
            temp(Exclude(:,i))=NaN;
            for x=1:j
                corrtemp(x)=NaN; %default
                if x==i
                    continue                    
                end
                
                temp2 = WSpredlogit(:,x);
                temp2(Exclude(:,x))=NaN;
                isnans = isnan(temp)|isnan(temp2);
                if sum(isnans)/length(isnans)==1
                    continue
                else
                    try
                    corrtemp(x) = corr(temp(~isnans),temp2(~isnans));
                    catch me
                    me.message
                    end
                end
                %error message: "Requires a data matrix X" - when?
                
            end
            tempmaxcorrtemp = max(corrtemp);
            if isnan(tempmaxcorrtemp)
                tempmaxcorrtemp=0;
            end
            maxcorr(i,1)=tempmaxcorrtemp; %%%
            end
        end
        
        clear diffhistmin diffhist
        for i=1:j
            diffhist(i,1) = sum(sum(abs(h1(i,:) - h1)));
        end
        for i=1:j %%%
            range=1:j; range(i)=[];
            %diffhistmin(i,1) = sum(min(abs(h1(range,:) - h1(i,:))));
            diffhistmin(i,1) = min(sum(abs(h1(range,:) - h1(i,:)),2));
        end
        
        subjn1 = zeros(length(acc),1) + n;
        Expn = zeros(length(acc),1) + Exp;
        subjn = Expn*10000 + subjn1;
        eegi = (1:j)';
        
        EEGRef_difffromtop = EEGRef - max(EEGRef);
        
        
        if ~exist('Tacc')
            Tacc = table();
            Tacc.acc=acc;
            Tacc.maxcorr=maxcorr;
            Tacc.subjn1=subjn1;
            Tacc.subjn=subjn;
            Tacc.Expn=Expn;
            Tacc.eegi=eegi;
            Tacc.meandiffZ=meandiffZ;
            Tacc.bimodalitycoef=bimodalitycoef;
            Tacc.EEGRef_difffromtop=EEGRef_difffromtop;
            Tacc.diffhist=diffhist;
            Tacc.noncentralness=noncentralness;
            Tacc.nonperipheralness=nonperipheralness;
            Tacc.meandiff=meandiff;
            Tacc.pooledSD=pooledSD;
            Tacc.mean_=mean_;
            Tacc.SD_=SD_;
            Tacc.skewness_=skewness_;
            Tacc.kurtosis_=kurtosis_;
            Tacc.diffhistmin=diffhistmin;
            Tacc.a1=fitdata(:,1);
            Tacc.a2=fitdata(:,2);
            Tacc.b1=fitdata(:,3);
            Tacc.b2=fitdata(:,4);
            Tacc.c1=fitdata(:,5);
            Tacc.c2=fitdata(:,6);
            Tacc.r2=fitdata(:,7);
            Tacc.FextremeOver10=FextremeOver10;
            Tacc.FextremeUnderNeg10=FextremeUnderNeg10;
            range = 1:size(Tacc,1);
        else
            I = Tacc.subjn==Exp*10000+n;
            Tacc(I,:)=[];
            range = size(Tacc,1)+1 : size(Tacc,1) + length(acc);
            Tacc{range,{'acc','maxcorr','subjn1','subjn','Expn','eegi','FextremeOver10','FextremeUnderNeg10','meandiffZ','bimodalitycoef','diffhist','EEGRef_difffromtop','nonperipheralness','noncentralness','pooledSD','meandiff','mean_','SD_','skewness_','kurtosis_','diffhistmin','a1','a2','b1','b2','c1','c2','r2'}}...
                       =[acc, maxcorr, subjn1,  subjn , Expn,  eegi , FextremeOver10, FextremeUnderNeg10, meandiffZ , bimodalitycoef , diffhist , EEGRef_difffromtop,  nonperipheralness,  noncentralness,  pooledSD,  meandiff,  mean_,  SD_,  skewness_, kurtosis_,   diffhistmin,fitdata(:,1),fitdata(:,2),fitdata(:,3),fitdata(:,4),fitdata(:,5),fitdata(:,6),fitdata(:,7)];
        end
        
        try
            Tacc.accPred(range,:)=predict(mdlAcc,Tacc(range,:));
            %Tacc.predAcc1(range,:) = predict(mdlAcc,Tacc(range,:));
            %accpred2 = predict(mdlAcc2,Tacc(range,:));
            TviewAcc = table(acc,Tacc.accPred(range,:));
            if plotson
                figure(99)
                
                for i=1:j
                    subplot(1,j,i)
                    xlabel(num2str(Tacc.accPred(range(i)),2));
                    ylabel('');
                    legend('off');
                    set(gca,'box','off','tickdir','out');
                end
                
                pause(0.2)
            end
            %pause
        catch me
            disp(me.message)
        end
        %     %%
        %     clear temp
        %     for i=1:size(BreathDataTable2,1)
        %         temp{i,1} = [settings.savename '_' num2str(n)];
        %     end
        %     BreathDataTable2.Subj = temp;
        %
        %     if exist('BreathDataFullTable')
        %         BreathDataFullTable = outerjoin(BreathDataFullTable,BreathDataTable2,'MergeKeys',true);
        %         %BreathDataFullTable = [BreathDataFullTable;BreathDataTable2];
        %     else
        %         BreathDataFullTable = [BreathDataTable2]; %start new table
        %     end
        %
        %     disp(['end: ' num2str(n)]);
        %     toc
        %pause(0.5)
        %%
        if savedatalargetable
            M=size(DataEventHypnog_Mat_ds,1);
                [~,bestEEGi] = max(Tacc.accPred(range));
                %[~,bestEEGi] = max(Tacc.acc(range));
                Time=DataEventHypnog_Mat_ds.Time;
                Pbeta=eval(['DataEventHypnog_Mat_ds.Pbetalogfilt' num2str(bestEEGi)]);
                Palpha=eval(['DataEventHypnog_Mat_ds.Palphalogfilt' num2str(bestEEGi)]);
                Ptheta=eval(['DataEventHypnog_Mat_ds.Pthetalogfilt' num2str(bestEEGi)]);
                Pdelta=eval(['DataEventHypnog_Mat_ds.Pdeltalogfilt' num2str(bestEEGi)]);
                %SpO2off;
                %SpO2ever;
                Epochs = DataEventHypnog_Mat_ds.Epochs;
                EventsAr = DataEventHypnog_Mat_ds.EventsAr;
                NoiseBinary = Tselect.NoiseBinary;                  % cheating here, should be separate noise for each EEG
                %Exclude = Exclude(:,bestEEGi);
                WakeNoAR = DataEventHypnog_Mat_ds.WakeNoAR;
                BestEEGi = ones(M,1)*bestEEGi;
                Subjn = ones(M,1)*subjn(1);
                
            BreathDataTable2 = table(Subjn,Time,BestEEGi,Epochs,EventsAr,NoiseBinary,WakeNoAR,SpO2off,SpO2ever,Pbeta,Palpha,Ptheta,Pdelta);
            
            if exist('BreathDataFullTable')
                BreathDataFullTable = outerjoin(BreathDataFullTable,BreathDataTable2,'MergeKeys',true);
                %BreathDataFullTable = [BreathDataFullTable;BreathDataTable2];
            else
                BreathDataFullTable = [BreathDataTable2]; %start new table
            end
        end
        if savedatalargetablebest
            M=size(DataEventHypnog_Mat_ds,1);
                %[~,bestEEGi] = max(Tacc.accPred(range));
                [~,bestEEGi] = max(Tacc.acc(range));
                Time=DataEventHypnog_Mat_ds.Time;
                Pbeta=eval(['DataEventHypnog_Mat_ds.Pbetalogfilt' num2str(bestEEGi)]);
                Palpha=eval(['DataEventHypnog_Mat_ds.Palphalogfilt' num2str(bestEEGi)]);
                Ptheta=eval(['DataEventHypnog_Mat_ds.Pthetalogfilt' num2str(bestEEGi)]);
                Pdelta=eval(['DataEventHypnog_Mat_ds.Pdeltalogfilt' num2str(bestEEGi)]);
                %SpO2off;
                %SpO2ever;
                Epochs = DataEventHypnog_Mat_ds.Epochs;
                EventsAr = DataEventHypnog_Mat_ds.EventsAr;
                NoiseBinary;
                %Exclude;
                WakeNoAR = DataEventHypnog_Mat_ds.WakeNoAR;
                BestEEGi = ones(M,1)*bestEEGi;
                Subjn = ones(M,1)*subjn(1);
                
            BreathDataTable2best = table(Subjn,Time,BestEEGi,Epochs,EventsAr,NoiseBinary,WakeNoAR,SpO2off,SpO2ever,Pbeta,Palpha,Ptheta,Pdelta);
            
            if exist('BreathDataFullTablebest')
                BreathDataFullTablebest = outerjoin(BreathDataFullTablebest,BreathDataTable2best,'MergeKeys',true);
                %BreathDataFullTable = [BreathDataFullTable;BreathDataTable2];
            else
                BreathDataFullTablebest = [BreathDataTable2best]; %start new table
            end
        end
        
        catch me1
            me1.message
        end
    end
end

%% test ability to detect best signal from statistics
%Tacc.subjnnom = nominal(Tacc.subjn);
Tacc.eeginom = nominal(Tacc.eegi);

NpatientsTotal = length(unique(Tacc.subjn))
if 1
    Itrain = (mod(Tacc.subjn,2)==1); %odd = training set
    Itest = ~Itrain;
else
    Itrain = logical(ones(size(Tacc,1),1));
    Itest = Itrain;
end

%PredVars = {'a1','a2','b1','b2','c1','c2','r2','EEGRef_difffromtop','diffhistmin','mean_','SD_','skewness_','kurtosis_','noncentralness','nonperipheralness'};

PredVars = {'maxcorr','bimodalitycoef','r2','diffhistmin','mean_','SD_','noncentralness','nonperipheralness'};

Tacc.acc(Tacc.acc<0.000001)=0.000001;
Tacc.acc(Tacc.acc>0.999999)=0.999999;
mdlAcc = fitglm(Tacc(Itrain,[PredVars {'acc'}]),'acc ~ nonperipheralness','ResponseVar','acc','link',1)

%first add good terms
Nsteps=30 %triple final number of terms, was equivalent to backwards elim
pthres = 0.5
mdlAcc = step(mdlAcc,'Criterion','sse','Penter',pthres,'Premove',pthres*1.01,'Nsteps',Nsteps,'upper','interactions')
mdlAcc.Rsquared.Ordinary
NtermsCurrent = length(mdlAcc.Coefficients.Estimate)-1

%reduce to 10 terms, can take a few loops
count=0;
while NtermsCurrent>10 || count>5
    count = count+1
NtermsCurrent = length(mdlAcc.Coefficients.Estimate)-1
Nsteps=NtermsCurrent-10%NtermsCurrent-50
pthres = 10^-200;
mdlAcc = step(mdlAcc,'Criterion','sse','Penter',pthres,'Premove',pthres*1.01,'Nsteps',Nsteps,'upper','interactions')
mdlAcc.Rsquared.Ordinary
NtermsCurrent = length(mdlAcc.Coefficients.Estimate)-1
end

%test set:
PredAtest = predict(mdlAcc,Tacc(Itest,:)); PredAtest(PredAtest<0)=0; PredAtest(PredAtest>1)=1;
I2 = ~isnan(PredAtest) & ~isnan(Tacc.acc(Itest));
temp = Tacc.acc(Itest);
SSE = nansum((temp(I2) - PredAtest(I2)).^2);
SStot = nansum((temp(I2) - nanmean(temp(I2))).^2);
Rsq = 1 - SSE/SStot

if 0
mdlAcc = compact(mdlAcc)
save mdlAcc mdlAcc -v7.3
end

%% Plot

Accpred1 = predict(mdlAcc,Tacc);
Accpred1(Accpred1>1)=1;
[Accpred1(find(Tacc.subjn==10008),:) Tacc.acc(find(Tacc.subjn==10008),:)]
[Accpred1(find(Tacc.subjn==20003),:) Tacc.acc(find(Tacc.subjn==20003),:)]
[Accpred1(find(Tacc.subjn==30029),:) Tacc.acc(find(Tacc.subjn==30029),:)]

size(Tacc,1)
% Plot Accuracy Prediction
figure(98); clf(98); set(gcf,'color',[1 1 1])
%hold('on')
xval = Accpred1;
yval = Tacc.acc;

yval = sqrttrans(yval)*100;
xval = sqrttrans(xval)*100;
scatter(xval,yval,4,[0.3 0.3 0.9],'filled','markerfacealpha',0.5);

set(gca,'xtick',100*sqrttrans([0 33 50 60 70 80 90 98 100]/100));
set(gca,'ytick',100*sqrttrans([0 33 50 60 70 80 90 98 100]/100));

axis([0 100 0 100])
set(gca,'xticklabels',sqrttransi(get(gca,'xtick')/100)*100)
set(gca,'yticklabels',sqrttransi(get(gca,'ytick')/100)*100)


% Plot Accuracy Prediction
figure(99); clf(99); set(gcf,'color',[1 1 1])
%hold('on')
xval = PredAtest;
yval = Tacc.acc(Itest);

yval = sqrttrans(yval)*100;
xval = sqrttrans(xval)*100;
scatter(xval,yval,10,[0.3 0.3 0.9],'filled','markerfacealpha',0.33);

set(gca,'xtick',100*sqrttrans([0 33 50 60 70 80 90 98 100]/100));
set(gca,'ytick',100*sqrttrans([0 33 50 60 70 80 90 98 100]/100));

temp = get(gca,'xtick')/100;

axis([0 100 0 100])
set(gca,'xticklabels',sqrttransi(temp)*100)
set(gca,'yticklabels',sqrttransi(get(gca,'ytick')/100)*100)


%%
%mdlAcc = compact(mdlAcc);
save('mdlAcc','mdlAcc','-v7.3')
save('Tacc','Tacc','-v7.3')

%% plot table data
%10033 fit is terrible

BreathDataFullTableBackup = BreathDataFullTable;
%%
% BreathDataFullTable

UniqueSubjList = unique(BreathDataFullTable.Subjn);

BreathDataFullTable.WakeNoAR=(BreathDataFullTable.Epochs==4)*1;


%%
if 0
    BreathDataFullTable(BreathDataFullTable.Subjn==30029,:)=[];
    BreathDataFullTableChosen(BreathDataFullTableChosen.Subjn==30029,:)=[];
    BreathDataFullTablebest(BreathDataFullTablebest.Subjn==30029,:)=[];
end
%%
BreathDataFullTableChosen = BreathDataFullTable;
BreathDataFullTablebest;

%%
% MADOX 92 acc, SODBOAT 93, RICCADSA 94 %
%% debug single patient data

%30029 is broken (events out of sync)
%20013 needs review
if 1
WStbl = BreathDataFullTableChosen;
else
WStbl = BreathDataFullTablebest;    
end

if 1
% e = [0 1] % [0 1 =0.9383] [0 1 =0.9382]
plotson=0;
for i=1:50%length(UniqueSubjList)%madox: 51:75
    
subjx = UniqueSubjList(i);

if plotson
figure(89); clf(89); set(gcf,'color',[1 1 1]);
end
I=WStbl.Subjn==subjx;

if plotson
ax89(1)=subplot(3,1,1);
plot(WStbl.Time(I),WStbl.Epochs(I),'r');
hold('on');
plot(WStbl.Time(I),WStbl.EventsAr(I));
set(gca,'box','off','xtick',[],'xcolor',[1 1 1])
ax89(2)=subplot(3,1,2);
hold('on');
plot(WStbl.Time(I),WStbl.SpO2off(I));

set(gca,'box','off','xtick',[],'xcolor',[1 1 1])
ax89(3)=subplot(3,1,3);
plot(WStbl.Time(I),WStbl.NoiseBinary(I)*5,'g');
hold('on');

box('off');
linkaxes(ax89,'x')
end
Tselect = WStbl(I,:);
% 
% ARscoringoff = ones(size(Tselect,1),1);
% li = find(Tselect.EventsAr>0,1,'first');
% ri = find(Tselect.EventsAr<1,1,'last');
% ARscoringoff(li:ri)=0;
% ARscoringoff=logical(ARscoringoff);

Exclude = Tselect.SpO2off==1 | Tselect.Pbeta==-Inf | isnan(Tselect.Epochs) | Tselect.Epochs<0 | Tselect.Epochs>4;
[WSpredlogit,Tselect,EEGRef] = PredWakeSleep(Tselect,mdlA,RefTable,Exclude);

if plotson
ax89(2)=subplot(3,1,2);
plot(WStbl.Time(I),Exclude,'r--');
ax89(3)=subplot(3,1,3);
plot(WStbl.Time(I),Tselect.WSpredlogit);
plot(WStbl.Time(I),Tselect.WakeNoAR);
end
Tselect.Total = Tselect.Pbeta_ref + Tselect.Palpha_ref + Tselect.Ptheta_ref + Tselect.Pdelta_ref;
Tselect.NoiseBinary = RemoveShortSegments(1*(predict(mdlNoise,Tselect)>NoisePredThres),MinNoiseTime,dT);

Exclude = e(1)&ARscoringoff | Tselect.NoiseBinary==1 | e(2)&Tselect.SpO2off==1 | Tselect.Pbeta==-Inf | isnan(Tselect.Epochs) | Tselect.Epochs<0 | Tselect.Epochs>4;
ExcludeAR = Exclude | Tselect.Epochs==4 & Tselect.EventsAr==0 | Tselect.Epochs~=4 & Tselect.EventsAr==1;

[WSpredlogit,Tselect,EEGRef] = PredWakeSleep(Tselect,mdlA,RefTable,Exclude);
if plotson
plot(WStbl.Time(I),Tselect.WSpredlogit);
ylim([-10 10])
end
% ax89(2)=subplot(3,1,2);
% plot(WStbl.Time(I),Tselect.Total,'g');

% WStbl.ExcludeAR(I)=ExcludeAR;
% WStbl.WSpredlogit(I)=WSpredlogit;
performance = PredictiveValue(1*(Tselect.WakeNoAR(~ExcludeAR)>0.5),1*(WSpredlogit(~ExcludeAR)>0),Tselect.WakeNoAR(~ExcludeAR));

acc1(i,1) = performance.Acc_sem_chance_p(1);
[i acc1(i,1)]
end
%
thres=0
performance = PredictiveValue(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>0.5),1*(WStbl.WSpredlogit(~WStbl.ExcludeAR)>thres),WStbl.WakeNoAR(~WStbl.ExcludeAR))

[minacc,minacc1]=min(acc1)

end


%%
% [608 699 701 982 992 994 1026 1044 1099 1131 1136 1153 1176 1183 
% 967 970] %967,970 missing Pbetalogfilt1 but ok now

%%
for x=5 %1:6
ExpforAcc = x;
UniqueSubjList = unique(BreathDataFullTable.Subjn);
acc1 = UniqueSubjList*NaN;
for i=1:length(UniqueSubjList)
    if ~(UniqueSubjList(i)>(10000*ExpforAcc) && UniqueSubjList(i)<=(10000*(ExpforAcc+1)))
        continue
    end
    if UniqueSubjList(i)>51185
        continue
    end
    I=find(Tacc.subjn==UniqueSubjList(i));
    if 1
        acc1(i,1) = max(Tacc.acc(I));
    else
        [~,tempi] = max(Tacc.accPred(I));
        acc1(i,1) = Tacc.acc(I(tempi));
    end
end
sum(~isnan(acc1))
accmediandataset(x)=nanmedian(acc1);
accmeandataset(x)=nanmean(acc1);
end

