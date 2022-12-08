%% BOOTSTRAP VERSION OF AHI AND PSG PARAMETERS
% run starthere.m first
%%
global settings AMasterSpreadsheet n

study='MESA';

t_start = clock;
[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');

[num,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
% use for mros1
MasterWorksheet(2907:end,:)=[]; num(2907:end,:)=[];
% use for mros2
% MasterWorksheet(1:2906,:)=[];

lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,1),'UniformOutput',false)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,1:2);

NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0

MaxM = size(num,1);

try
    M = max(settings.Mrange);
catch
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

if exist('MrangeOverride')
    Mrange=MrangeOverride;
else
    Mrange=settings.Mrange;
end
% intialize variables to zeros/NaN incase computation failed or no
% analysis file present
success = zeros(numel(Mrange),1);

AnalyzedDirectory=settings.AnalyzedDirectory;
savename=settings.savename;
ConvertedDirectory=settings.ConvertedDirectory;
%%
for n=Mrange
   
    n
    
    AnalyzedDir = [AnalyzedDirectory savename '_' num2str(n) '.mat'];
    ConvDir=[ConvertedDirectory Filenames{n,1} '.mat'];
      
    if exist(AnalyzedDir)==2 && exist(ConvDir)==2
    try
        clear W2 Epochs EpochsT EpochsT2 RespT3pa Evts SigT Time
        
        W1=load(AnalyzedDir);
        
        settings.useCentralPositionDatabase=1;
        settings.Boot=50;
              
        
        W2=load(ConvDir);
        if isfield(W2,'DataEventHypnog_Mat')
            SigT=array2table(W2.DataEventHypnog_Mat);
            SigT.Properties.VariableNames = W2.ChannelsList;
            W2 = rmfield(W2, {'DataEventHypnog_Mat','ChannelsList'});
            W2.SigT = SigT;
        else
            SigT=W2.SigT;
        end
        
        Time=SigT.Time;
        dt = Time(2,1) - Time(1,1);
        Fs = 1/dt;
        
        %% make epochs table
 
        switch study
            case 'MrOS'
                try
                    Epochs = Evts.Hypnogram;
                catch
                    Evts=EvtsData{1, 1};
                    Epochs = Evts.Hypnogram;
                    disp('EvtsData present');
                end
                
            case 'MESA' % use updated converted file. 2017 files are too old to use
             
                if ~settings.ignoreCPAPdata
                    [CPAPoff,CPAP] = getCPAP(SigT,settings,ChannelsList);
                else %a priori knowledge that CPAP is not administered
                    CPAPoff=ones(height(SigT),1);
                    CPAP=zeros(height(SigT),1);
                end
                % get Position code
                try
                    settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
                    settings.supinepositioncode = settings.positioncodes(1);
                catch
                    % above code failing if settings.protocol is char array.
                    settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(n,:));
                    settings.supinepositioncode = settings.positioncodes(1);
                end
                
                Evts=W2.Evts;
                Epochs =Evts.Hypnogram;
                Evts.RespT = EventRespSigtoRespT(SigT.EventsResp,SigT.Time);
%                 [Evts]=EventRespTable(SigT,Evts,SigT.Properties.VariableNames);
                [Evts]=getAHIAll(SigT,CPAPoff,Evts,SigT.Properties.VariableNames,settings,'EventsResp','EventsAr');
        end
    
        
        EpochsT = table(Epochs);
        EpochsT.EpochStart = Evts.Hypnogram_t;
        
        zerocheck = std(diff(EpochsT.EpochStart))
        
        EpochsT.EpochNum=(1:size(EpochsT,1))';
        EpochsT=[EpochsT(:,3) EpochsT(:,1:2)];
        
        %% make events table | TO DO - get Events.RespT from Analyzed data - only events meeting 3pa
        RespT3pa = Evts.RespT;
        try
        RespT3pa=movevars(RespT3pa,'InclAHI3a','After',width(RespT3pa));
        RespT3pa=movevars(RespT3pa,'InclAHI4','After','InclAHI3a');
        RespT3pa=movevars(RespT3pa,'InclAHI3','After','InclAHI3a');
        end
        
        Ind0=RespT3pa.InclAHI3a==0;
        temp = RespT3pa{:,2:end-3};
        temp(Ind0,:)=NaN;
        RespT3pa{:,2:end-3} = temp;
        
        
        I = sum(RespT3pa.EventCodes==[2:6],2)>0;
        RespT3pa = RespT3pa(I,:);
        
        
        RespT3pa.EpochNumber = floor((RespT3pa.EventStart - EpochsT.EpochStart(1))/30 + 1);
        
        EpochsT.Nevents = histc(RespT3pa.EpochNumber,(1:height(EpochsT))');
        
        %% stages
        EpochsT.N1 = 1*(EpochsT.Epochs==2);
        EpochsT.N2 = 1*(EpochsT.Epochs==1);
        EpochsT.N3 = 1*(EpochsT.Epochs==0);
        EpochsT.REM = 1*(EpochsT.Epochs==3);
        EpochsWake=1*(EpochsT.Epochs==4);
        TST=size(EpochsT.Epochs,1)-sum(EpochsWake);   % TST as number of epochs
        % sum(EpochsT.N1)+sum(EpochsT.N2)+sum(EpochsT.N3)+sum(EpochsT.REM)
        
        %% Position : Flateral, Arousals
        
        EpochsT.NArousals = NaN(size(EpochsT.Epochs,1),1);
        ArNum= SigT.EventsAr;
        
        EpochsT.ODI3 = NaN(size(EpochsT.Epochs,1),1);
        EpochsT.ODI4 = NaN(size(EpochsT.Epochs,1),1);
        % Remove SpO2 Artifacts
        SpO2Raw=SpO2ArtifactReject(SigT.SpO2,1/Fs);
        
        
        EpochsT.PosEpoch=NaN(size(EpochsT.Epochs,1),1);
        EpochsT.Flateral = NaN(size(EpochsT.Epochs,1),1);
        EpochsT.FSupine = NaN(size(EpochsT.Epochs,1),1);
        PosRaw = SigT.Position;
        
        if ~(isfield(settings,'useCentralPositionDatabase') && settings.useCentralPositionDatabase==1)
            [~,~,settings.poscodesdatabase] = xlsread(AMasterSpreadsheet,3,'B2:K55');
        else
            PosDatabaseDir = [settings.codedir '\Position\'];
            PosDatabaseSpreadsheet = [PosDatabaseDir, 'PositionDatabase.xlsx']; %
            [~,~,settings.poscodesdatabase] = xlsread(PosDatabaseSpreadsheet,1,'B2:K55');
        end
        settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
        [~,~,settings.protocol] = xlsread(AMasterSpreadsheet,1,'AF4:AF10003');
        settings.protocol(2907:end)= [];
        
        
        
        positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
        
        for xx=1:size(EpochsT.Epochs,1)
            Idx1=find(Time==EpochsT.EpochStart(xx));
            Idx2=find(Time==(EpochsT.EpochStart(xx)+30))-1;
            temp = PosRaw(Idx1:Idx2);
            temp2 = PositionTranslator(positioncodes,settings.positioncodesout,temp);
            EpochsT.PosEpoch(xx) = mode(temp2);
            EpochsT.FSupine(xx) = nanmean(PositionSelector(temp2,'Supine'));
            EpochsT.Flateral(xx)=1-EpochsT.FSupine(xx);
            
            % Arousal
            clear temp I1
            temp=ArNum(Idx1:Idx2);
            I = diff([NaN;temp]);
            I1 = find(I==1);
            %     I2 = find(I==-1);
            %     [I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
            %     lengthsN = I2N-I1N;
            EpochsT.NArousals(xx)=size(I1,1);
            
            % ODIs -- remove
            
            %     SpO2temp=SpO2Raw(Idx1:Idx2);
            %     Include=1-(SpO2temp==0); % include only when sao2 is non-zero
            %     [ODI3temp,ODI4temp]=CalcODIEpochVer(SpO2temp,dt,Include);
            %     EpochsT.ODI3(xx) = ODI3temp;
            %     EpochsT.ODI4(xx)=ODI4temp;
            
        end
        
        %% COMB 30
        
        EpochsT2=EpochsT;
        
        combW = 30;
        I = 1:size(EpochsT2,1);
        for j=1:3
            swap=j-1;
            clear I2 EpochT2Iteration criteriacomb temp2
            
            
            %     combW = 30;
            %     epochnumber = [1:400]';
            Etime = (EpochsT.EpochStart - EpochsT.EpochStart(1))/60;
            %     swap = 1;
            temp2 = mod(Etime + combW*swap,2*combW)/(2*combW);
            criteriacomb = temp2<0.5;
            if j==3
                criteriacomb = logical(ones(size(EpochsT2,1),1));
            end
            %
            %         I2 = datasample(I,length(I))';
            EpochT2Iteration=EpochsT2(criteriacomb,:);
            
            EpochT2Iteration(EpochT2Iteration.Epochs>3|isnan(EpochT2Iteration.Epochs),:)=[];
            %stages
            TST(j)=size(EpochT2Iteration.Epochs,1); %TST is in number of epochs
            FN1temp(j,1)=(sum(EpochT2Iteration.N1==1))/TST(j);
            FN2temp(j,1)=(sum(EpochT2Iteration.N2==1))/TST(j);
            FN3temp(j,1)=(sum(EpochT2Iteration.N3==1))/TST(j);
            FRemtemp(j,1)=(sum(EpochT2Iteration.REM==1))/TST(j);
            %position
            FLateraltemp(j,1)=nanmean(EpochT2Iteration.Flateral==1);
            %AHI
            TotHr=TST(j)*30/3600;
            AHItemp(j,1)=(sum(EpochT2Iteration.Nevents))/TotHr;
            %ArI
            ArItemp(j,1)=(sum(EpochT2Iteration.NArousals))/TotHr;
            %TST Total
            TSTtemp(j,1)=TST(j)*30/60; % in minutes
            %ODI
            Timetemp(j,1)=nanmean(EpochT2Iteration.EpochStart);
        end
        
        %%
        subj=n*(ones(3,1));
        J = [1 2 3]';
        Type = {'Odd','Even','Total'}';
        CombTtemp=table(subj,J,Type,FN1temp,FN2temp,FN3temp,FRemtemp,FLateraltemp,AHItemp,ArItemp,TSTtemp,Timetemp);
        CombTtemp.Properties.VariableNames={'Subject','TypeID','Type','FN1','FN2','FN3','FREM','FLateral','AHI',...
            'ArI','TST','MeanTime'};
        
        if ~isempty(CombTtemp)
            if ~exist('CombTPSGParam')
                CombTPSGParam = CombTtemp;
            else
                CombTPSGParam = [CombTPSGParam;CombTtemp];
            end
        end
        clear CombTtemp FN1temp FN2temp FN3temp FRemtemp FLateraltemp AHItemp ArItemp;
        % end subject analysis
        
        %% BootStrap
        booton=1;
        if booton
            % remove wake epochs
            EpochsT2=EpochsT;
            EpochsT2(logical(EpochsWake),:)=[];
            
            I = 1:size(EpochsT2,1);
            for j=1:settings.Boot
                clear I2 EpochT2Iteration
                I2 = datasample(I,length(I))';
                EpochT2Iteration=EpochsT2(I2,:);
                %stages
                TST(j)=size(EpochT2Iteration.Epochs,1);
                FN1temp(j,1)=(sum(EpochT2Iteration.N1==1))/TST(j);
                FN2temp(j,1)=(sum(EpochT2Iteration.N2==1))/TST(j);
                FN3temp(j,1)=(sum(EpochT2Iteration.N3==1))/TST(j);
                FRemtemp(j,1)=(sum(EpochT2Iteration.REM==1))/TST(j);
                %position
                FLateraltemp(j,1)=nanmean(EpochT2Iteration.Flateral==1);
                %AHI
                TotHr=TST(j)*30/3600;
                AHItemp(j,1)=(sum(EpochT2Iteration.Nevents))/TotHr;
                %ArI
                ArItemp(j,1)=(sum(EpochT2Iteration.NArousals))/TotHr;
                %TSTtemp
                TSTtemp(j,1)=TST(j);
                %ODI
                
                
            end
            
            %%
            clear BootTtemp
            subj=n*(ones(settings.Boot,1));
            BootTtemp=table(subj,FN1temp,FN2temp,FN3temp,FRemtemp,FLateraltemp,AHItemp,ArItemp,TSTtemp);
            BootTtemp.Properties.VariableNames={'Subject','FN1','FN2','FN3','FREM','FLateral','AHI',...
                'ArI','TSTNEpochs'};
            
            if ~isempty(BootTtemp)
                if ~exist('BootTPSGParam')
                    BootTPSGParam = BootTtemp;
                else
                    BootTPSGParam = [BootTPSGParam;BootTtemp];
                end
            end
            clear BootTtemp FN1temp FN2temp FN3temp FRemtemp FLateraltemp AHItemp ArItemp TSTtemp;
            % end subject analysis
        end
    catch
    end
    end
end

fname=[settings.workdir,'Summary' filesep 'AhiTstBootT_Oct2021_2.mat'];
save(fname,'CombTPSGParam','BootTPSGParam','settings');
