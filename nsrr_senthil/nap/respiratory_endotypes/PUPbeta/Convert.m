function Convert(n1,settings1,AMasterSpreadsheet1,HarmonizedChannelSpreadsheet1)
%% global variables and settings
global ChannelsList ChannelsFs Evts settings AMasterSpreadsheet n Info HarmonizedChannelSpreadsheet

% PUPbetaGUI is a global, figure handle for the PUPbeta window
% AmasterSpreadheet is coming in from PUPbeta gui, keep as is.
% handletext is coming in from PUPbeta gui, keep as is.

if exist('settings1')
    settings=struct();
    settings=settings1;
end
if exist('AMasterSpreadsheet1')
    AMasterSpreadsheet=AMasterSpreadsheet1;
end

if exist('HarmonizedChannelSpreadsheet1')
    HarmonizedChannelSpreadsheet=HarmonizedChannelSpreadsheet1;
end

n=n1;


%%
disp(' '); % add row space for visual clarity in command window
displaytext=['Starting Convert to .mat for: ' settings.Filenames{n}(1:end-4)];
disp(displaytext); 
% set(handletext,'String',displaytext); drawnow;
settings.PosFigName = [settings.ConvertedDirectory, settings.Filenames{n}(1:end-4) 'pos'];
[SigT,ChannelsFs,Evts,Info] = GeneralImport(n);

%%
try
    if isfield(settings,'SaveNameIsSubjID_NotEDF') && settings.SaveNameIsSubjID_NotEDF==1 % use Subj ID, not EDF name as output filename (DLM, was for PACE data)
        Filenames_ = settings.MasterWorksheet(:,11); % reassigning this every time isn't efficient
        fname=[Filenames_{n}(1:end-4) '_XHz' '.mat'];
    else
        fname=[settings.Filenames{n}(1:end-4) '_XHz' '.mat'];
    end
    displaytext=['Saving data to: ' fname];
    disp(displaytext); 
%     set(handletext,'String',displaytext); drawnow;
    save([settings.ConvertedDirectory, fname],'SigT','ChannelsFs','Evts','Info','-v7.3');
    displaytext='Finished saving data';
    disp(displaytext);
%     set(handletext,'String',displaytext); drawnow;
catch me
    disp(me.message);
    % disp(getReport(me));
end


%% GeneralImport: Currently works for Minerva, Profusion, Spike
function [SigT,ChannelsFs,Evts,Info] = GeneralImport(n)
global settings HarmonizedChannelSpreadsheet AMasterSpreadsheet
% handletext --removed display to GUI for parallel processing
Filenames=settings.Filenames;
ChannelNumbers=settings.ChannelNumbers;
system = char(Filenames(n,7));

try
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    switch system
        case {'ProfusionXML','ProfusionXML_APN','Minerva','ProfusionXML_MrOs','BrainRT',...
                'RemLogic','RemLogicXML','RemLogicTextFr','RemLogicText',...
                'Alice','AliceTaiwan','AliceG3','AliceG3Machine','AliceG3Taiwan','Deltamed','SomnoStarEpochsXls','GrassTwinTxt',...
                'SandmanPlusLabview','Sandman','NKxdf','NKxdfB','EdfUnscored','Danalyzer','Michele','ApneaLink','AnnotEannot',...
                'NSRR','MignotAutostage','Nox','ZMachineSynergy'}
            SignalFormat='EDF';
            
            try % causing error when ChannelNumbers are 1 row while using harmonized list, so moved to try-catch--RMA
            if isnan(ChannelNumbers(n,1))
                disp('Flow channel number is NaN, skipping')
                Info=NaN;
                ColumnHeads=NaN;
                DataEventHypnog_Mat=NaN;
                ChannelsList=NaN;
                ColumnHeadsList=NaN;
                return
            end
            end
            
        case 'Spike'
            SignalFormat='Spike';
        case 'SpikeDise'
            SignalFormat='SpikeDise';
        case 'NoxMat'
            SignalFormat='NoxMat';
        case 'NoxMatSAS'
            SignalFormat='NoxMatSAS';
        case 'NoxMatT3'
            SignalFormat='NoxMatT3';
        
%             % DLM commented "case 'Nox'" because it already appears in list above
%         case 'Nox' % EDF/XLS format from Nox
%             SignalFormat='Nox';
            
    end
    
    %% Signals
    %This section imports respiratory flow and other polysomnography signal data from the EDF files
    displaytext='Loading signals';
    disp(displaytext); 
%     set(handletext,'String',displaytext); drawnow;
    
    %Default Signals List
    ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12','Position','Pmask','EKG','EKG2','LOC','ROC'};
    %SS: first channel renamed to Flow from Pnasal, otherwise later the function will not function
    
    disp(['SignalFormat: ' SignalFormat]);
    
    if strcmp(SignalFormat,'EDF')
        %% EDF
        ChannelsFs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        EDFfilenamedir = [directory fname];
        
        %% Use Harmonized Channel Labels
        if isfield(settings,'UseHarmonizedChannelNumbers')&& settings.UseHarmonizedChannelNumbers==1
            ChannelNumbers=NaN(length(Filenames),24);
            if isfield (settings,'NSRRSave') && settings.NSRRSave==1 % removing pupstart folder in nsrr
                 HChannelSpreadsheet=[settings.HarmonizedChannelSpreadsheetDir 'HarmonizedChannelLabels.xlsx'];
            else             

            HChannelSpreadsheet=[settings.workdir 'PUPStart' filesep 'HarmonizedChannelLabels.xlsx'];
            end
            ChannelNumbers(n,:) = getChannelNumbersHarmonizedLabels(EDFfilenamedir,HChannelSpreadsheet);
            Info.ChannelsListSource=ChannelNumbers;
             if isnan(ChannelNumbers(n,1))
                disp('Flow channel number is NaN, skipping')
                Info=NaN;
                ColumnHeads=NaN;
                DataEventHypnog_Mat=NaN;
                ChannelsList=NaN;
                ColumnHeadsList=NaN;
                return
            end
        end
        
        
        
        % Look to see if there is HarmonizedChannelsFile (if have this,
        % then skip getEDFchannelnames earlier)
        % Read in EDFLabels from EDF (use getEDFchannelNAmes approach)
        % Make channelnamesoptions structure variable by searching individual signals from HarmonizedChannelsFile for actual channel in EDFlabels (as per Spike below)
        % If found, record channel number for output in ChannelNumbers array
        % Keep "Info.ChannelsListSource" - should match length of ChannelsList for source data channels. OK not to update as we add calculated channels.
        
        [Header,signalHeader] = blockEdfLoad(EDFfilenamedir); signalHeader=struct2table(signalHeader); %used if blockEDFLoad methods are used
        
        if isfield(settings,'BlockEDFload') && settings.BlockEDFload==2 %fails if channels are labelled the same; also fails sometimes for unknown reasons
            tempCHN=ChannelNumbers(n,:);
            tempCHN(tempCHN>height(signalHeader))=NaN; %error check
            I = ~isnan(tempCHN);
            siglabs = signalHeader.signal_labels(tempCHN(I));
            disp('Loading EDF in Block Mode');
            [~,~,signalCellX] = blockEdfLoad(EDFfilenamedir,siglabs');
            I2 = cumsum(I); I2(isnan(tempCHN))=NaN;
        end
        
        for i=1:length(ChannelsList)
            displaytext=['Collecting channel:' ChannelsList{i}];
            try
                disp(displaytext);
%                 set(handletext,'String',displaytext); drawnow;
                if isfield(settings,'BlockEDFload') && settings.BlockEDFload==1
                    [~,~,signalCell] = blockEdfLoad(EDFfilenamedir,signalHeader.signal_labels(ChannelNumbers(n,i)));
                    eval([ChannelsList{i} ' = signalCell{1};']);
                    ChannelsFs(i)=signalHeader.samples_in_record(ChannelNumbers(n,i));
                    LabelTemp=signalHeader.signal_labels{ChannelNumbers(n,i)};
                    Filt=signalHeader.prefiltering{ChannelNumbers(n,i)};
                elseif isfield(settings,'BlockEDFload') && settings.BlockEDFload==2
                    eval([ChannelsList{i} ' = signalCellX{I2(i)};']);
                    ChannelsFs(i)=signalHeader.samples_in_record(ChannelNumbers(n,i));
                    LabelTemp=signalHeader.signal_labels{ChannelNumbers(n,i)};
                    Filt=signalHeader.prefiltering{ChannelNumbers(n,i)};
                else
                    eval(['[' ChannelsList{i} ',ChannelsFs(i),~,~,LabelTemp,~,Filt,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
                end
                %e.g. [x,Fs,Start_date,Start_time,Label,Dimension,Filt,Coef,Nmb_chans,N] = readedfrev3(fname,ch,t1,t2)
                displaytext=['   Found channel labelled:' LabelTemp ' at ' num2str(ChannelsFs(i)), ' Hz'];
                disp(displaytext); 
%                 set(handletext,'String',displaytext); drawnow;
                ChannelFound(i)=1;
                if isfield(settings,'reportEDFprefilters') && settings.reportEDFprefilters==1
                    % list prefilters for channel
                    outstr = [LabelTemp, ' - EDF PreFilter:', Filt]; disp(outstr);
                end
            catch me
                displaytext=['   No channel:' ChannelsList{i}];
                if 0
                    displaytext=['   No channel:' ChannelsList{i} ' ' me.message];
                end
                disp(displaytext); 
%                 set(handletext,'String',displaytext); drawnow;
                %eval([ChannelsList{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
            end
        end
        
    elseif strcmp(SignalFormat,'Spike')
        %% Spike
        %Overwrite default ChannelsList
        ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','LOC','ROC',...
            'Position','Pmask','EKG','alphaFlow','kFlow','Pes','Edi','Pepi','GGpmax','Pnasal','SnoreDB','NoxAudio',...
            'Chin','StimAmpl'};
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        %         MultiplyPositionChannelby=5; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat','Sao2ear','Sao2fing','SaO2fing'};
        channelnameoptions.Position={'Position','position','Pos','pos','NoxPos', 'Body'};
        channelnameoptions.pCO2={'PCO2','CO2_Ana','CO2_anal', 'CO2','pCO2','ETCO2'};
        channelnameoptions.pO2={'pO2','O2_Ana','O2_anal','PO2','ETO2'};
        channelnameoptions.EEG1={'EEG1','EEG_C3_A2_clean','EEG_C3_A2','C3_A2','C3-A2','EEG_C3_A','C3','C3_M2'};
        channelnameoptions.EEG2={'EEG2','EEG_C4_A1_clean','EEG_C4_A1','C4_A1','C4-A1','C4','C4_M1'};
        channelnameoptions.EEG3={'EEG3','EEG_F3_A2_clean','EEG_F3_A2','F3_A2','F3-A2','F3','F3_M2', 'EEG_A1_A2'};
        channelnameoptions.EEG4={'EEG4','EEG_O2_A1_clean','EEG_O2_A1','O2_A1','O2-A1','EEG_O2_A','O2','O2_M1'};
        channelnameoptions.EEG5={'EEG5','EEG_F4_A1_clean','EEG_F4_A1','F4_A1','F4-A1','F4','F4_M1'};
        channelnameoptions.EEG6={'EEG6','O1','O1_M2'};
        channelnameoptions.LOC={'LOC','EOG_L','LEOG','E1_M2', 'EOG_LOC_A'};
        channelnameoptions.ROC={'ROC','EOG_R','REOG','E2_M2', 'EOG_ROC_A'};
        channelnameoptions.EKG={'EKG','ECG','ECG1_ECG2','ECG_a_ECG_b', 'ECG_I'};
        channelnameoptions.Thorax={'ThNoxRIP','Thorax','RC','Chest','CHEST','Belt2','Thor', 'Effort_TH'};
        channelnameoptions.Abdomen={'AbNoxRIP','Abdomen','ABD','Abdom','ABDM','Belt1','Abdo', 'Effort_AB'};
        channelnameoptions.alphaFlow={'alphaFlow'};
        channelnameoptions.kFlow={'kFlow'};
        channelnameoptions.Pes={'Pes','pes','Pes_clean'};
        channelnameoptions.Edi={'Edi','edi'}; %'EMGdi'
        channelnameoptions.Pepi={'Epi','epi','Pepi','pepi'};
        %channelnameoptions.FlowEdi={'FlowEdi'};
        %channelnameoptions.FlowPes={'FlowPes'};
        channelnameoptions.GGpmax={'GGpmax','GGPmax'};
        channelnameoptions.Flow={'Vflow','Flow','flow','FLOW','NoxPnasal','Pnasal','PNasal','PNASAL','Pmask','PMask','Pressure','Press_Pat'}; %'Vflow'
        channelnameoptions.Pnasal={'Pnasal','PNasal','Pressure'}; %secondary Pnasal signal if there is also Flow
        channelnameoptions.Pulse={'Pulse', 'PulseRate'}; 
        channelnameoptions.Pleth={'Pleth'}; 
        channelnameoptions.Pmask={'CPAP','Pmask','PMask'};
        channelnameoptions.Snore={'cSnore'};
        channelnameoptions.SnoreDB={'SnoreDB'};
        channelnameoptions.NoxAudio={'NoxAudio','Audio_Volume'};
        channelnameoptions.Chin={'CHIN','Chin','ChinEMG'};
%         channelnameoptions.StimAmpl={'StimAmpl'};

        
        channelnamestemp=fieldnames(channelnameoptions);
        for i=1:length(channelnamestemp)
            temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp) %go through each possible channel label
                %Does it exist?
                for j=1:length(w)
                    %w(j).name
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';'])
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        ChannelsFs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        ChannelsStart = NaN*ones(length(ChannelsList),1);
        
        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    ChannelsFs(i)=1./eval([ChannelsList{i} '.interval']);
                    ChannelsStart(i)=eval([ChannelsList{i} '.start']);
                    if ChannelsStart(i) < 0.1 % then it likely was not imported between cursor regions, so call it 0 for clarity
                        ChannelsStart(i) = 0;
                    end
                    % this IF was added by DLM. For some reason, alphaFlow
                    % channel in 1161 had an interval of 0.0013, however
                    % the data length was the same as other channels,
                    % implying that the interval should have been the same
                    % as the other channels (i.e. 0.008).
                    if 0 % this is not required with new spike to matlab
                        % exports, and should be removed now.
                        if (ChannelsFs(i) > 740) && (ChannelsFs(i) < 741)
                            str = ['---- CAUTION ----   Manual set Fs for ', ChannelsList{i}]; disp(str);
                            ChannelsFs(i) = 125;
                            pause;
                        end
                    end
                catch me
                    ChannelsFs(i)=settings.Fs; %assumed
                end
                temp = eval([ChannelsList{i} '.values']);
                try
                    ChannelTitlesoriginal{i}=eval([ChannelsList{i} '.title']);
                catch me
                    ChannelTitlesoriginal{i}=[];
                end
                eval(['clear ' ChannelsList{i}]);
                eval([ChannelsList{i} '=temp;']);
            else
                % display to command window, but not to PUPbeta gui
                disp(['strewth, no ', ChannelsList{i}, ' found']);
            end
        end
        
        %% add code to import external Edi file here
        try
            Edifilename = [Filenames{n,4} Filenames{n,1}(1:end-4) '_Edi.mat'];
            if exist(Edifilename)==2
                w2 = matfile(Edifilename);
                List = {'EMGdi1','EMGdi2','EMGdi3','EMGdi4','EMGdi5','EMGdiBest','EMGdiCh','EMGdiX'};
                ListRename = {'Edi1','Edi2','Edi3','Edi4','Edi5','EdiBest','EdiCh','Edi'};
                %             load(Edifilename,List{:});
                
                for i=1:length(List)
                    try
                        eval([ListRename{i} ' = w2.' List{i} '*1000'';']);
                        'successfully loaded external Edi signals'
                    catch me
                        'failed to load Edi signal'
                    end
                end
                
                %             Edi1 = w.EMGdi1;
                %             Edi2 = w.EMGdi2;
                %             Edi3 = w.EMGdi3;
                %             Edi4 = w.EMGdi4;
                %             Edi5 = w.EMGdi5;
                %             EdiBest = w.EMGdiBest;
                %             EdiCh = w.EMGdiCh;
                %             Edi = w.EMGdiX;
                
                EdiChannelFound=zeros(length(ListRename),1);
                for i=1:length(ListRename)
                    EdiChannelFound(i)=(exist(ListRename{i})==1)*1;
                end
                List = List(EdiChannelFound==1);
                ListRename = ListRename(EdiChannelFound==1);
                
                %eventually export this also into ConvertedInfo
                try
                    EdiInfo = mergestructs(w2.BestEdiInfo,w2.settings);
                catch me
                    'fix EdiInfo'
                end
                ChannelsList = [ChannelsList ListRename];
                i = find(strcmp(ChannelsList,'Flow'));
                ChannelsFs = [ChannelsFs;repmat(ChannelsFs(i),length(ListRename),1)];
                ChannelFound = [ChannelFound;EdiChannelFound];
            end
        catch
            disp('failed to import external Edi');
        end
   
    elseif strcmp(SignalFormat,'SpikeDise')
        %% Spike files with DISE data - different labels
        
        ChannelsList={'Flow','FlowNasal','FlowOral','Snore','SnoreDB','SpO2'...
            };
        
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat','Sao2ear','Sao2fing','SaO2fing'};
        channelnameoptions.Flow={'Vflow','Flow','FlowCombined','flow','Pnasal','PNasal','Pmask','PMask'}; %'Vflow'
        channelnameoptions.FlowNasal = {'FlowNasal'};
        channelnameoptions.FlowOral = {'FlowOral'};
        channelnameoptions.Snore={'cSnore'};
        channelnameoptions.SnoreDB={'SnoreDB'};
        
        channelnamestemp=fieldnames(channelnameoptions);
        for i=1:length(channelnamestemp)
            temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';'])
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        ChannelsFs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        ChannelsStart = NaN*ones(length(ChannelsList),1);
        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    ChannelsFs(i)=1./eval([ChannelsList{i} '.interval']);
                    ChannelsStart(i)=eval([ChannelsList{i} '.start']);
                    if ChannelsStart(i) < 0.1 % then it likely was not imported between cursor regions, so call it 0 for clarity
                        ChannelsStart(i) = 0;
                    end
                    % this IF was added by DLM. For some reason, alphaFlow
                    % channel in 1161 had an interval of 0.0013, however
                    % the data length was the same as other channels,
                    % implying that the interval should have been the same
                    % as the other channels (i.e. 0.008).
                    if 0 % this is not required with new spike to matlab
                        % exports, and should be removed now.
                        if (ChannelsFs(i) > 740) && (ChannelsFs(i) < 741)
                            str = ['---- CAUTION ----   Manual set Fs for ', ChannelsList{i}]; disp(str);
                            ChannelsFs(i) = 125;
                            pause;
                        end
                    end
                catch me
                    ChannelsFs(i)=settings.Fs; %assumed
                end
                temp = eval([ChannelsList{i} '.values']);
                try
                    ChannelTitlesoriginal{i}=eval([ChannelsList{i} '.title']);
                catch me
                    ChannelTitlesoriginal{i}=[];
                end
                eval(['clear ' ChannelsList{i}]);
                eval([ChannelsList{i} '=temp;']);
            else
                % display to command window, but not to PUPbeta gui
                disp(['strewth, no ', ChannelsList{i}, ' found']);
            end
        end
        
    elseif strcmp(SignalFormat,'NoxMat')
        %% Nox
        %Overwrite default ChannelsList
        ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
            'X_Axis','Y_Axis','Z_Axis','Pmask','EKG','Pes','TcPCO2','AudioVolume','Pleth','Pulse'...
            };
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        MultiplyPositionChannelby=1; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SpO2'};
        channelnameoptions.X_Axis={'X_Axis'};
        channelnameoptions.Y_Axis={'Y_Axis'};
        channelnameoptions.Z_Axis={'Z_Axis'};
        
        channelnameoptions.EEG1={'C3'};
        channelnameoptions.EEG2={'M2'};
        channelnameoptions.EEG3={'C4'};
        channelnameoptions.EEG4={'M1'};
        channelnameoptions.EEG5={'F3'};
        channelnameoptions.EEG6={'M2'};
        channelnameoptions.EEG7={'O1'};
        channelnameoptions.EEG8={'M2'};
        channelnameoptions.EEG9={'F4'};
        channelnameoptions.EEG10={'M1'};
        channelnameoptions.EEG11={'O2'};
        channelnameoptions.EEG12={'M1'};
        
        channelnameoptions.LOC={'E1'};
        channelnameoptions.ROC={'E2'};
        
        channelnameoptions.EKG={'EKG','ECG','ecg'};
        channelnameoptions.Thorax={'Inductance_Thorax'};
        channelnameoptions.Abdomen={'Inductance_Abdomen'};
        
        channelnameoptions.Pes={'Millar_Pes','channel_3','Channel_3','Channel_5'};
        
        channelnameoptions.Flow={'Channel_7','Pneumotach_Flow','pneumoflow','Pneumoflow','Mask_Pressure','Pression_masque'}; %'Vflow'
        
        channelnameoptions.Pmask={'Mask_Pressure','Pression_masque'};
        
        channelnameoptions.TcPCO2={'TcCO2'};
        channelnameoptions.NoxAudio={'NoxAudio','Audio_Volume'};
        channelnameoptions.Pleth={'Pleth'};
        channelnameoptions.Pulse={'Pulse','Pouls'};
        channelnamestemp=fieldnames(channelnameoptions);
        
        for i=1:length(channelnamestemp)
            temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';'])
                        ChannelTitlesTemp{i} = char(temp(nn));
                        %eval([channelnamestemp{i} '.title=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        ChannelsFs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        %ChannelStartTime = NaN*ones(length(ChannelsList),1);
        Signals = filehandle.Signals;
        
        TimeInfo = filehandle.TimeInfo;
        
        %%
        %get sampling rates
        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    %I=find(strcmp(channelnamestemp,ChannelsList{i}));
                    %I2=find(strcmp(Signals.labelsNoSpace,ChannelTitlesTemp{I}));
                    %ChannelsFs(i)=eval([ChannelsList{i} '.fs']);
                    ChannelsFs(i)=eval([ChannelsList{i} '.fsTrue']);
                catch me
                    ChannelsFs(i)=settings.Fs; %assumed
                    disp(['missing sampling rate assuming settings.Fs']);
                end
            else
                % display to command window, but not to PUPbeta gui
                disp(['strewth, no ', ChannelsList{i}, ' found']);
            end
        end
        
        for i=1:length(ChannelsList)
            if ChannelFound(i)==1
                temp = eval([ChannelsList{i} '.data']);
                eval([ChannelsList{i} '=temp(:);']);
            end
        end
        
        %NoxPositionSignals
        Pos_AccMagnitude = sqrt(X_Axis.^2+Y_Axis.^2+Z_Axis.^2);
        %supine: Y=-1,X=0,Z=0; code:1/(pi/2)*atan2(0,1)
        %upright: Y=-0,X=0,Z=-1; code:1/(pi/2)*atan2(0,1)
        Position = 1/(pi/2)*atan2(X_Axis./Pos_AccMagnitude,-Y_Axis./Pos_AccMagnitude); %quarter turns
        %Incline = 1/(pi/2)*atan2(-Z_Axis./Pos_AccMagnitude,sqrt(Y_Axis.^2+X_Axis.^2)./Pos_AccMagnitude); %quarter turns
        Incline = -Z_Axis./Pos_AccMagnitude;
        % Supine, position = 0 incline = 0 for x_axis = 0  y_axis = -1  z_axis = 0
        % lateral (left?), postion = -1 incline = 0 for x_axis = -1  y_axis = 0  z_axis = 0
        % lateral (right?), postion = 1 incline = 0 for x_axis = 1  y_axis = 0  z_axis = 0

        I = find(strcmp(ChannelsList,'X_Axis')==1)
        %Fs_temp = ChannelsFs(I);
        ChannelFound(I+2)=0; %assumes X,Y,Zaxis are together
        ChannelsList{I}='Position';
        ChannelsList{I+1}='Incline';
        
        %calibrate Millar Pes:
        if exist('Pes')==1
            Pes = Pes/0.0073556; %([0.25/25 V/mmHg] / 1.35951 mmHg/cmH2O = 0.007356 V/cmH2O); 0.25/25/1.35951 = 0.007356
        end
        
        
    elseif strcmp(SignalFormat,'NoxMatSAS')
        %% Nox
        %Overwrite default ChannelsList
        ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','LOC','ROC','REF1','REF2',...
            'X_Axis','Y_Axis','Z_Axis','EKG','AudioVolume','Pleth','Pulse','SpO2_B_B','AFZ','Thorax_Fast','Abdomen_Fast'...
            'ECG_LA','ECG_LF','ECG_RA','SnoreDB','Snore'...
            };
        %Note EEG1-4 are referenced against E3,4l EEG5-8 are repeated signals that are unreferenced
        %Consider adding SpO2_B_B Thorax_Fast Abdomen_Fast AFZ, other EKGs, (ECG_LA ECG_LF ECG_RA), key impedances
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        MultiplyPositionChannelby=1; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SpO2'};
        channelnameoptions.SpO2_B_B={'SpO2_B_B'};
        channelnameoptions.X_Axis={'X_Axis'};
        channelnameoptions.Y_Axis={'Y_Axis'};
        channelnameoptions.Z_Axis={'Z_Axis'};
        %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal'};
        %channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
        channelnameoptions.AFZ={'AFZ'};
        channelnameoptions.EEG1={'AF3'}; %needs to be AF3 - (E3+E4)/2;
        channelnameoptions.EEG2={'AF4'};
        channelnameoptions.EEG3={'AF7'};
        channelnameoptions.EEG4={'AF8'};
        channelnameoptions.EEG5={'AF3'}; %unreferenced (actually referenced to AF7,AF8 average)
        channelnameoptions.EEG6={'AF4'};
        channelnameoptions.EEG7={'AF7'};
        channelnameoptions.EEG8={'AF8'};
        
        if isfield(settings, 'mouthEEG') && settings.mouthEEG == 1
            channelnameoptions.EEG9={'ECG_LF'};
            channelnameoptions.EEG10={'ECG_RA'};
        end
        
        channelnameoptions.LOC={'E1'}; %needs to be E1-E4
        channelnameoptions.ROC={'E2'}; %needs to be E2-E3
        channelnameoptions.REF1={'E3'};
        channelnameoptions.REF2={'E4'};
        channelnameoptions.EKG={'EKG','ECG','ecg'};
        channelnameoptions.EKG_LA={'ECG_LA'};
        channelnameoptions.EKG_LF={'ECG_LF'};
        channelnameoptions.EKG_RA={'ECG_RA'};
        channelnameoptions.Thorax={'Inductance_Thorax'};
        channelnameoptions.Thorax_Fast={'Thorax_Fast'};
        channelnameoptions.Abdomen={'Inductance_Abdomen'};
        channelnameoptions.Abdomen_Fast={'Abdomen_Fast'};
        channelnameoptions.Flow={'Mask_Pressure','Pression_masque','Nasal_Pressure'}; %'Vflow'
        %channelnameoptions.AudioVolume={'Audio_Volume'};
        channelnameoptions.NoxAudio={'NoxAudio','Audio_Volume'};
        channelnameoptions.Pleth={'Pleth','Pulse_Wave__Pleth_'};
        channelnameoptions.Pulse={'Pulse','Pouls'};
        channelnameoptions.Snore={'Snore'};
        channelnameoptions.SnoreDB={'SnoreDBPhone'}; %For smart phone recorded sound
        
        channelnamestemp=fieldnames(channelnameoptions);
        
        for i=1:length(channelnamestemp)
            temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';'])
                        ChannelTitlesTemp{i} = char(temp(nn));
                        %eval([channelnamestemp{i} '.title=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        ChannelsFs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        %ChannelStartTime = NaN*ones(length(ChannelsList),1);
        Signals = filehandle.Signals;
        
        TimeInfo = filehandle.TimeInfo;
        
        %% code to import external audio file
        try
            Audiofilename = [Filenames{n,4} Filenames{n,1}(1:end-4) 'Snore.mat'];
            if exist(Audiofilename)==2
                w3 = matfile(Audiofilename);
                List = {'Snore'};
                
                for i=1:length(List)
                    try
                        eval([List{i} ' = w3.' List{i} ';']);
                        disp('successfully loaded external Audio signals')
                    catch me
                        disp('failed to load Audio signal')
                    end
                end
                
                % TEMPORARY
                if isfield(Snore,'fs')
                    Snore.data = Snore.data';
                    Snore.fsTrue = Snore.fs;
                end
                
                
                %                 SnoreChannelFound=zeros(length(List),1);
                %                 for i=1:length(List)
                %                     SnoreChannelFound(i)=(exist(List{i})==1)*1;
                %                 end
                %                 List = List(SnoreChannelFound==1);
                %
                %                 %eventually export this also into ConvertedInfo
                %                 try
                %                     SnoreInfo = mergestructs(w3.settings);
                %                 catch me
                %                     'fix SnoreInfo'
                %                 end
                %                 ChannelsList = [ChannelsList List];
                %                 i = find(strcmp(ChannelsList,'Flow'));
                %                 ChannelsFs = [ChannelsFs;repmat(ChannelsFs(i),length(List),1)];
                %                 ChannelFound = [ChannelFound;SnoreChannelFound];
            end
        catch
            disp('failed to import external Audio');
        end
        
        
        %get sampling rates
        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    %I=find(strcmp(channelnamestemp,ChannelsList{i}));
                    %I2=find(strcmp(Signals.labelsNoSpace,ChannelTitlesTemp{I}));
                    %ChannelsFs(i)=eval([ChannelsList{i} '.fs']);
                    ChannelsFs(i)=eval([ChannelsList{i} '.fsTrue']);
                catch me
                    ChannelsFs(i)=settings.Fs; %assumed
                    disp(['missing sampling rate assuming settings.Fs']);
                end
                %                 temp = eval([ChannelsList{i} '.data']);
                %                 try
                %                     ChannelTitlesoriginal{i}=ChannelTitlesTemp{i};
                %                 catch me
                %                     ChannelTitlesoriginal{i}=[];
                %                 end
                %                 if 1
                %                 eval(['clear ' ChannelsList{i}]);
                %                 eval([ChannelsList{i} '=temp;']);
                %                 end
            else
                % display to command window, but not to PUPbeta gui
                disp(['strewth, no ', ChannelsList{i}, ' found']);
            end
        end
        
        for i=1:length(ChannelsList)
            if ChannelFound(i)==1
                temp = eval([ChannelsList{i} '.data']);
                eval([ChannelsList{i} '=temp(:);']);
            end
        end
        
        %NoxPositionSignals
        PosLength = max([length(X_Axis) length(Y_Axis) length(Z_Axis)]);
        try
            if length(X_Axis)<PosLength
                nandiff = nan(PosLength-length(X_Axis),1);
                X_Axis = [X_Axis; nandiff];
            end
            if length(Y_Axis)<PosLength
                nandiff = nan(PosLength-length(Y_Axis),1);
                Y_Axis = [Y_Axis; nandiff];
            end
            if length(Z_Axis)<PosLength
                nandiff = nan(PosLength-length(Z_Axis),1);
                Z_Axis = [Z_Axis; nandiff];
            end
        catch me
        end
        
        Pos_AccMagnitude = sqrt(X_Axis.^2+Y_Axis.^2+Z_Axis.^2);
        %supine: Y=-1,X=0,Z=0; code:1/(pi/2)*atan2(0,1)
        %upright: Y=-0,X=0,Z=-1; code:1/(pi/2)*atan2(0,1)
        Position = 1/(pi/2)*atan2(X_Axis./Pos_AccMagnitude,-Y_Axis./Pos_AccMagnitude); %quarter turns
        %Incline = 1/(pi/2)*atan2(-Z_Axis./Pos_AccMagnitude,sqrt(Y_Axis.^2+X_Axis.^2)./Pos_AccMagnitude); %quarter turns
        Incline = -Z_Axis./Pos_AccMagnitude;
        
        I = find(strcmp(ChannelsList,'X_Axis')==1);
        %Fs_temp = ChannelsFs(I);
        ChannelFound(I+2)=0; %assumes X,Y,Zaxis are together
        ChannelsList{I}='Position';
        ChannelsList{I+1}='Incline';
        
        % Referencing Electrodes:
        disp('Referencing Electrodes [added 20200407, consider repeating prior conversions]');
        try
            tempref = (REF1+REF2)/2; %(E3+E4)/2
            EEG1 = EEG1 - tempref;
            EEG2 = EEG2 - tempref;
            EEG3 = EEG3 - tempref;
            EEG4 = EEG4 - tempref;
        catch me
            disp('failed EEG ref')
        end
        
        try
            LOC = LOC - REF2; %E1 - E4
            ROC = ROC - REF1; %E2 - E3
        catch me
            disp('failed EOG ref')
        end
        
          
    elseif strcmp(SignalFormat,'NoxMatT3')
        %% Nox
        %Overwrite default ChannelsList
        ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','LOC','ROC','REF1','REF2',...
            'X_Axis','Y_Axis','Z_Axis','EKG','AudioVolume','Pleth','Pulse','SpO2_B_B','AFZ','Thorax_Fast','Abdomen_Fast'...
            'ECG_LA','ECG_LF','ECG_RA','Snore'...
            };
        %Note EEG1-4 are referenced against E3,4l EEG5-8 are repeated signals that are unreferenced
        %Consider adding SpO2_B_B Thorax_Fast Abdomen_Fast AFZ, other EKGs, (ECG_LA ECG_LF ECG_RA), key impedances
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        MultiplyPositionChannelby=1; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SpO2'};
        channelnameoptions.SpO2_B_B={'SpO2_B_B'};
        channelnameoptions.X_Axis={'X_Axis'};
        channelnameoptions.Y_Axis={'Y_Axis'};
        channelnameoptions.Z_Axis={'Z_Axis'};
        %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal'};
        %channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
        channelnameoptions.AFZ={'AFZ'};
        channelnameoptions.EEG1={'AF3','ECG'}; 
        channelnameoptions.EEG2={'AF4'};
        channelnameoptions.EEG3={'AF7'};
        channelnameoptions.EEG4={'AF8'};
        channelnameoptions.EEG5={'AF3'}; %unreferenced (actually referenced to AF7,AF8 average)
        channelnameoptions.EEG6={'AF4'};
        channelnameoptions.EEG7={'AF7'};
        channelnameoptions.EEG8={'AF8'};
        
        if isfield(settings, 'mouthEEG') && settings.mouthEEG == 1
            channelnameoptions.EEG9={'ECG_LF'};
            channelnameoptions.EEG10={'ECG_RA'};
        end
        
        channelnameoptions.LOC={'E1'}; %needs to be E1-E4
        channelnameoptions.ROC={'E2'}; %needs to be E2-E3
        channelnameoptions.REF1={'E3'};
        channelnameoptions.REF2={'E4'};
        channelnameoptions.EKG={'EKG'};
        channelnameoptions.EKG_LA={'ECG_LA'};
        channelnameoptions.EKG_LF={'ECG_LF'};
        channelnameoptions.EKG_RA={'ECG_RA'};
        channelnameoptions.Thorax={'Inductance_Thorax','Thorax_RIP'};
        channelnameoptions.Thorax_Fast={'Thorax_Fast'};
        channelnameoptions.Abdomen={'Inductance_Abdomen','Abdomen_RIP'};
        channelnameoptions.Abdomen_Fast={'Abdomen_Fast'};
        channelnameoptions.Flow={'Mask_Pressure','Nasal_Pressure'}; %'Vflow'
        channelnameoptions.AudioVolume={'Audio_Volume'};
        channelnameoptions.Pleth={'Pleth','Pulse_Wave__Pleth_'};
        channelnameoptions.Pulse={'Pulse'};
        channelnameoptions.Snore={'Snore'};
        channelnamestemp=fieldnames(channelnameoptions);
        
        for i=1:length(channelnamestemp)
            temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';'])
                        ChannelTitlesTemp{i} = char(temp(nn));
                        %eval([channelnamestemp{i} '.title=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        ChannelsFs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        %ChannelStartTime = NaN*ones(length(ChannelsList),1);
        Signals = filehandle.Signals;
        
        TimeInfo = filehandle.TimeInfo;
        
        %% code to import external audio file
        try
            Audiofilename = [Filenames{n,4} Filenames{n,1}(1:end-4) 'Snore.mat'];
            if exist(Audiofilename)==2
                w3 = matfile(Audiofilename);
                List = {'Snore'};
                
                for i=1:length(List)
                    try
                        eval([List{i} ' = w3.' List{i} ';']);
                        disp('successfully loaded external Audio signals')
                    catch me
                        disp('failed to load Audio signal')
                    end
                end
                
                % TEMPORARY
                if isfield(Snore,'fs')
                    Snore.data = Snore.data';
                    Snore.fsTrue = Snore.fs;
                end
                
                
                %                 SnoreChannelFound=zeros(length(List),1);
                %                 for i=1:length(List)
                %                     SnoreChannelFound(i)=(exist(List{i})==1)*1;
                %                 end
                %                 List = List(SnoreChannelFound==1);
                %
                %                 %eventually export this also into ConvertedInfo
                %                 try
                %                     SnoreInfo = mergestructs(w3.settings);
                %                 catch me
                %                     'fix SnoreInfo'
                %                 end
                %                 ChannelsList = [ChannelsList List];
                %                 i = find(strcmp(ChannelsList,'Flow'));
                %                 ChannelsFs = [ChannelsFs;repmat(ChannelsFs(i),length(List),1)];
                %                 ChannelFound = [ChannelFound;SnoreChannelFound];
            end
        catch
            disp('failed to import external Audio');
        end
        
        
        %get sampling rates
        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    %I=find(strcmp(channelnamestemp,ChannelsList{i}));
                    %I2=find(strcmp(Signals.labelsNoSpace,ChannelTitlesTemp{I}));
                    %ChannelsFs(i)=eval([ChannelsList{i} '.fs']);
                    ChannelsFs(i)=eval([ChannelsList{i} '.fsTrue']);
                catch me
                    ChannelsFs(i)=settings.Fs; %assumed
                    disp(['missing sampling rate assuming settings.Fs']);
                end
                %                 temp = eval([ChannelsList{i} '.data']);
                %                 try
                %                     ChannelTitlesoriginal{i}=ChannelTitlesTemp{i};
                %                 catch me
                %                     ChannelTitlesoriginal{i}=[];
                %                 end
                %                 if 1
                %                 eval(['clear ' ChannelsList{i}]);
                %                 eval([ChannelsList{i} '=temp;']);
                %                 end
            else
                % display to command window, but not to PUPbeta gui
                disp(['strewth, no ', ChannelsList{i}, ' found']);
            end
        end
        
        for i=1:length(ChannelsList)
            if ChannelFound(i)==1
                temp = eval([ChannelsList{i} '.data']);
                eval([ChannelsList{i} '=temp(:);']);
            end
        end
        
        %NoxPositionSignals
        PosLength = max([length(X_Axis) length(Y_Axis) length(Z_Axis)]);
        try
            if length(X_Axis)<PosLength
                nandiff = nan(PosLength-length(X_Axis),1);
                X_Axis = [X_Axis; nandiff];
            end
            if length(Y_Axis)<PosLength
                nandiff = nan(PosLength-length(Y_Axis),1);
                Y_Axis = [Y_Axis; nandiff];
            end
            if length(Z_Axis)<PosLength
                nandiff = nan(PosLength-length(Z_Axis),1);
                Z_Axis = [Z_Axis; nandiff];
            end
        catch me
        end
        
        Pos_AccMagnitude = sqrt(X_Axis.^2+Y_Axis.^2+Z_Axis.^2);
        %supine: Y=-1,X=0,Z=0; code:1/(pi/2)*atan2(0,1)
        %upright: Y=-0,X=0,Z=-1; code:1/(pi/2)*atan2(0,1)
        Position = 1/(pi/2)*atan2(X_Axis./Pos_AccMagnitude,-Y_Axis./Pos_AccMagnitude); %quarter turns
        %Incline = 1/(pi/2)*atan2(-Z_Axis./Pos_AccMagnitude,sqrt(Y_Axis.^2+X_Axis.^2)./Pos_AccMagnitude); %quarter turns
        Incline = -Z_Axis./Pos_AccMagnitude;
        
        I = find(strcmp(ChannelsList,'X_Axis')==1);
        %Fs_temp = ChannelsFs(I);
        ChannelFound(I+2)=0; %assumes X,Y,Zaxis are together
        ChannelsList{I}='Position';
        ChannelsList{I+1}='Incline';
        
        % Referencing Electrodes:
        disp('Referencing Electrodes [added 20200407, consider repeating prior conversions]');
        try
            tempref = (REF1+REF2)/2; %(E3+E4)/2
            EEG1 = EEG1 - tempref;
            EEG2 = EEG2 - tempref;
            EEG3 = EEG3 - tempref;
            EEG4 = EEG4 - tempref;
        catch me
            disp('failed EEG ref')
        end
        
        try
            LOC = LOC - REF2; %E1 - E4
            ROC = ROC - REF1; %E2 - E3
        catch me
            disp('failed EOG ref')
        end
   
   
  %%   
        
    end % end signal import
        
    ChannelsList(ChannelFound==0)=[];
    ChannelsFs(ChannelFound==0)=[];
    try
        ChannelTitlesoriginal(ChannelFound==0)=[];
    end
    try
        Info.ChannelTitlesoriginal=ChannelTitlesoriginal;
    end

    %% Bring in optional signals here
    % currently code is sitting near end of fn


    % 


    %% Replace Flow channel with profusion exported ASCII unfiltered nasal pressure channel
    % EAS added on 2021-12-15
    % looks for patient_nasal_ascii.txt files where patient is the patient
    % ID in AMasterSpreadsheet
    % See SyncResp folder for code if synchronisation of signals needed due
    % to different hardward recording devices
    if isfield(settings,'ProfusionASCIINasalPressureImport') && settings.ProfusionASCIINasalPressureImport==1
        importNasalASCII = false;
        nasalfname=[Filenames{n,1}(1:end-4) '_nasal_ascii.txt'];
        if exist([directory nasalfname],'file')==2
            importNasalASCII = true;
        else  
            disp('Warning: ProFusion Nasal Pressure ASCII file suspected but not used, please label as filename_nasal_ascii.txt')
        end
        if importNasalASCII
            displaytext = 'Replacing flow data'; disp(displaytext); 
            M = readmatrix([directory nasalfname]);
            oldFlow = -1 * Flow; % inverted because EDF nasal pressure is usually the opposite polarity to ASCII nasal pressure
            Flow = Flow * 0 + NaN;
            if (length(M) < length(oldFlow))
                Flow(1:length(M)) = M;
            else
                Flow = M;
            end
            % plot entire old flow and new flow signals
            fig = figure(11); 
            clf(fig); hold on;
            fig.WindowState = 'maximized';
            title([settings.Filenames{n}(1:end-4) ' Old Flow and Imported Nasal ASCII Flow', ' Entire']);
            plot(oldFlow);
            plot(Flow);
            legend('Old Flow', 'New Flow (Nasal ASCII)');
            exportgraphics(fig, [settings.ConvertedDirectory, settings.Filenames{n}(1:end-4), '_nasal_ascii_entire.png']);
            % plot the last minute of old flow and new flow signals
            flowIndex = find(strcmp(ChannelsList,'Flow') == 1);
            flowFs = ChannelsFs(flowIndex);
            index = length(Flow) - (flowFs * 60);
            if (index >= 1)            
                xlim([index, length(Flow)]);
                title([settings.Filenames{n}(1:end-4) ' Old Flow and Imported Nasal ASCII Flow', ' Last Minute']);
                exportgraphics(fig, [settings.ConvertedDirectory, settings.Filenames{n}(1:end-4), '_nasal_ascii_last_minute.png']);
            end
            % plot the first minute of old flow and new flow signals
            index = flowFs * 60;
            if (index <= length(Flow))
                xlim([1, index]);
                title([settings.Filenames{n}(1:end-4) ' Old Flow and Imported Nasal ASCII Flow', ' First Minute']);
                exportgraphics(fig, [settings.ConvertedDirectory, settings.Filenames{n}(1:end-4), '_nasal_ascii_first_minute.png']);
            end
            close(fig);
        end
    end
    
    %% add patient disconnects
    % after signal import, and before any further processing, add patient disconnect data.
    % Occassionally, data (e.g. Cincinatti Childrens Hospital) has periods with patient disconnects,
    % where the edf is on 'hold' while real time just keeps on trucking on, so
    % the exported signals are out of sync with annotations after this period
    if isfield(settings,'AddDisconnectPeriods') && settings.AddDisconnectPeriods==1
        % read disconnect data for study from corresponding text file
        % (usually there is only one disconnect per study, but not always...
        textfilename=[directory Filenames{n,1}(1:end-4) '_disconnect.txt'];
        if exist(textfilename,'file')==2
            displaytext = 'Adding patient disconnect periods';
            disp(displaytext);
%             set(handletext,'String',displaytext); drawnow;
            % get the disconnect information: (in text file, named 'yourstudyname_disconnect.txt')
            % first column is onsetTime, for eg. 5:25:21.35, as HH:MM:SS.ss from start of study, and
            % second column is holdTime, for e.g. 15.21 as SS.ss, and
            % each row is each disconnect.
            fid = fopen( textfilename );
            cac = textscan(fid,'%s%f','Delimiter','\t');
            fclose( fid );
            
            for disconnectNum=1:size(cac{1,1},1)
                % get the disconnect time (from start of study)
                % [~,~,~,theH,theM,theS]=datevec( cac{disconnectNum,1}, 'HH:MM:SS.FFF' );
                [~,~,~,theH,theM,theS]=datevec( cac{1,1}{disconnectNum,1}, 'HH:MM:SS.FFF' );
                % get the disconnect duration
                % holdTime = cac{disconnectNum,2};  %cac{i,2} is the duration of disconnect
                holdTime = cac{1,2}(disconnectNum);  %cac{i,2} is the duration of disconnect
                % now, loop through each signal/channel,
                % get the ChannelsFs to work out at what point to 'inject' a pause
                % work out how many samples the pause should be at Channel_Fs
                for ChannelLoop=1:length(ChannelsList)
                    CurrentFs = ChannelsFs(ChannelLoop); % get fs
                    holdSamples = floor(holdTime*CurrentFs); % how many samples in holdTime
                    % work out how many samples the hold period starts at
                    holdStartsAt = (theH*60*60*CurrentFs) + (theM*60*CurrentFs) + (floor(theS*CurrentFs));
                    % inject hold period into channel
                    % ideally NaN, but some signals are filtered, and NaN's
                    % go bad, so using zeros instead.
                    eval([ChannelsList{ChannelLoop} ' = [' ChannelsList{ChannelLoop} '(1:' num2str(holdStartsAt) '); zeros(' num2str(holdSamples) ',1); ' ChannelsList{ChannelLoop} '(' num2str(holdStartsAt) '+1:end)];' ]);
                end
            end
            Info.DisconnectPeriods = cac; % also save these periods into the Info struct for access during Analysis
        end
    end
    % the second stage of 'AddPatientDisconects' is during Analysis, where
    % we reject analysis from windows with disconnect periods.
    
    %% Calibrate any signal using text
    for i=1:length(ChannelsList)
        textfilename=[directory Filenames{n,1}(1:end-4) '_' ChannelsList{i} '_cal.txt'];
        if exist(textfilename,'file')==2
            if ~(exist(ChannelsList{1},'var')==1)
                continue
            end
            displaytext = ['Recalibrating: ' ChannelsList{i}];
            disp(displaytext); 
%             set(handletext,'String',displaytext); drawnow;
            [M] = dlmread(textfilename); %col1 = M(:,1); col2 = M(:,2);
            eval([ChannelsList{i} '=M(1).*' ChannelsList{i} '+M(2);']);
        end
    end
    
    %% Fix FsDrift Nox C1 vs A1 (make this generic in future)
    if strcmp(SignalFormat,'NoxMat') && exist('Pes')==1 && exist('Flow')==1
        FsErrorList = {'Flow','Pes'};
        [Fserror,delay] = FixFsDriftC1A1f(Flow,Pmask,ChannelsFs,ChannelsList)
        
        for i=1:length(FsErrorList)
            try
                I = find(strcmp(ChannelsList,FsErrorList{i}));
                ChannelsFs(I) = ChannelsFs(I)*(1+Fserror);
                if 1 %incorporated C1A1 delay also
                    TimeTemp=(0:(1/ChannelsFs(I)):0+(length(eval(FsErrorList{i}))-1)*(1/ChannelsFs(I)))'; % This is the time vector associated with the 100Hz Flow data.
                    temp = interp1(TimeTemp,eval(FsErrorList{i}),TimeTemp-delay,'pchip');
                    eval([FsErrorList{i} '= temp;']);
                end
            catch me
            end
        end
    end
    
    
    %% Time information
    %base the study duration on Flow signal (presumably present)
    N_Flow=length(Flow);
    Fs_Flow=ChannelsFs(find(strcmp(ChannelsList,'Flow')==1));
    N_timeXHz = round((N_Flow/Fs_Flow*settings.Fs));
    
    displaytext=['Get recording start time'];
    disp(displaytext);
    
    if strcmp(SignalFormat,'EDF')
        
        if isfield(settings,'BlockEDFload') && settings.BlockEDFload>0
            StartTimeText=Header.recording_starttime; %new from blockEDFLoad
        else
            fid = fopen([directory fname],'r');
            fseek(fid,176,-1);
            StartTimeText = char(fread(fid,8,'char')')
            fclose(fid); % Close file
        end
        
        try
            StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
            if StartTime<43200; StartTime=StartTime+86400; end
            disp(['Found StartTime: ' num2str(StartTime) ', ' StartTimeText]);
        catch me
            StartTime=0;
            disp('Failed to import EDF Start Time, set to 0 sec');
        end
        
    elseif strcmp(SignalFormat,'Spike') || strcmp(SignalFormat,'SpikeDise')
        
        StartTime=0;
        try 
            ChannelsStartNotNan = ChannelsStart(~isnan(ChannelsStart));
            StartTime=ChannelsStartNotNan(1);
        catch
            
        end
%         try
%             eval(['StartTime=filehandle.StarttimeSpike']);
%         catch me
%             %add CED function here in try-catch
%             %if StartTime<43200, StartTime=StartTime+86400; end
%             %
%             % Will need to load the matson library to use the CED functions
%         end
    elseif strcmp(SignalFormat,'NoxMat') || strcmp(SignalFormat,'NoxMatSAS') || strcmp(SignalFormat,'NoxMatT3')
        StartTime = TimeInfo.StartTime;
       
         
    end

    
    %% Create Time
    % Time array is based on desired sampling rate here:
    Time=(StartTime:(1/settings.Fs):StartTime+(N_timeXHz-1)*(1/settings.Fs))'; % This is the time vector associated with the _XHz Flow data.
    EndTime=Time(end);
    ChannelsList = ['Time',ChannelsList];
    ChannelsFs = [settings.Fs;ChannelsFs];
    
    %If needed for handling flow etc:
    TimeFlow=(StartTime:(1/Fs_Flow):StartTime+(N_Flow-1)*(1/Fs_Flow))'; % This is the time vector associated with the _XHz Flow data.
    
    
    %%
    if isfield(settings,'MergeWithLabChart') && settings.MergeWithLabChart==1
        'To do: write code to import labchart, then sync using flow, and interpolate signals'
        fnameadicht = [fname(1:find(fname=='.',1,'last')) 'adicht'];
        if exist([directory fnameadicht])==2
            disp( 'Found LabChart adicht file, importing signals')
        end
        run ImportAdichtAll2Struct; %Output T
        
        TimeFlow = StartTime + (0:(1/Fs_Flow):(1/Fs_Flow*(length(Flow)-1)))';
        TimeFlowSync = StartTimeLCeq + (0:(1/Fs):(1/Fs*(length(T.FlowSync)-1)))';
        
        FlowSyncds = interp1(TimeFlowSync,T.FlowSync,TimeFlow);
        %
        Lagaverage2=SyncFlowExact(Flow,FlowSyncds,TimeFlow,0.9);
        delay = Lagaverage2(1);
        %Time2 = TimeFlow + Lagaverage2; %corrected Time is not used
        
        Fserror = ((1/Fs_Flow)/median(diff(TimeFlow + Lagaverage2))) - 1; %Fserror = (dt/dt2) - 1; FsNew = FsOld*(1+Fserror)
        
        TimeFlowSync2 = StartTimeLCeq + delay + (0:(1/(Fs*(1+Fserror))):(1/(Fs*(1+Fserror))*(length(T.FlowSync)-1)))';
        %check:
        %LagaverageCheck=SyncFlowExact(Flow,interp1(TimeFlowSync2,T.FlowSync,TimeFlow),TimeFlow,0.9);
        
        %make new signals that will sync with EDF StartTime
        Tfields = fieldnames(T);
        FsNew = Fs*(1+Fserror);
        %TimeFlowSync3 = TimeFlowSync2 + (StartTime - TimeFlowSync2(1));
        TimeTemp = (StartTime:(1/FsNew):TimeFlow(end))';
        for i=1:length(Tfields)
            signaltemp = getfield(T,Tfields{i});
            temp = interp1(TimeFlowSync2,signaltemp,TimeTemp,'pchip');
            temp(TimeTemp<TimeFlowSync2(1))=0;
            temp(TimeTemp>TimeFlowSync2(end))=0;
            T = setfield(T,Tfields{i},temp);
        end
        if 0
            figure(99); clf(99);
            TimeTemp=StartTime + (0:(1/FsNew):0+(length(T.FlowSync)-1)*(1/FsNew))';
            ax99(1)=subplot(2,1,1); plot(TimeFlow,Flow);
            %             ax99(2)=subplot(2,1,2); plot(TimeFlowSync2,T.FlowSync);
            ax99(2)=subplot(2,1,2); plot(TimeTemp,T.FlowSync);
            linkaxes(ax99,'x');
        end
        
        %Warning do not call things <varname>LC
        Tfields = fieldnames(T);
        isVariable=[];
        for i=1:length(Tfields)
            isVariable(i)=exist(Tfields{i});
            if isVariable(i)==1
                newfield = [Tfields{i} 'LC'];
                T = setfield(T,newfield,getfield(T,Tfields{i}));
                T = rmfield(T,Tfields{i});
            end
        end
        Tfields = fieldnames(T);
        for i=1:length(Tfields)
            isVariable(i)=exist(Tfields{i});
            if isVariable(i)==0
                eval([Tfields{i} ' = T.' Tfields{i} ';']);
            end
        end
        clear T
        
        ChannelsList = [ChannelsList Tfields'];
        ChannelsFs = [ChannelsFs;repmat(FsNew,length(Tfields),1)];
        
        if 1
            Flow = FlowSyncds;
        end
    end
    
    %%
    %figure(99); plot(TimeFlow,Flow)
    
    %% Flow Sanity Qc check
    try
        if ~(exist('Info')==1)
            Info=struct();
        end
        settings.fname=fname;
        % 1. find if flow is inverted & overall noise level
        [PrUpright,FnoiseAll,Info]=FlowInvertedDetectTool(Flow,TimeFlow,Info);
        Info.FlowQ.PrUpright = PrUpright;
        Info.FlowQ.FnoiseAll = FnoiseAll;
        
        if FnoiseAll>=0.1
            disp('noisy or absent flow signal')
        end
        
        
        % 2. find if a high pass/low pass filter is applied to flow signal
        verbose=0;
        ploton=0;
        FlowRS = interp1(TimeFlow,Flow,Time,'linear');
        FlowFilterDetect = FlowFilterDetector(FlowRS,Time,ploton,verbose);
        
        clear FlowRS;
        Info.FlowQ.FlowFilterDetect = FlowFilterDetect;
        
        try 
            if FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1)<8
                disp('Warning: Flow appears smoothed or downsampled');
            end

            if FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1)>0.011
                disp('Warning: Flow appears distorted by baseline adjustment (high pass)');
            end
        end
        
        % 3. find if flow is clipped or not
        [~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[]);
        
        Ttemp.FclippedI=FclippedI;
        Ttemp.FclippedE=FclippedE;
        Ttemp.FclippedTotal=FclippedE+FclippedI;
        
        Info.FlowQ.FclippedTotal = Ttemp.FclippedTotal;
        
        if Ttemp.FclippedTotal>0.005
            disp('Warning: Flow appears clipped');
            disp(['   Clipping fraction: ' num2str(100*Ttemp.FclippedTotal,2) ' %']);
            
        end

    catch
        disp('FAILED QC PROCESS')
    end
    
    
    %%
    if strcmp(system,'SandmanPlusLabview') %where Pmask is the hidden sync channel:
        I=find(strcmp(ChannelsList,'Pmask'));
        if ~isempty(I)
            ChannelsList{I}='SynchSM';
            eval([ChannelsList{I} '=Pmask;']);
            clear Pmask
        end
        [AllSigLbView]=LabviewReadAllChannels(settings,n,Filenames,SynchSM);
        % 1	Time 2-sync_DAQ 3-BP 4-CO2 5-ETCo2 6-blank 7-TCD 8-Sao2 w/f 9-% Sao2
        % 10-13-raw fnirs data1; 14-17: raw data2; 18-20:channels
        % removed; 21-24: hbo; 25-28:hbr; 29-avg hbo; 30-avg hbr;
        if settings.fnirs
            ChannelsListLbView={'SynchLbView','BPWave','Co2Wave','EtCo2','fNirsMarker','CBFVWave','SaO2WaveLbView','SaO2perLbView',...
                'Raw690_1','Raw690_2','Raw690_3','Raw690_4',...
                'Raw830_1','Raw830_2','Raw830_3','Raw830_4',...
                'RemvdChFlag1','RemvdChFlag2','RemvdChFlag3',...
                'Hbo1','Hbo2','Hbo3','Hbo4','Hb1','Hb2','Hb3','Hb4',...
                'AvgHbO','AvgHbR'};
            
            ChannelsNumbersLbView=[2:30];
        else
            ChannelsListLbView={'SynchLbView','BPWave','Co2Wave','EtCo2','fNirsMarker','CBFVWave','SaO2WaveLbView','SaO2perLbView'};
            ChannelsNumbersLbView=[2:9];
        end
        %         LbviewFs=1/(AllSigLbView(2,1)-AllSigLbView(1,1)); % looks weird. double check?
        ChannelsList=[ChannelsList,ChannelsListLbView];
        for i=1:length(ChannelsListLbView)
            eval([ChannelsListLbView{i} '=AllSigLbView(:,ChannelsNumbersLbView(i));']);
            ChannelsFs=[ChannelsFs;128];
        end
        clear AllSigLbView;
        if nanmean(Co2Wave)<2 && prctile(Co2Wave,95)<2
            Co2Wave=Co2Wave.*100;
            EtCo2=EtCo2.*100;
        end
        if nanmean(CBFVWave)<2 && prctile(CBFVWave,95)<2
            CBFVWave=CBFVWave.*100;
        end
        if nanmean(SaO2WaveLbView)<2 && prctile(SaO2WaveLbView,95)<2
            SaO2WaveLbView=SaO2WaveLbView.*100;
            SaO2perLbView=SaO2perLbView.*100;
        end
        
        Pos_indx = find(strcmp(ChannelsList,'Position')==1);
        if isempty(Pos_indx) % if it's not there for some reason
            try
                Position=getPosSandman(Filenames,Time,StartTime,settings);
                ChannelsList = [ChannelsList,'Position'];
                ChannelsFs = [ChannelsFs;settings.Fs];
            catch
                Position = Time*0 + NaN;
                ChannelsList = [ChannelsList,'Position'];
                ChannelsFs = [ChannelsFs;settings.Fs];
            end
        end
    end
    
    %% If exist('BPwave') generate SBP / DBP / MAP
    if exist('BPWave')
        Fs_BP=ChannelsFs(find(strcmp(ChannelsList,'BPWave')));
        SysDiaTrace=[];
        [SysDiaTrace(:,1),SysDiaTrace(:,2),SysDiaTrace(:,3)]=BPFeatures(BPWave,Time,ChannelsList,ChannelsFs);
        ChannelsListBP={'SysBP','DiaBP','Map'};
        ChannelsNumbersBP=[1,2,3];
        ChannelsList=[ChannelsList,ChannelsListBP];
        for i=1:length(ChannelsListBP)
            eval([ChannelsListBP{i} '=SysDiaTrace(:,ChannelsNumbersBP(i));']);
            ChannelsFs=[ChannelsFs;Fs_BP];
        end
        
        clear SysDiaTrace
        
        if nanmean(BPWave)<3 && prctile(BPWave,95)<3
            BPWave=BPWave.*100;
        end
    end
    %
    
    %% If exist('CBFVwave') generate SBF / DBF / MBF --Blood flow signal    
    if exist('CBFVWave')
        Fs_CBFV=ChannelsFs(find(strcmp(ChannelsList,'CBFVWave')));
        SysDiaTraceCBFV=[];
        [SysDiaTraceCBFV(:,1),SysDiaTraceCBFV(:,2),SysDiaTraceCBFV(:,3)]=CBFVFeatures(CBFVWave,Time,ChannelsList,ChannelsFs);
        ChannelsListCBFV={'SysCBFV','DiaCBFV','MeanCBFV'};
        ChannelsNumbersCBFV=[1,2,3];
        ChannelsList=[ChannelsList,ChannelsListCBFV];
        for i=1:length(ChannelsListCBFV)
            eval([ChannelsListCBFV{i} '=SysDiaTraceCBFV(:,ChannelsNumbersCBFV(i));']);
            ChannelsFs=[ChannelsFs;Fs_CBFV];
        end
        
        clear SysDiaTraceCBFV
        
    end
    
    
    %% Check SpO2 is 0-100 not 0-1
    if exist('SpO2', 'var') % Dan added this because DISE data does not have SpO2
        temp = SpO2;
        temp(temp<0.4)=NaN;
        if nanmedian(temp)<1 && prctile(temp,95)<1 && prctile(temp,95)>0.5 && nanmedian(temp)>0.5
            disp('SpO2 is from 0-1, multiplying by 100');
            SpO2=SpO2*100;
        end
    end
    
    %% EEG length check
    %all EEGs must be same length
    %Make channels the same length (sometimes these are off by a few samples)
    disp('Starting EEG analysis');
    EEGList = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12'};
    for i=1:length(EEGList)
        try
            EEGlength(i)=length(eval(EEGList{i}));
        catch me
            EEGlength(i)=NaN;
        end
    end
    EEGlength = round(nanmedian(EEGlength));
    for i=1:length(EEGList)
        try
            while length(eval([EEGList{i} ]))~=EEGlength
                if length(eval([EEGList{i} ]))<EEGlength
                    eval([EEGList{i} '(end+1)=' EEGList{i} '(end);']); %add a sample to the end if the channel is too short by 1
                elseif length(eval([EEGList{i} ]))>EEGlength
                    eval([EEGList{i} '(end)=[];']); %delete a sample from the end if the channel is too long by 1
                end
            end
        catch me
        end
    end
    clear EEGList EEGlength
    
    if settings.EverySecondEEGisRef(n) %up to 6 EEG signals (12 leads)
        if exist('EEG1','var')&&exist('EEG2','var')
            EEG1=EEG1-EEG2;
        end
        clear EEG2;
        Itemp=find(strcmp(ChannelsList,'EEG2')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
        if exist('EEG3','var')&&exist('EEG4','var')
            EEG3=EEG3-EEG4;
        end
        clear EEG4;
        Itemp=find(strcmp(ChannelsList,'EEG4')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
        if exist('EEG5','var') && exist('EEG6','var')
            EEG5=EEG5-EEG6;
        end
        clear EEG6;
        Itemp=find(strcmp(ChannelsList,'EEG6')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
        if exist('EEG7','var')&&exist('EEG8','var')
            EEG7=EEG7-EEG8;
            
        end
        clear EEG8;
        Itemp=find(strcmp(ChannelsList,'EEG8')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
        if exist('EEG9','var')&&exist('EEG10','var')
            EEG9=EEG9-EEG10;
        end
        clear EEG10;
        Itemp=find(strcmp(ChannelsList,'EEG10')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
        if exist('EEG11','var')&&exist('EEG12','var')
            EEG11=EEG11-EEG12;
        end
        clear EEG12;
        Itemp=find(strcmp(ChannelsList,'EEG12')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
    end
    
    %Automatically subtracts EKG signals if two EKG channels are provided
    if exist('EKG')&&exist('EKG2')
        EKG=EKG-EKG2;
        clear EKG2;
        Itemp=find(strcmp(ChannelsList,'EKG2')==1);
        ChannelsList(Itemp)=[];
        ChannelsFs(Itemp)=[];
    end
    
    %% Filter EEG if settings say so
    if isfield(settings,'FilterEEG')
        
        EEGList = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12'};
        
        count=1;
        EEGi=[];
        for i=1:length(EEGList)
            I = find(strcmp(ChannelsList',EEGList(i)));
            if ~isempty(I)
                EEGi(count)=I;
                count=count+1;
            end
        end
        
        filter_LFcutoff_butter0 = settings.FilterEEG(1);
        filter_HFcutoff_butter0 = settings.FilterEEG(2);
        filter_order0 = 2;
        
        %backup = EEG1;
        
        for i=1:length(EEGi)
            [B_butter0,A_butter0] = butter(filter_order0,[filter_LFcutoff_butter0 filter_HFcutoff_butter0]/(ChannelsFs(EEGi(i))/2));
            temp = eval(ChannelsList{EEGi(i)});
            temp2 = filtfilt(B_butter0,A_butter0,temp); %filtfilt, otherwise flow signal is right-shifted
            eval([ChannelsList{EEGi(i)} '=temp2;']);
        end
        
        %     figure(2)
        %     plot((0:(1/ChannelsFs(EEGi(1))):0+(length(backup)-1)*(1/ChannelsFs(EEGi(1)))),[backup EEG1])
        %
    end
    
    %% EEG multiplier
    if isfield(settings,'EEGmultiplier')
        
        EEGList = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12','LOC','ROC'};
        count=1;
        EEGi=[];
        for i=1:length(EEGList)
            I = find(strcmp(ChannelsList',EEGList(i)));
            if ~isempty(I)
                eval([EEGList{i} '=' EEGList{i} '*settings.EEGmultiplier;']);
            end
        end
    end
    
    %% Events and Hypnogram
    % moved up in case it changes StartTime
    try
        
        % DLM: should add something here to skip annotations/epochs file
        % finding, and EventsAndHypnogram code for when running EDF only
        
        if exist([Filenames{n,5} Filenames{n,2}])~=2
            disp(['   *** Warning: Can not find annotations file: ' [Filenames{n,5} Filenames{n,2}]] )
        end
        if exist([Filenames{n,6} Filenames{n,3}])~=2
            disp(['   *** Warning: Can not find epochs file: ' [Filenames{n,6} Filenames{n,3}]] )
        end
        
        
        %[Evts,StartTimeNew] =
        % EventsAndHypnogram(system,Filenames(n,:),StartTime,EndTime,handletext); handletext doesn't exist dso removed here, SS: 7/16/2021
        disp ('Starting Events and Hypnogram');
        [Evts,StartTimeNew] = EventsAndHypnogram(system,Filenames(n,:),StartTime,EndTime,[],settings); %% Epochs and events
        
        Info.StartTimeInfo.StartTimeOriginal=StartTime;
        Info.StartTimeInfo.StartTime=StartTime;
        
        if StartTimeNew~=StartTime
            Time = Time - Time(1) + StartTimeNew; %update Time array
            StartTime = StartTimeNew; %update for below
            EndTime = Time(end); %not used, but just good practice.
            Info.StartTimeInfo.StartTimeFromEandH=StartTimeNew; %overwrite
            Info.StartTimeInfo.StartTime=StartTimeNew; %overwrite
        end
    catch
        if ~exist('Evts')
            Evts=struct();
        end
    end
    %TimeinsecToDatetime(Evts.Table1.EventStart)
    
    %% Correct StartTime For Signals Manually (run exactly once)
    % Note this is for cases where start time is wrong relative to annotations files
    % Do not use for adding a start time to both signals and events.
    
    textfilename=[directory Filenames{n,1}(1:end-4) '_StartTime.txt'];
    if exist(textfilename,'file')==2
        displaytext = 'Correcting StartTime For Signals';
        disp(displaytext); 
%         set(handletext,'String',displaytext); drawnow;
        %Position.values_original=Position.values;
        %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
        % textread was not recommended, upgraded to dlmread
        [M] = dlmread(textfilename); col1 = M(:,1);
        StartTimeOrig = StartTime;
        StartTime = col1;
        Info.StartTimeInfo.StartTimeTextOverwrite=StartTime; %overwrite
        Info.StartTimeInfo.StartTime=StartTime; %overwrite
        %StartTime = StartTime - col1;
        if abs(Time(1)-StartTime)>0.001
            Time = Time + StartTime-Time(1); %update Time array
            EndTime = Time(end); %not used, but just good practice.
        end
    end
    
    %% Sampling Frequency Errors, Time Drift
    % settings.FsLeak=-0.00008; %lag of EDF time vs actual in proportion
    if settings.FsLeak~=0
        disp('Warning: fixing time leak');
        ChannelsFs2 = ChannelsFs;
        ChannelsFs_fix = ChannelsFs*0 + 1;
        I=find(sum(string(ChannelsList)==string({'Epochs','EventsAr','EventsResp'})'));
        ChannelsFs_fix(I)=0;
        ChannelsFs2(ChannelsFs_fix==1) = ChannelsFs2(ChannelsFs_fix==1)*(1-settings.FsLeak);
        ChannelsFs = ChannelsFs2;
    end
    
    
    %% Profusion Position Manual Scoring/Editing Import
    importposfromprofusion=0;
    posfname=[Filenames{n,1}(1:end-4) '_posP.txt'];
    if exist([directory posfname],'file')==2
        importposfromprofusion=1;
    elseif isfield(settings,'ProfusionPositionUpdate') && settings.ProfusionPositionUpdate && exist([directory Filenames{n,1}(1:end-4) '.txt'],'file')==2
            % ProfusionPositionUpdate must be off by default otherwise this
            % will conflict with other system event import routines. 
            disp('Warning: Manual Profusion position file suspected but not used, please label as filename_posP.txt')
            % importposfromprofusion=1; %disabled
    end
    if importposfromprofusion
        %  posfname=[Filenames{n,1}(1:end-4) '.txt'];
        %  if exist([directory posfname],'file')==2
        if ~(exist('Position','var')==1)
            Position = Time*0 + NaN;
        end
        displaytext = 'Correcting Profusion position data'; disp(displaytext); 
        %Position.values_original=Position.values;
        %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
        % textread was not recommended, upgraded to dlmread,
        % starting with row offset of 1 because it has headers
        % column offset of 0.
        [M] = dlmread([directory posfname],',',1,0);
        
        % readtable also works, but is much slower, and will
        % have to convert column names if they dont meet reqs
        %summarytable = readtable([directory posfname]);
        
        % back up the old position data, because we overwrite it
        Position_old = Position;

        %settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
        % convert this epoch length data into time length data
        NewBodyPos = M(:,4); % column 4 is the corrected position
        Position_t = StartTime+(0:30:(30*(length(NewBodyPos)-1)))';
        Position = interp1(Position_t+15,NewBodyPos,Time,'nearest','extrap');
        
        TimeTemp=(StartTime:(1/16):StartTime+(length(Position_old)-1)*(1/16))'; % This is the time vector associated with the 100Hz Flow data.
        %SS query: why is dt hardcoded here as 16 Hz for TimeTemp?
        Temp = interp1(TimeTemp,Position_old,Time,'linear');
        Position_old_fulllength = Temp;
        
        figure(10); clf(figure(10));
        fig = gcf;
        fig.Color = [1 1 1];
        fig.Units = 'Inches';
        %fig.Position = [ 20.5 2 15 8];
        fig.Position = [ 2 2 10 6];
        
        plot(Position_old_fulllength+0.02, 'r-'); hold on;
        plot(Position-0.02, 'b-');
        legend('old', 'new');
        box off
        
        %for plotting:
        yticknumall = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n})';
        yticktextall = settings.poscodesdatabase(1,2:end)';
        PosTable=table(yticknumall,yticktextall);
        PosTable.Properties.VariableNames={'PosCode','Pos'};
        PosTable=sortrows(PosTable);
        PosTable(isnan(PosTable{:,1}),:)=[];

        yticks(PosTable.PosCode); %[0:1:4]
        yticklabels(PosTable.Pos); % {'Right','Supine','Left','Prone','Upright'} note these labels apply to the newer grael model compumedics output

        ylim([min(PosTable.PosCode)-0.1 max(PosTable.PosCode)+0.1]);
        grid on
        
        xlimits = xlim();
        xlim([-100000 xlimits(2)]);
        xticks([0 length(Position)])
        
        titlestr_=['Position data for ', Filenames{n}(1:end-4)];
        titlestr = regexprep(titlestr_,'_','');
        title(titlestr);
        
        %savestr = ['C:\Users\dwayne\OneDrive\OneDrive - The University of Queensland\Projects\PuPwithCorrectedPosition\Figures\', titlestr];
        %print(fig, savestr, '-dtiff', '-r300');
        
        saveas(fig, settings.PosFigName, 'png');
        close(figure(10));
        
        Pos_indx = find(strcmp(ChannelsList,'Position')==1);
        if isempty(Pos_indx) % if it's not there for some reason
            ChannelsList = [ChannelsList,'Position'];
            ChannelsFs = [ChannelsFs;settings.Fs];
        else % update what's there
            ChannelsFs(Pos_indx) = settings.Fs;
        end
    end
    
    %% should be a position signal by here, warn user if not
    if ~exist('Position','var')
        disp('*** Warning, no position signal, unclear if convert process will finish ***')
    end
    
    %% Events and Hypnogram was here
    
    %% Use PositionT table from events import (if exists). Overwrites Position signal data.
    % Uses the "Default" codes 1=Supine,2=Left,3=Right,4=Prone,5=Unknown,6=Upright
    % SystemPos in settings.protocol(n)='Default' to use this
    if isfield(Evts,'PositionT')
        disp('loading Position information from Position Events List:')
        if exist('Position')==1 && ~exist('PositionOrig')
            PositionOrig = Position;
        end
        
        %Position = .... PositionT
        Position = 0*Time + 5; %default 5 (unknown)
        for i=1:height(Evts.PositionT)
            I = Time>Evts.PositionT.Time(i);
            Position(I) = Evts.PositionT.Codes(i);
        end
        
        figure(10101); clf(10101);
        set(gcf,'color',[1 1 1]);
        %plot(Time,Position);
        plot(Time/86400,Position);
        datetick('x','HH:MMPM');
        set(gca,'ytick',[1:6],'yticklabels',{'Supine','Left','Right','Prone','Unknown','Upright'},'tickdir','out','box','off');
        
        if sum(ChannelsList=="Position")==0
            ChannelsList = [ChannelsList {'Position'}];
            ChannelsFs = [ChannelsFs;ChannelsFs(1)]; %Using sampling rate of Time, which comes from Flow
        else
            ChannelsFs(find(ChannelsList=="Position"))=settings.Fs; %Position has new sampling rate. Fixed 5/21/2021.
        end
    end
    
    %% Position channel fix using text
    % Fix Pos data if necessary
    % assumes times are in time since start recording in sec
    % e.g. for System Pos = 'Default', Supine = 1, Left = 2, Right = 3
    % Time is in native units as per converted file, see ManualPositionTool.xlsx
    textfilename=[directory Filenames{n,1}(1:end-4) '_pos.txt'];
    if exist(textfilename,'file')==2
        if ~(exist('Position','var')==1)
            Position = Time*0 + NaN;
        end
        displaytext = 'Correcting position data';
        disp(displaytext);
%         set(handletext,'String',displaytext); drawnow;
        %Position.values_original=Position.values;
        %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
        % textread was not recommended, upgraded to dlmread
        [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
        
        reftime = Time(1); %original code uses time since start recording in seconds as default
        if isfield(settings,'PositionCorrectionTimeIsSinceStartRecording') && settings.PositionCorrectionTimeIsSinceStartRecording==0
            reftime = 0;
        end
        
        for i=1:(length(col1)-1)
            Position(Time-reftime>=col1(i)&Time-reftime<col1(i+1))=col2(i);
        end
        Position(Time-reftime>=col1(end))=col2(end);
        
        if sum(ChannelsList=="Position")==0
            ChannelsList = [ChannelsList {'Position'}];
            ChannelsFs = [ChannelsFs;settings.Fs]; %Using sampling rate of Time, which comes from Flow
        else
            ChannelsFs(find(ChannelsList=="Position"))=settings.Fs; %Position has new sampling rate. Fixed 5/21/2021.
        end
        if 0
            figure(99);plot(Time,Position)
            %sum(diff(Position)~=0);
            %length(M)-1;%should be equal to sum(diff(Position)~=0);
        end
    end
    
    
    %% Plot Hypnogram
    if 0 %Plot Hypnogram
        
        figure(2);
        
        xtime=(Time(1):30:(Time(1)+30*(length(Evts.Hypnogram)-1)));
        
        ax2(1)=subplot(2,1,1);
        stairs(xtime/86400,Evts.Hypnogram);
        %datetick('x','HH:MMPM');
        
        ax2(2)=subplot(2,1,2);
        stairs(Evts.Hypnogram_t/86400,Evts.Hypnogram);
        %datetick('x','HH:MMPM');
        
        %             linkaxes(ax2,'x');
    end
    
    
    %% Make Continuous Epochs,EventsResp,EventsAr
    try
        EpochLength=30; %default
        
        if ~isfield(Evts,'Hypnogram')
            Evts.Hypnogram_t = [Time(1):30:Time(end)]';
            Evts.Hypnogram = Evts.Hypnogram_t*0+8;
        end

        if ~isnan(Evts.Hypnogram)  % Dan added to deal with converting data without event scoring
            EpochLength=Evts.Hypnogram_t(2)-Evts.Hypnogram_t(1);
            Epochs = interp1(Evts.Hypnogram_t+EpochLength/2,Evts.Hypnogram,Time,'nearest','extrap');
            % any data before Evts.Hypnogram_t(1) == 8;
            Epochs(Time<Evts.Hypnogram_t(1)) = 8;
            % any data after Evts.Hypnogram_t(end)+30 == 8;
            Epochs(Time>(Evts.Hypnogram_t(end)+EpochLength)) = 8;
        else
            disp('warning: no sleep stages detected');
            Epochs=NaN(size(Time));
        end
        
        EventsResp=0*Time;
        try
        for m=1:length(Evts.RespT.EventStart) %Resp Events
            lefti=round((Evts.RespT.EventStart(m)-StartTime)*settings.Fs)+1;
            righti=lefti+round((Evts.RespT.EventDuration(m))*settings.Fs);
            if lefti<1, lefti=1; end
            if righti>length(Time), righti=length(Time); end
            if Evts.RespT.EventCodes(m)>1
                EventsResp(lefti:righti)=Evts.RespT.EventCodes(m);
            end
        end
        end
        
        EventsAr=0*Time;
        try
        for m=1:length(Evts.ArT.EventStart) %Arousals
            lefti=round((Evts.ArT.EventStart(m)-StartTime)*settings.Fs)+1;
            righti=lefti+round((Evts.ArT.EventDuration(m))*settings.Fs);
            if lefti<1, lefti=1; end
            if righti>length(Time), righti=length(Time); end
            EventsAr(lefti:righti)=1;
        end
        end
        
        ChannelsList = [ChannelsList,{'Epochs','EventsAr','EventsResp'}];
        ChannelsFs = [ChannelsFs; settings.Fs;settings.Fs;settings.Fs];
        
        %% Evts
        % if an empty Evts struct is made above (i.e. drop out from EventsAndHypnogram function,
        % then we fail at this point because no fields for RespT or ArT.
        % DLM added line below to make empty structs for RespT and ArT        
        if ~isfield(Evts,'RespT')
            Evts.RespT.EventStart = [];
            % Evts.RespT.EventCodes = []; not sure if this is req'd
        end
        if ~isfield(Evts,'ArT')
            Evts.ArT.EventStart = [];
        end
        
        Evts.RespT.starttimesi=round((Evts.RespT.EventStart-StartTime)*settings.Fs+1);
        Evts.ArT.starttimesi=round((Evts.ArT.EventStart-StartTime)*settings.Fs+1);
        
        I = Evts.RespT.starttimesi > length(Time); %Time(end)*(1/(Time(2)-Time(1)));
        if sum(I)>0
            disp('warning: events scored after end of recording');
        end
        Evts.RespT(I,:)=[];
        
        I = Evts.ArT.starttimesi > Time(end)*(1/(Time(2)-Time(1)));
        if sum(I)>0
            disp('warning: arousals scored after end of recording');
        end
        Evts.ArT(I,:)=[];
       
        %% check AHI and TST are as expected
        if ~isempty(Evts.RespT.EventStart) | ~(all(Evts.Hypnogram==8) | isnan(Evts.Hypnogram)) % this was "if 1", DLM changed to skip this if Evts are non-existent; DV added nan condition
            [AHI.Total,Evtstemp,~]=getAHIEvtSubset(Evts,Epochs,Time,0*Time,1+0*Time);
            displaytext=['Total AHI: ' num2str(AHI.Total(58),4) ' per hr'];
            disp(displaytext); 
%             set(handletext,'String',displaytext); drawnow;
            displaytext=['Total Sleep Time: ' num2str(sum(Evtstemp.Hypnogram<4)*EpochLength/60) ' min'];
            disp(displaytext);
%             set(handletext,'String',displaytext); drawnow;
        else % DLM added this else statement for when using EDF without scoring
            AHI.Total=[];
            displaytext='Not checking AHI and TST because no events or hypnog found';
            disp(displaytext);
        end
        
    catch EvtsAndHypnog_fail
        if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
            keyboard
        else
            disp(' ************************ Failed Events/Hypnogram Routine ************************ ');
            if isfield(settings,'verbose') && settings.verbose==1
                disp(EvtsAndHypnog_fail.getReport);
            else
                disp(EvtsAndHypnog_fail.message);
            end
        end
    end
    
    %% Plot Debug if needed
    if 0
        figure(99);
        plot(TimeFlow,(Flow-nanmean(Flow))./nanstd(Flow));
        hold on
        plot(Time,EventsResp);
    end
    
    %% code to condition the Stimulation Amplitude signal (specific for DVs HGNS data)
    % load StimAmpl signal
    if strcmp(SignalFormat,'EDF') % DLM added this get out of jail free card
        % don't attempt to do StimAmp stuff is we are dealing with EDF conversion
    else
    % TRANSPORTED CODE TO STIM ON STUFF        
    end
    
    %% EEG power analysis, EKG
    displaytext=['EEG processing'];
    disp(displaytext);
%     set(handletext,'String',displaytext); drawnow;
    
    EEGsignals = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12'};
    temp=[];
    for i=1:length(EEGsignals)
        if exist(EEGsignals{i},'var')
            temp{length(temp)+1}=EEGsignals{i};
        end
    end
    EEGsignals=temp; %overwrite
    
    if ~isempty(EEGsignals)
        %Sampling rate of each EEG
        Fs_EEG=[];
        for i=1:length(temp)
            Fs_EEG(i)=ChannelsFs(find(strcmp(ChannelsList,EEGsignals{i})==1));
        end
        Fs_EEG_ = mode(Fs_EEG); %average EEG sampling rate
        
        % build in option to resample or discard EEG signals that have a
        % different sampling rate to the average
        % Resample channels
        
        TimeEEG=(StartTime:(1/Fs_EEG_):StartTime+(length(eval(EEGsignals{1}))-1)*(1/Fs_EEG_))'; % This is the time vector associated with the 100Hz Flow data.
        
        for i=1:length(EEGsignals)
            if Fs_EEG(i)~=Fs_EEG_
                displaytext=['Resampling: ' EEGsignals{i} ' from ' num2str(Fs_EEG(i)) ' to ' num2str(Fs_EEG_) 'Hz'];
                disp(displaytext); 
%                 set(handletext,'String',displaytext); drawnow;
                
                %eval([EEGsignals{i} ' = resample(' EEGsignals{i} ',Fs_EEG_,ChannelsFs(i));']); %only works if Fsamp / Fs are integer multiples of each other
                TimeTemp=(StartTime:(1/Fs_EEG(i)):StartTime+(length(eval(EEGsignals{i}))-1)*(1/Fs_EEG(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(EEGsignals{i}),TimeEEG,'linear');
                eval([EEGsignals{i} '= Temp;']);
                Fs_EEG(i)=Fs_EEG_;
                ChannelsFs(find(strcmp(EEGsignals{i},ChannelsList)))=Fs_EEG_;
            end
        end
        
        TimeEEG=(StartTime:(1/Fs_EEG_):StartTime+(length(eval(EEGsignals{1}))-1)*(1/Fs_EEG_))'; % This is the time vector associated with the 100Hz Flow data.
        
        fftlengthoptions = 2.^[6:9];
        [~,bestoptioni] = min(abs(Fs_EEG_ - fftlengthoptions));

        if exist('EKG','var')
            %EKG_original=EKG; %not needed
            Fs_EKG = ChannelsFs(find(strcmp(ChannelsList,'EKG')==1));
            if Fs_EKG~=Fs_EEG_
                displaytext=['Resampling EKG from ' num2str(Fs_EKG) ' to ' num2str(Fs_EEG_) 'Hz'];
                disp(displaytext);
%                 set(handletext,'String',displaytext); drawnow;
                %EKGeeg = resample(EKG,Fs_EEG_,Fs_EKG); %untested
                TimeTempEKG=(StartTime:(1/Fs_EKG):StartTime+(length(EKG)-1)*(1/Fs_EKG))'; % This is the time vector associated with the 100Hz Flow data.
                EKGeeg = interp1(TimeTempEKG,EKG,TimeEEG,'linear','extrap');
            else
                EKGeeg = EKG;
            end
            signalslist=[EEGsignals,'EKGeeg'];
        else
            signalslist=[EEGsignals];
        end
        
        %Make channels the same length (sometimes these are off by a few samples)
        for i=1:length(signalslist)
            if ~isempty(signalslist{i})
                while length(eval([signalslist{i} ]))~=length(TimeEEG)
                    if length(eval([signalslist{i} ]))<length(TimeEEG)
                        eval([signalslist{i} '(end+1)=' signalslist{i} '(end);']); %add a sample to the end if the channel is too short by 1
                    elseif length(eval([signalslist{i} ]))>length(TimeEEG)
                        eval([signalslist{i} '(end)=[];']); %delete a sample from the end if the channel is too long by 1
                    end
                end
            end
        end
        
        BetaPowerSettings.polyorder=1;
        BetaPowerSettings.suf1 = ''; %'.values'
        BetaPowerSettings.timeeegstr = 'TimeEEG';
        BetaPowerSettings.timestr = 'Time';
        BetaPowerSettings.flowstr = 'Flow';
        BetaPowerSettings.spo2str = 'SpO2';
        BetaPowerSettings.arstr = 'EventsAr';
        BetaPowerSettings.epochsstr = 'Epochs';
        BetaPowerSettings.fft_length = fftlengthoptions(bestoptioni); %256 %Best: 1710,1723:256, 1343:512
        BetaPowerSettings.scoredarousalsinwake = settings.scoredarousalsinwake; %1 for Lauren's scoring.
        % Note Epochs WakeSleep==8)=NaN; %unknown sleep is NaN
        
        %EKG analysis
        if isfield(settings, 'skipEKG')
            skipEKG=settings.skipEKG;
        else
            skipEKG=0;
        end
        if exist('EKG','var')&&(1-skipEKG)
            try
                displaytext=['EKG processing']; disp(displaytext); %                 
                if ~isfield(settings,'EKGplot'); settings.EKGplot=0; end
                [ECG_peak_i,HR] = EKGpeakdetection(EKGeeg,TimeEEG,1./Fs_EEG_,0,settings.EKGplot);
                ChannelsList = [ChannelsList,'HR'];
                ChannelsFs = [ChannelsFs;Fs_EEG_];
            catch
                disp('Failed EKG; check if signal is empty');
                ECG_peak_i=NaN;
            end
        else
            ECG_peak_i=NaN;
            %settings.processEEG=0; %overwrite
        end
        
        %     figure(9)
        %     plot(TimeEEG,EKGeeg)
        
        %settings.allpowers = 1;
        try
            disp('EEG beta power analysis');
            [WakeSleep,WakeSleepInfo,EEGsignals,EEGsignalsOut,Powersignals,PowersignalsOut,PowersignalsOutAll] = ...
                BetaPowerRun(EEGsignals,1/Fs_EEG_,BetaPowerSettings,ECG_peak_i);
            
            %Clean EEG signals
            for i=1:length(EEGsignals)
                eval([EEGsignals{i} '=EEGsignalsOut{i};'])
            end
            
            %uses info from Beta Sigmoid (phased out from here)
            try
                if 1
                    %best EEG
                    disp(WakeSleepInfo.AUC_M);
                    
                    [~,besti]=max(WakeSleepInfo.AUC_M);
                    EEG = eval([EEGsignals{besti}]);
                    ChannelsList = [ChannelsList,{'EEG'}];
                    ChannelsFs = [ChannelsFs;Fs_EEG_];
                    %add best EEG Powers
                    ChannelsList = [ChannelsList,'WakeSleep',Powersignals];
                    ChannelsFs = [ChannelsFs ;settings.Fs; settings.Fs*ones(length(Powersignals),1)];
                    %Bring in Powersignals
                    for i=1:length(Powersignals)
                        eval([Powersignals{i} '=PowersignalsOut{i};'])
                    end
                end
            catch AssignEEG_fail                
                AssignEEG_fail.message
            end
            
            if settings.allpowers
                Powersignalsalllist=[];
                for j=1:length(EEGsignals)
                    for i=1:length(PowersignalsOutAll{j})
                        eval([Powersignals{i} num2str(j) '=PowersignalsOutAll{j}{i};']);
                        Powersignalsalllist{length(Powersignalsalllist)+1} = [Powersignals{i} num2str(j)];
                    end
                end
                PowersignalsalllistFs = settings.Fs * ones(length(Powersignalsalllist),1);
                
                ChannelsList = [ChannelsList,Powersignalsalllist];
                ChannelsFs = [ChannelsFs;PowersignalsalllistFs];
                
            end
            
            
            clear EEGsignalsOut PowersignalsOut
            
        catch EEGbetaPowerAnalysis_fail
            disp(EEGbetaPowerAnalysis_fail.message);
            EEG = eval([EEGsignals{1}]);
            ChannelsList = [ChannelsList,{'EEG'}];
            ChannelsFs = [ChannelsFs;Fs_EEG_];
            WakeSleepInfo.AUC_M=NaN;
        end
        
    else %No EEG present: Empty list of EEG channels
        disp('Warning: No EEG');
        WakeSleepInfo.AUC_M=NaN;
        if exist('EKG','var')
            displaytext=['EKG processing'];
            disp(displaytext);
%             set(handletext,'String',displaytext); drawnow;
            try
                Fs_EKG = ChannelsFs(find(strcmp(ChannelsList,'EKG')==1));
                TimeTempEKG=(StartTime:(1/Fs_EKG):StartTime+(length(EKG)-1)*(1/Fs_EKG))'; % This is the time vector associated with the 100Hz Flow data.
                if ~isfield(settings,'EKGplot'); settings.EKGplot=0; end
                [ECG_peak_i,HR] = EKGpeakdetection(EKG,TimeTempEKG,1./Fs_EKG,0,settings.EKGplot);
                ChannelsList = [ChannelsList,'HR'];
                ChannelsFs = [ChannelsFs;Fs_EKG];
            catch me
                disp('failed EKG analysis');
            end
        end
    end
    
    clear TimeTempEKG    
    Info.WakeSleepInfo = WakeSleepInfo;
    
    %% Debug by plotting (need same Fs)
    if 0
        figure(1);
        %ax(1)=subplot(5,1,1); plot(Time,Flow);
        ax(2)=subplot(5,1,2); plot(Time,[Epochs EventsAr EventsResp]);
        ax(3)=subplot(5,1,3); plot(TimeEEG,EEG);
        ax(4)=subplot(5,1,4); plot(Time,EventsAr);
        ax(5)=subplot(5,1,5); plot(Time,WakeSleep);
        linkaxes(ax,'x');
    end
    
    clear TimeEEG
    
    %% Pes decontaminate
    try
        if 0 %&not leaving clean signal, needs work
            if exist('Pes','var')&&exist('EKG','var')
                disp('Pes EKG artifact removal')
                Fs_EKG = ChannelsFs(find(strcmp(ChannelsList,'EKG')==1));
                if Fs_EKG~=settings.Fs
                    'resampling EKG to match Resp'
                    temp = resample(EKG,settings.Fs,Fs_EKG);
                    deltalength = length(temp)-length(Time);
                    if deltalength>0
                        'removing samples in EKG to match Resp'
                        temp(end-deltalength+1:end)=[];
                    elseif deltalength<0
                        'adding samples in EKG to match Resp'
                        temp(end:end-deltalength)=temp(end);
                    end
                    Pes_clean = PesRemoveEKG(Pes,temp,1/settings.Fs,ECG_peak_i);
                else
                    Pes_clean = PesRemoveEKG(Pes,EKG,1/settings.Fs,ECG_peak_i);
                end
                Pes = Pes_clean;
            end
        end
    catch PesDecontaminate_fail
        disp('failed Pes decontaminate');
        disp(PesDecontaminate_fail.message);
    end
    
    %% Flow decontaminate
    try
        if isfield(settings,'FlowRemoveEKG') && settings.FlowRemoveEKG %&not leaving clean signal, needs work
            if exist('Flow','var')&&exist('EKG','var')
                disp('Flow EKG artifact removal')
                if Fs_EKG~=settings.Fs
                    'resampling EKG to match Resp'
                    temp = resample(EKG,settings.Fs,Fs_EKG);
                    deltalength = length(temp)-length(Time);
                    if deltalength>0
                        'removing samples in EKG to match Resp'
                        temp(end-deltalength+1:end)=[];
                    elseif deltalength<0
                        'adding samples in EKG to match Resp'
                        temp(end:end-deltalength)=temp(end);
                    end
                else
                    temp = EKG;
                end
                if settings.Fs~=Fs_EEG_
                    temp2 = round(ECG_peak_i*settings.Fs/Fs_EEG_); %check this works
                else
                    temp2 = ECG_peak_i;
                end
                PesRemoveEKG(Flow,temp,1/settings.Fs,temp2,[0.5 2.5],1,[],0.2);
                Flow_clean = PesRemoveEKG(Flow,temp,1/settings.Fs,temp2,[0.5 2.2],0,100,0.2);
                
                %             Pes = Pes_clean;
            end
        end
    catch FlowDecontaminate_fail
        disp('failed Flow decontaminate');
        disp(FlowDecontaminate_fail.message);
        
    end
    
    %% Remove from channel list
    ChannelsRemoveList={'REF1','REF2'};
    if 0
        ChannelsRemoveList={'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12','REF1','REF2'};
    end
    
    for i=length(ChannelsList):-1:1
        for j=length(ChannelsRemoveList):-1:1
            if strcmp(ChannelsRemoveList{j},ChannelsList{i})==1
                eval(['clear ' ChannelsList{i}]);
                ChannelsList(i)=[];
                ChannelsFs(i)=[];
            end
        end
    end
    
    %% Remove artifact in all signals using text files
    % UPDATE: _mask and _art
    % artefact -> remove periods that are complete rubbish and distracting
    % mask -> do not include these periods in analysis (moreover, summary analysis)
    %
    % this would allow re-analysis of data without re-running entire process
    %
    % also, write artifact as separate fn, with capacity to handle various
    % time formats (seconds, seconds+offset, time of day)
    %
    if 1
        for j=1:length(ChannelsList)
            textfilename=[directory fname(1:end-4) '_' ChannelsList{j} '_art.txt'];
            dt = 1/ChannelsFs(j);
            if exist(textfilename,'file')==2
                displaytext = ['Removing artifact from ' ChannelsList{j}];
                disp(displaytext); 
%                 set(handletext,'String',displaytext); drawnow;
                %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
                % textread was not recommended, upgraded to dlmread
                
                offsettime=StartTime;
                [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
                for i=1:length(col1)
                    %dttemp=eval([channels{j} '.interval']);
                    if isfield(settings,'DoNotUseOffsetTimeInArtefactRemoval') && settings.DoNotUseOffsetTimeInArtefactRemoval==1
                        lefti=round(col1(i)/dt+1);
                        righti=round(col2(i)/dt+1);
                    else
                        lefti=round((col1(i)-offsettime)/dt+1);
                        righti=round((col2(i)-offsettime)/dt+1);
                    end
                    if lefti<1, lefti=1; end
                    if righti>length(eval(ChannelsList{j})), righti=length(eval(ChannelsList{j})); end
                    eval([ChannelsList{j} '(lefti:righti)=NaN;']);
                end
                Percent_removed = sum(isnan(eval(ChannelsList{j})))/length(eval(ChannelsList{j}))*100;
                displaytext = [num2str(Percent_removed), '% of ', ChannelsList{j}, ' removed as artifact'];
                disp(displaytext); 
                % set(handletext,'String',displaytext); drawnow;
                % figure(); plot(Time,eval([ChannelsList{j}]));
            end
        end
    end
    
    
    %% SpO2 artifact (set to zeros, NaN makes future artifact detection fail)
    if 1
        try
            if exist('SpO2', 'var')
                I = SpO2>101|SpO2<40;
                Ileft = [0;I(1:end-1)];
                Iright = [I(2:end);0];
                SpO2(I==1)=0;
                SpO2(Ileft==1)=0;
                SpO2(Iright==1)=0;
                %         figure()
                %         plot(SpO2)
            else
                disp('SpO2 signal not detected - filling with NaN')
                SpO2 = nan(size(Flow));
                ChannelsList{end+1} = 'SpO2';
                ChannelsFs(end+1) = ChannelsFs(strcmp(ChannelsList,'Flow'));
            end
        catch SpO2artifact_fail
            disp('failed SpO2 artifact detection')
            SpO2artifact_fail.message
        end
    end
    
    %% Get Additional Signals from EDF file    
    if isfield(settings,'BringInNewSignal')
        
        %HarmonizedChannelSpreadsheet contains options for signal labels
        [~,~,HarmonizedChannelLabels] = xlsread(HarmonizedChannelSpreadsheet,1,'A1:X10000');
        LabelHeaders=HarmonizedChannelLabels(1,:); % headers
        HarmonizedChannelLabels(1,:)=[]; % label options
        
        if 0
            [LabelH,~,FsH] = EDFChannelLabels([directory '\' fname]); % from EDF
            LabelH=LabelH';
        else
            LabelH = signalHeader.signal_labels';
            FsH = signalHeader.samples_in_record';
        end
        
        for kk=1:length(settings.BringInNewSignal)
            SigtoLoad=settings.BringInNewSignal{kk};
            try
                %get all the possible options from Excel Sheet
                SigtoLoadLabelOptions=HarmonizedChannelLabels(:,find(strcmp(SigtoLoad,LabelHeaders)));
                
                % compare it with EDF label list and get channel number
                ChannelMatch=[];
                for ll=1:length(SigtoLoadLabelOptions)
                    if isempty(ChannelMatch)
                        ChannelMatch=find(strcmp(SigtoLoadLabelOptions(ll),strtrim(LabelH))); % gives position in LabelMatch
                    end
                end
                
                if ~isempty(ChannelMatch) % load the signal
                    
                    if isfield(settings,'BlockEDFload') && settings.BlockEDFload>0
                        [~,~,signalCell] = blockEdfLoad(EDFfilenamedir,signalHeader.signal_labels(ChannelMatch));
                        data = signalCell{1};
                        Fsdata=signalHeader.samples_in_record(ChannelMatch);
                        Labeldata=signalHeader.signal_labels{ChannelMatch};
                    else
                        [data,Fsdata,~,~,Labeldata,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelMatch-1,0,Inf); %the '-1' is because EDF channel numbers start from 0.
                    end
                    Fsdata = Fsdata*(1-settings.FsLeak);
                    
                    if strcmp(SigtoLoad,'Pulse')  % for pulse signal
                        data(data<20)=NaN;
                        data(data>200)=NaN;
                    end
                    
                    eval([SigtoLoad '=data;']);
                    displaytext=[SigtoLoad ':Found channel labelled:' Labeldata ' at ' num2str(Fsdata), ' Hz'];
                    disp(displaytext);
%                     set(handletext,'String',displaytext); drawnow;
                    ChannelsList{end+1} = SigtoLoad;
                    ChannelsFs(end+1) = Fsdata;
                    
                else
                    disp (['signal not found:' SigtoLoad])
                end
                
            catch
                disp([SigtoLoad ':loading signals failed'])
            end
        end
    end
    
    %%    
    if isfield(settings,'PPGisCalledPulse') && settings.PPGisCalledPulse==1
        try
            PPG=Pulse; %careful this will overwrite PPG if it exists.
            
            clear Pulse;
            I = find(ChannelsList=="Pulse");
            ChannelsList(I)=[];
            ChannelsFs(I)=[];
        end
    end
    
    %% Resample channels
    % SpO2_ = SpO2;
    % ChannelsFs = ChannelsFs;
    for i=2:length(ChannelsList) %skip "Time"
        if ChannelsFs(i)~=settings.Fs
            displaytext=['Resampling: ' ChannelsList{i} ' from ' num2str(ChannelsFs(i)) ' to ' num2str(settings.Fs) ' Hz'];
            disp(displaytext); 
%             set(handletext,'String',displaytext); drawnow;
            if strcmp(ChannelsList(i),'SpO2')||strcmp(ChannelsList(i),'Position') || strcmp(ChannelsList(i),'Pulse')  %use nearest rather than interp
                TimeTemp=(StartTime:(1/ChannelsFs(i)):StartTime+(length(eval(ChannelsList{i}))-1)*(1/ChannelsFs(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(ChannelsList{i}),Time,'nearest');
                eval([ChannelsList{i} '= Temp;']);
                tol(i)=NaN;
            elseif 1%(settings.Fs/ChannelsFs(i))>2
                TimeTemp=(StartTime:(1/ChannelsFs(i)):StartTime+(length(eval(ChannelsList{i}))-1)*(1/ChannelsFs(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(ChannelsList{i}),Time,'linear');
                eval([ChannelsList{i} '= Temp;']);
                tol(i)=NaN;
            else %unused
                tol(i)=0.00000001;
                while 1
                    [Q,P] = rat(ChannelsFs(i)/settings.Fs,tol(i));
                    if (Q*P)>2^31
                        tol(i)=tol(i)*10;
                    else
                        break
                    end
                end
                %eval([ChannelsList{i} ' = resample(' ChannelsList{i} ',round(settings.Fs),round(ChannelsFs(i)));']); %only works if Fsamp / Fs are integer multiples of each other
                eval([ChannelsList{i} ' = resample(' ChannelsList{i} ',P,Q);']); %only works if Fsamp / Fs are integer multiples of each other
            end
            ChannelsFs(i)=settings.Fs;
        end
    end
    
    %%  Make channels the same length (sometimes these are off by a few samples)
    for i=1:length(ChannelsList)
        if ~isempty(ChannelsList{i})
            tempX=abs(length(eval([ChannelsList{i} ]))-length(Time));
            if tempX>10
                disp(['Warning: channel length difference = ' num2str(tempX)]);
            end
            while length(eval([ChannelsList{i} ]))~=length(Time)
                if length(eval([ChannelsList{i} ]))<length(Time)
                    eval([ChannelsList{i} '(end+1)=' ChannelsList{i} '(end);']); %add a sample to the end if the channel is too short by 1
                elseif length(eval([ChannelsList{i} ]))>length(Time)
                    eval([ChannelsList{i} '(end)=[];']); %delete a sample from the end if the channel is too long by 1
                end
            end
        end
    end
    
    %% Channel data at XHz: Combining data into large matrix
    if 0
        displaytext = ['Setting ColumnHeads'];
        disp(displaytext);
%         set(handletext,'String',displaytext); drawnow;
        %These PSG signals are expected (hardcoded) in this order:
        %If signals are not available, we need to either (1) make a fake channel (above) or (2) edit the Analysis program to handle the absence
        ColumnHeadsList = {'Time','Thorax','Flow','Epochs','EventsResp','EventsAr','Abdomen','SpO2','EEG','Position'};
        
        ColumnHeads=NaN*zeros(1,length(ChannelsList));
        I=zeros(1,length(ChannelsList));
        for i=1:length(ColumnHeadsList)
            for j=1:length(ChannelsList)
                if strcmp(ColumnHeadsList{i},ChannelsList{j})==1
                    ColumnHeads(i)=j;
                    I(j)=1;
                end
            end
        end
        
        ExtraChannels = ChannelsList(I==0);
        ColumnHeadsList = [ColumnHeadsList ExtraChannels];
        
        ColumnHeads=NaN*zeros(1,length(ChannelsList));
        I=zeros(1,length(ChannelsList));
        for i=1:length(ColumnHeadsList)
            for j=1:length(ChannelsList)
                if strcmp(ColumnHeadsList{i},ChannelsList{j})==1
                    ColumnHeads(i)=j;
                    I(j)=1;
                end
            end
        end
    end
    
    %% Combining data into large matrix
    displaytext='Combining data into large matrix';
    disp(displaytext); 
%     set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat = zeros(length(Time),length(ChannelsList));
    for x=1:length(ChannelsList)
        eval(['DataEventHypnog_Mat(:,x)=' ChannelsList{x} ';']);
    end
    %ChannelsList(find(cellfun(@isempty,ChannelsList)))=[];
    
    if 0
        run PlotXHzData;
    end
    
    figure(99); close(99);
    figure(5); close(5);
    figure(4); close(4);
    fclose all;
    
    %% Rename DataEventHypnog_Mat as SigT
    SigT=array2table(DataEventHypnog_Mat);
    
    try
        if sum(find(cell2mat(cellfun(@(x)any(~isnan(x)),ChannelsList,'UniformOutput',false))))>1
            SigT.Properties.VariableNames = ChannelsList;
        end
    end
    
%     %% Combine snore into large matrix - obsolete
%     if isfield(settings,'AnalyzeSnore') && settings.AnalyzeSnore==1
%         displaytext='Combining SNORE data into large matrix';
%         disp(displaytext); set(handletext,'String',displaytext); drawnow;
%         
%         savenamesnore = [directory(1:end-7),'Converted\',fname(1:end-4),'_snore.mat'];
%         save(savenamesnore,'SnoreInterpStruct','SnoreStruct', 'SnoreChannelsList', 'SnoreChannels_Fs', '-v7.3')
%     end
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               ------  Additional analyses  -------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %% find stable breathing periods
    if isfield(settings,'findStableBreathing') && settings.findStableBreathing==1
        % Here, we are attempting to find continuous periods of stable breathing.
        disp('Detecting periods of stable breathing');
        try
            % This adds full length stableBB_mins to DataEventHypnogMat.
            % Showing the number of minutes in each block of stable breathing
            % as minutes (integer from 0-10). The summary table is also saved
            % in the Info struct.
            [StableBB_mins, StableBB_Table] = getStableBreathingPeriods(SigT,ChannelsList);
            SigT.StableBreathing=StableBB_mins;
            ChannelsList{end+1} = 'StableBreathing';
            ChannelsFs(end+1) = settings.Fs;
            Info.StableBB = StableBB_Table;
        catch FaultFindingStableBB
            disp('Notice: Failed detection of stable breathing');
            if settings.verbose
                disp(FaultFindingStableBB.message);
            end
        end
    end
    
    
    %% WS analysis    
    try
        if isfield(settings,'WSanalysis') && settings.WSanalysis==1
            
            if sum(strcmp(ChannelsList,'Epochs'))==0
                %No Epochs Signal, add one
                Epochs = 2*ones(size(SigT,1),1); % DM changed from 1 (N2) to 2 (N1)
                SigT.Epochs = Epochs;
                ChannelsList = [ChannelsList {'Epochs'}];
            end
            
            if sum(strcmp(ChannelsList,'EventsAr'))==0
                %No Epochs Signal, add one
                EventsAr = 0*ones(size(SigT,1),1);
                SigT.EventsAr = EventsAr;
                ChannelsList = [ChannelsList {'EventsAr'}];
            end
            
            %Remove signals if they already exist (good for debug, etc)
            WSOutputList =  {'WPr','ArPr','EventsArWS','WPrB','ArPrB','EventsArWSB'};
            strmatch = any(string(ChannelsList)==string(WSOutputList(:)),1);
            SigT(:,strmatch)=[];
            ChannelsList(strmatch)=[];

            disp('WS analysis');
            load mdlA
            load mdlNoise
            load RefTable
            load mdlAcc
            
            if ~isfield(settings,'WSArVersion')
                settings.WSArVersion=2; %default to version 2 from 1/7/2020 2:10pm.
            end
            
            if settings.WSArVersion==1
                load mdlAR
            elseif settings.WSArVersion==2
                load mdlARwsbalance
            end
            
            [WPr,ArPr,WPrB,ArPrB,WSinfo] = RunWS(SigT,ChannelsList,(1-settings.scoredarousalsinwake),mdlA,mdlNoise,RefTable,mdlAcc,mdlAR);
            minarlength=2.5;
            EventsArWS = RemoveShortSegments(ArPr>0.5,minarlength,1./settings.Fs,1); % downsampling first would lead to some speed up here.
            EventsArWSB = RemoveShortSegments(ArPrB>0.5,minarlength,1./settings.Fs,1);

            if 1
                try 
                    figure(33); clf(33);
                    ax1(1)=subplot(3,1,1);
                    plot(SigT(:,1),SigT.Epochs);
                    hold on;
                    plot(SigT(:,1),SigT.EventsAr,'g','linewidth',2);
                    plot(SigT(:,1),max([WPr ArPr]')','color',[0.5 0.5 0.5]);
                    %plot(DataEventHypnog_Mat(:,1),max([WPrB ArPrB]')');
                    plot(SigT(:,1),WPr,'color',[0.5 0.5 0.9]);
                    plot(SigT(:,1),EventsArWS,'r','linewidth',1);
                    %plot(DataEventHypnog_Mat(:,1),ArSig);
                    ax1(2)=subplot(3,1,2);
                    ch = find(strcmp(ChannelsList,'EEG1')==1);
                    ch2 = find(strcmp(ChannelsList,'EEG3')==1);
                    plot(SigT(:,1),[SigT.EEG1 SigT.EEG3+nanstd(SigT.EEG1)]);
                    ax1(3)=subplot(3,1,3);
                    plot(SigT(:,1),SigT.Flow(:,ch));
                    linkaxes(ax1,'x')
                    
                    diff(get(gca,'xlim'));
                    temp = get(gca,'xlim');
                    set(gca,'xlim',temp);
                catch plot_fail
                end
            end
            
            
            if 1
                %% AA arousal intensity
                if ~isfield(settings,'skiparousalintensity') || settings.skiparousalintensity==1 %default is skip
                    %do nothing
                else 
                    try
                        TimeEEG=(StartTime:(1/Fs_EEG_):StartTime+(length(eval(EEGsignals{1}))-1)*(1/Fs_EEG_))'; % This is the time vector associated with the _XHz Flow data.
                        
                        try
                            if settings.runarousalintensity(1)==1
                                [SigT.ArIntensityOrig,SigT.ArIntensity]=ArousalIntensityRun(EEGsignals,EventsAr,Epochs,Time,TimeEEG,WSinfo.Acc);
                                ChannelsList = [ChannelsList,'ArIntensityOrig','ArIntensity'];
                                ChannelsFs = [ChannelsFs ; settings.Fs ;settings.Fs ];
                                I = find(diff(EventsAr)==1)+1;
                                I(Epochs(I)>=4)=[];
                                SigT.ArIntensityOrig(I);
                                Info.ArIntensity.ArIntensityOrigMedian = nanmedian(SigT.ArIntensityOrig(I));
                                Info.ArIntensity.ArIntensityOrigN = sum(~isnan(SigT.ArIntensityOrig(I)));
                                
                            end
                        catch AI_run1
                            AI_run1.message
                        end
                        
                        try
                            if settings.runarousalintensity(2)==1
                                [SigT.ArIntensityOrigWS,SigT.ArIntensityWS]=ArousalIntensityRun(EEGsignals,EventsArWS,Epochs,Time,TimeEEG,WSinfo.Acc);
                                ChannelsList = [ChannelsList,'ArIntensityOrigWS','ArIntensityWS'];
                                ChannelsFs = [ChannelsFs ; settings.Fs ;settings.Fs ];
                                I = find(diff(EventsArWS)==1)+1;
                                I(Epochs(I)>=4)=[];
                                SigT.ArIntensityOrigWS(I);
                                Info.ArIntensity.ArIntensityOrigWSMedian = nanmedian(SigT.ArIntensityOrigWS(I));
                                Info.ArIntensity.ArIntensityOrigWSN = sum(~isnan(SigT.ArIntensityOrigWS(I)));
                            end
                        catch AI_run2
                            AI_run2.message
                        end
                        
                        try
                            if settings.runarousalintensity(3)==1
                                [SigT.ArIntensityOrigWSB,SigT.ArIntensityWSB]=ArousalIntensityRun(EEGsignals,EventsArWSB,Epochs,Time,TimeEEG,WSinfo.Acc);
                                ChannelsList = [ChannelsList,'ArIntensityOrigWSB','ArIntensityWSB'];
                                ChannelsFs = [ChannelsFs ; settings.Fs ;settings.Fs ];
                                I = find(diff(EventsArWSB)==1)+1;
                                I(Epochs(I)>=4)=[];
                                SigT.ArIntensityOrigWSB(I);
                                Info.ArIntensity.ArIntensityOrigWSBMedian = nanmedian(SigT.ArIntensityOrigWSB(I));
                                Info.ArIntensity.ArIntensityOrigWSBN = sum(~isnan(SigT.ArIntensityOrigWSB(I)));
                            end
                        catch AI_run3
                            AI_run3.message
                        end
                        
                    catch meAI
                        if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
                            keyboard
                            meAI.getReport
                        else
                            disp('Failed Arousal Intensity')
                        end
                    end
                end
            end
                        
            %% Add signals from WS analyses as Channels
            SigT = [SigT table(WPr,ArPr,EventsArWS,WPrB,ArPrB,EventsArWSB)];
            ChannelsList = [ChannelsList,'WPr','ArPr','EventsArWS','WPrB','ArPrB','EventsArWSB'];
            
            ChannelsFs = [ChannelsFs;settings.Fs;settings.Fs;settings.Fs;settings.Fs;settings.Fs;settings.Fs];
            Info.WSinfo = WSinfo;
            %add signals back to DataEventHypnog_Mat ChannelsList ChannelsFs
            
            % Try REMEOG
            if isfield(settings,'phasicREMfromEOG') && settings.phasicREMfromEOG==1
                try
                    disp('starting EOG');
                    [SigT,ChannelsList,ChannelsFs]=getPhasicREM(SigT,ChannelsList,ChannelsFs);
                catch me
                    disp('failed EOG');
                    disp(me.message);
                end
            end
            
            % Try WSanalysis autostaging for when there is no human hypnogram
            if isfield(settings,'useWSforHypnogram') && settings.useWSforHypnogram>=1
                %SigT.WArPrB = max([SigT.WPrB' ; SigT.ArPrB'])';
                EpochLength=Evts.Hypnogram_t(2)-Evts.Hypnogram_t(1);
                logit = @(p) log(p./(1-p));
                logitinverse = @(x) 1./(1+exp(-x));
                clear TempHypnogram HypnogramEst
                for i=1:length(Evts.Hypnogram_t)
                    li = round((Evts.Hypnogram_t(i)-SigT.Time(1))/dt + 1);
                    if i<length(Evts.Hypnogram_t)
                        ri = round((Evts.Hypnogram_t(i+1)-SigT.Time(1))/dt );
                    else
                        ri = li+round(EpochLength/dt)-1;
                    end
                    if ri>height(SigT), ri=height(SigT); end
                    temp = SigT.WPrB(li:ri);        %%%%% choice of signal for this purpose
                    HypnogramEst(i,1)=nanmedian(temp);
                    if mean(isnan(temp))>0.5
                        HypnogramEst(i,1)=NaN;
                    end
                    TempHypnogram(i,1)=(HypnogramEst(i)>0.5) *2 + 2;
                    if isnan(TempHypnogram(i,1))
                        TempHypnogram(i,1)=8;
                    end
                    Evts.HypnogramWS=TempHypnogram;
                    
                end
                if settings.useWSforHypnogram==2 %actually use it, to be tested
                    Evts.Hypnogram = Evts.HypnogramWS;
                    Epochs = interp1(Evts.Hypnogram_t+EpochLength/2,Evts.Hypnogram,Time,'nearest','extrap');
                    % any data before Evts.Hypnogram_t(1) == 8;
                    Epochs(Time<Evts.Hypnogram_t(1)) = 8;
                    % any data after Evts.Hypnogram_t(end)+30 == 8;
                    Epochs(Time>(Evts.Hypnogram_t(end)+EpochLength)) = 8;
                    SigT.Epochs=Epochs;
                end
            end
        end 
    catch WSfail
        if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
            keyboard
            WSfail.getReport
        else
            disp('failed WS analysis')
        end
    end
    
    
    %% For HGNS data, make a signal that denotes stim on and off
    if isfield(settings,'MakeStimAmpl') && settings.MakeStimAmpl
        % Make Stime Ampl data
        for j=1:length(w)
            if strcmp(w(j).name,'StimAmpl')
                StimAmpl=filehandle.StimAmpl;
            end
        end

        if exist('StimAmpl','var')
            StimAmpl.values(:,2) = []; % remove column of 0s

            % next check if stim amplitude signal was not started at the
            % beginning of the study (if not enter 0s)
            if StimAmpl.times(1)>0
                timesadd = (1:StimAmpl.times(1)-1)';
                stimsadd = zeros(size(timesadd));

                StimAmpl.times = [timesadd; StimAmpl.times];
                StimAmpl.values = [stimsadd; StimAmpl.values];
            end
            StimTime = 1:floor(Time(end))';
            StimVals = interp1(StimAmpl.times,StimAmpl.values,StimTime,'previous');

            if 1
                figure(71), set(gcf, 'Position', [550 528 570 259])
                plot(StimAmpl.times,StimAmpl.values,StimTime, StimVals,'r--'), hold on
            end
        else
            StimTime = 1:floor(Time(end))';
            StimVals = zeros(size(StimTime));
        end
        StimAmpl = interp1(StimTime,StimVals,Time,'previous');
        SigT.StimAmpl = StimAmpl;
    end
    
    % Load text file with indices indicated HGNS on and off
    textfilename=[directory fname(1:end-4) '_stim_on_off.txt'];
    clear StimON
    if exist(textfilename,'file')==2 % using this as a flag for running stim on/off routine, generally
        displaytext = 'Finding HGNS stimulation on and off';
        disp(displaytext); 
        
        % Get original Chin signal for now, can fix later
        for j=1:length(w)
            for mm = 1:length(channelnameoptions.Chin)
                ChinChanName = channelnameoptions.Chin{mm};
                if strcmp(w(j).name,ChinChanName)
                    ChinLong=filehandle.(ChinChanName);
                    dtChin = ChinLong.interval;
                    ChinLong = ChinLong.values;
                end
            end
        end
        
         %%
        [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
%         dt = 1/ChannelsFs(strcmp(ChannelsList,'Chin'));\\
        StimSelect = false(size(Time));
        for i=1:length(col1)
            % Make channel that indicates which regions were selected for
            % analysis
            dt = Time(2)-Time(1);
            lefti=round(col1(i)/dt+1);
            righti=round(col2(i)/dt+1);
            
            % make sure indices don't go out of bounds
            if lefti<1, lefti=1; end
            if righti>length(Time), righti=length(Time); end
            
            % make region = 1
            StimSelect(lefti:righti) = true;

            
            [StimONtemp, TimeStimON] = FindStimON(ChinLong,col1(i),col2(i),dtChin,StartTime,'on_off');
            if ~exist('StimON','var')
                StimON = StimONtemp;
            else
                StimON = StimON | StimONtemp;
            end
        end
        
        %% Adjust stim amplitiude based on notes     
        textfilename=[directory fname(1:end-4) '_stim_amplitudes.txt'];
        if exist(textfilename,'file')==2
            displaytext = ['Adjusting stimulation amplitudes'];
            disp(displaytext);
            %                 set(handletext,'String',displaytext); drawnow;
            %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
            % textread was not recommended, upgraded to dlmread
            offsettime = StartTime;
            dt=StimTime(2) - StimTime(1);
            [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2); col3 = M(:,3);
            for i=1:length(col1)
                %dttemp=eval([channels{j} '.interval']);
                if isfield(settings,'DoNotUseOffsetTimeInArtefactRemoval') && settings.DoNotUseOffsetTimeInArtefactRemoval==1
                    lefti=round(col1(i)/dt);
                    righti=round(col2(i)/dt);
                else
                    lefti=round((col1(i)-offsettime)/dt);
                    righti=round((col2(i)-offsettime)/dt);
                end
                if lefti<1, lefti=1; end
                if righti>length(StimVals), righti=length(StimVals); end
                StimVals(lefti:righti)=col3(i);
            end
            %                 set(handletext,'String',displaytext); drawnow;
            %figure(); plot(Time,eval([ChannelsList{j}]));
            figure(71)
            plot(StimTime, StimVals,'k--'), hold off

            StimAmpl = interp1(StimTime,StimVals,Time,'previous');
            SigT.StimAmpl = StimAmpl;
%             if any(strcmp(ChannelsList,'StimAmpl'))
%                 ChannelsFs(strcmp(ChannelsList,'StimAmpl')) = 1/dt;
%             else
%                 ChannelsList = [ChannelsList,'StimAmpl'];
%                 ChannelsFs(strcmp(ChannelsList,'StimAmpl')) = 1/dt;
%             end
            clear StimVals
        end
    end
    
    % NEED TO MAKE IT POSSIBLE TO MAKE THIS PART OF CODE NOT DEPENDENT ON
    % ABOVE. OR FIND A WAY TO RUN ABOVE IF M IS EMPTY
    textfilename2=[directory fname(1:end-4) '_hgns_on.txt'];
    if exist(textfilename2,'file')==2
        HGNSon = false(size(StimON));
        M_ = dlmread(textfilename2); col1_ = M_(:,1); col2_ = M_(:,2);
        for ii = 1:length(col1_)
            ONindices = TimeStimON >= col1_(ii)  & TimeStimON <= col2_(ii);
            HGNSon(ONindices) = true;
        end
    end
    
    if exist('StimON','var')
        StimONinterp = interp1(TimeStimON,double(StimON),Time,'Previous');
        StimON = StimONinterp;
        
        SigT.StimON = StimON;
        ChannelsList = [ChannelsList,'StimON'];
        FsStimON = round(1/(Time(2)-Time(1)));
        ChannelsFs = [ChannelsFs; FsStimON]; 
        
        SigT.StimSelect = StimSelect;
        ChannelsList = [ChannelsList,'StimSelect'];
        FsStimSelect = round(1/(Time(2)-Time(1)));
        ChannelsFs = [ChannelsFs; FsStimSelect]; 
    end
    
    if exist('HGNSon','var')
        HGNSonInterp = interp1(TimeStimON,double(HGNSon),Time,'Previous');
        HGNSon = HGNSonInterp;
        
        SigT.HGNSon = HGNSon;
        ChannelsList = [ChannelsList,'HGNSon'];
        FsHGNSon = round(1/(Time(2)-Time(1)));
        ChannelsFs = [ChannelsFs; FsHGNSon]; 
    end
    
    % Plot 
    if exist('StimONinterp','var')
        figure(43), set(gcf,'Position',[498 188 937 650])
        ax(1) = subplot(3,1,3);
        TimeChin = linspace(0,Time(end),length(Chin));
        WArPr = max([ArPr,WPr],[],2);
        
        plot(TimeChin,Chin,'Color',[0.5 0.5 0.5]), hold on;
        plot(Time,StimONinterp*3,'g',Time,HGNSonInterp*3,'y',Time,StimAmpl,'r',...
            Time, WArPr,'b',Time,StimSelect*3,'k'), hold off
        ylim([-5 10])
        ylabel('Stim info')

        ax(2) = subplot(3,1,1);
        plot(TimeChin,EEG3)
        ylabel('EEG')

        ax(3) = subplot(3,1,2);
        plot(TimeFlow,Flow)
        ylabel('Flow')
        xlabel('Time (s)')
        linkaxes(ax,'x');
        
        ax(1).XAxis.Exponent = 0;
        xtickformat('%0.5f')
    end
    % black = HGNSon, green = individual stims, red = Stim amplitude
    

    %% get CPAP data    
    if ~settings.ignoreCPAPdata
        [CPAPoff,CPAP] = getCPAP(SigT,settings,ChannelsList);
    else %a priori knowledge that CPAP is not administered
        CPAPoff=ones(height(SigT),1);
        CPAP=zeros(height(SigT),1);
    end
    
    
    %% get Position code
    try
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
        settings.supinepositioncode = settings.positioncodes(1);
    catch
        % above code failing if settings.protocol is char array.
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(n,:));
        settings.supinepositioncode = settings.positioncodes(1);
    end
    
    %% SpO2 analysis (uses CPAPoff info from above)
    try
        spo2Chan=find(strcmp(ChannelsList,'SpO2')==1);
        if ~isempty(spo2Chan) && sum(isnan(SigT.SpO2)) ~= size(SigT,1)
            [SpO2data,SpO2fixed]=SpO2analysis(SigT,ChannelsList,CPAPoff);
            if 1 %replace with artifact-rejected SpO2
                SigT.SpO2 = SpO2fixed;
                disp('Using fixed Spo2 (artifact-removed SpO2)')
                clear SpO2fixed;
            end
            Evts.SpO2 = SpO2data;
        else
            Evts.SpO2=NaN;
        end
    catch faultSpO2analysis
        Evts.SpO2=NaN;
        
        if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
            keyboard
        else
            if isfield(settings,'verbose') && settings.verbose==1
                disp(faultSpO2analysis.getReport);
            else
                disp(faultSpO2analysis.message);
            end
        end 
    end
    
    %% Get Arousal Table
    try
        ArChan=find(strcmp(SigT.Properties.VariableNames,'EventsAr')==1);
        if ~isempty(ArChan) && sum(SigT.EventsAr)>0
            %Regular Scoring
            [Evts.ArT,Evts.ArTinfo]=getArT(SigT,SigT.Properties.VariableNames,"EventsAr",Evts.Hypnogram,Evts.Hypnogram_t);
        else
            disp('Skipped Arousal Table, no arousals')
        end
    catch ArousalT_fail
        if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
            keyboard
        else
            if isfield(settings,'verbose') && settings.verbose==1
                disp(ArousalT_fail.getReport);
            else
                disp('failed ArT, possibly no scored arousals');
                disp(ArousalT_fail.message);
            end
        end
    end
    
    try
        if (isfield(settings,'WSanalysis') && settings.WSanalysis==1)
            %Autoscoring 1
            [Evts.EvtsArAuto.ArT,Evts.EvtsArAuto.ArTinfo]=getArT(SigT,SigT.Properties.VariableNames,"EventsArWS",Evts.Hypnogram,Evts.Hypnogram_t);
            try
                Evts.EvtsArAuto.RespT = Evts.RespT;
                ROImask=CPAPoff; %update
                [Evts.EvtsArAuto]=getAHIAll(SigT,ROImask,Evts.EvtsArAuto,SigT.Properties.VariableNames,settings,'EventsResp','EventsArWS');
            end
            
            
            %Autoscoring 2
            [Evts.EvtsArAutoB.ArT,Evts.EvtsArAutoB.ArTinfo]=getArT(SigT,SigT.Properties.VariableNames,"EventsArWSB",Evts.Hypnogram,Evts.Hypnogram_t);
            try
                Evts.EvtsArAutoB.RespT = Evts.RespT;
                ROImask=CPAPoff; %update
                [Evts.EvtsArAutoB]=getAHIAll(SigT,ROImask,Evts.EvtsArAutoB,SigT.Properties.VariableNames,settings,'EventsResp','EventsArWSB');
            end
        end
        
    catch AutoArT_fail
        if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
            keyboard
        else
            if isfield(settings,'verbose') && settings.verbose==1
                disp(AutoArT_fail.getReport);
            else
                disp('failed creating Autoscored Arousal Table');
                disp(AutoArT_fail.message);
            end
        end
    end
    
    % nanmedian(Evts.EvtsArAuto.ArT.WSBalanceMax3(Evts.EvtsArAuto.ArT.AASMarousal==1))
    % nanmedian(Evts.EvtsArAuto.ArT.ArBalanceMax3(Evts.EvtsArAuto.ArT.AASMarousal==1))
    % nanmedian(Evts.EvtsArAuto.ArT.ArInt(Evts.EvtsArAuto.ArT.AASMarousal==1))
    % nanmedian(Evts.EvtsArAuto.ArT.ArIntOr(Evts.EvtsArAuto.ArT.AASMarousal==1))
    
    %% AutoEditArousal
    if isfield(settings,'EditManuallyScoredArousalDurations')
        % if there are no manually scored Ar, then this will fail
        try
            [Evts,SigT]=ReplaceHumanArousalTimingWithAuto(Evts,SigT);
        catch AutoEditArousal_fail
            disp('failed to edit manually scored Ar durations using autoscored Ar');
        end
    end
    
    %% Get AHI, EventsRespT    
    try
        if  (~isfield(settings, 'SiteOfCollapse') || (isfield(settings, 'SiteOfCollapse') && settings.SiteOfCollapse == 0))
            % skip this if just doing DISE stuff.
            % let's setup a flag to skip AHI-type analyses -SS
            disp('Original Scoring:');
            [Evts]=getAHIAll(SigT,CPAPoff,Evts,SigT.Properties.VariableNames,settings);
        end
    catch faultAHI
        disp(faultAHI.getReport);
    end
    
    %% Get Epoch Table-- epoch level features
    try
        [Evts]=EpochsTable(SigT,Evts,AMasterSpreadsheet,settings,n);
    catch EpochsT_fail
        disp(EpochsT_fail.getReport);
    end
    
    %% Heart Rate Response per Azarbarzin
    try
        Evts=getEvtPRMainFn(Evts,SigT);
        disp('success: Heart rate and/or pulse rate added to Evts.RespT/Evts.ArT');
    catch
        disp('failed adding heart rate and/or pulse rate to Evts.RespT/Evts.ArT');
    end
    
   
    
    %% PSG Report   
    if isfield(settings,'PsgReport') && settings.PsgReport==1
        settings.onlykeepdataafterCPAPswitchedoff=1;
        fname=extractBefore(char(Filenames(n,1)),'.');
        Evts.PSGreportSummaryT=getPSGreport(SigT,Evts,ChannelsList,settings,fname,CPAPoff);
    end

    %% AHI Confidence Interval--MIKE P Algorithm
    % double check with getAHI3pa..not getting exact results...weird!!
    %  [AHI.ThreePA,Evts,~] = getDesatArSubset(Evts,Epochs,Time,Position,CPAPoff,desatlist(1));

    if isfield(settings,'AHIConfidenceIntervals') && settings.AHIConfidenceIntervals==1
        TSTAhiCI=8; % default 8- to set the tst for all patients as same.
        alpha_level=5; % 95% CI
        iterations=3000;
        [Evts.AHICIs.AHICI,Evts.AHICIs.CILow,Evts.AHICIs.CIHigh]=Bootstrap_CI_Auto(Evts,TSTAhiCI,alpha_level,iterations);
    end
    
    %% NSRR Saving Evts.RespT
    if isfield(settings,'NSRRSave') &&  settings.NSRRSave==1
         try
%         ID=string(extractBefore(fname,'-rs1-harm-base.edf'));
%         catch
           ID= settings.nsrrid;
        end
        Outfoldername=settings.napfolder;
        if ~(exist(Outfoldername, 'dir') == 7)
            mkdir(Outfoldername);
        end
        NsrrEvtsT=Evts.RespT;
        % using time in seconds from start time
        NsrrEvtsT.EventStart=NsrrEvtsT.EventStart-StartTimeNew;
        NsrrEvtsT.EventEnd=NsrrEvtsT.EventEnd-StartTimeNew;
        fnamensrr=[Outfoldername 'resp_eventstable_M.txt'];
        NsrrEvtsT.Properties.VariableNames=upper(NsrrEvtsT.Properties.VariableNames);
        writetable(NsrrEvtsT, fnamensrr,'WriteVariableNames',true,'Delimiter','\t')
        
        % Flow QC parameters-- trying to round to 2 decimal points
        ID=string(ID);
        try
            filter_lp_predict=Info.FlowQ.FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
            filter_lp_predict=round(filter_lp_predict*100)/100;
            filter_hp_predict=Info.FlowQ.FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
            filter_hp_predict=round(filter_hp_predict*10000)/10000;
            fraction_clipping=round(Info.FlowQ.FclippedTotal*10000)/10000;
            prob_upright=round(Info.FlowQ.PrUpright*100)/100;
            fraction_noise=round(Info.FlowQ.FnoiseAll*100)/100;
            snr_window=round(Info.FlowQ.SNRwindow*100)/100;
            fs_flow=Fs_Flow;
            FlowQcT = table(ID,fs_flow,prob_upright,fraction_noise,snr_window,filter_lp_predict,filter_hp_predict,fraction_clipping);
        catch
            FlowQcT=table(ID,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
            FlowQcT.Properties.VariableNames={'ID','fs_flow','prob_upright','fraction_noise','snr_window','filter_lp_predict','filter_hp_predict',...
                'fraction_clipping'};
        end
        fnamensrr=[Outfoldername 'resp_qc_flow.txt'];
        FlowQcT.Properties.VariableNames=upper(FlowQcT.Properties.VariableNames);
        writetable(FlowQcT,fnamensrr,'WriteVariableNames',true,'Delimiter','\t')
    end
 
    
catch GeneralImport_fail
    if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
        keyboard
    else
        if isfield(settings,'verbose') && settings.verbose==1
            disp(GeneralImport_fail.getReport);
        else
            disp(GeneralImport_fail.message);
        end
    end
%     set(handletext,'String',displaytext); drawnow;
end
