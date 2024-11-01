function ConvertToMat_DVSnore()
clc
diaryfilename = 'ConvertLog.txt';
convertfilename = 'AMasterSpreadsheet.xlsx';

%% test for specific machine
[~, str] = system('hostname'); % alternate cmd is getenv('COMPUTERNAME')
targetstr = 'GS528-3792';
if strncmp(str, targetstr, 10)
    close all
    clear global ExportDataSpreadsheet
    clear
    clc
    disp('Detected DLM uni computer - overwriting selected variables');
    diaryfilename = 'ConvertLog_DLM.txt';
    %convertfilename = 'ConvertDataSpreadsheet_OA.xlsx';
    convertfilename = 'ConvertDataSpreadsheet_PuPtest.xlsx';
end

%% Turn diary logging on
diary(diaryfilename);
diary on

%% Ismac
global particle
particle='\';
    if ismac()
        particle='/';
    end

%% global variables and settings
global AMasterSpreadsheet handletext F_samp settings
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = convertfilename; % appending _OA for oral appliance dataset;
end

%% add path to necessary inbuilt functions from older matlab versions
% mydir  = pwd;
%     addpath([mydir particle 'Builtin_MatlabFns_confirmed_necessary' particle]);
    %addpath([mydir(1:find(mydir=='\',1,'last')) 'Builtin_MatlabFns\confirmed_necessary\']);

%% Start processing
t_start = clock;
displaytext='Starting up the Convert to .mat process';
disp(displaytext); set(handletext,'String',displaytext); drawnow;

% read spreadsheet (Files worksheet)
[~,~,raw] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
disp('Warning: XLS cols BD, BE now point to LOC, ROC (eyes); BF=Every2ndEEGref; update xls if needed');

ChannelNumbers_ = raw(:,17:end-1);
ChannelNumbers = NaN*ones(size(ChannelNumbers_));
for i = 1:size(ChannelNumbers_,1)
    for j = 1:size(ChannelNumbers_,2)
        if ischar(ChannelNumbers_{i,j})
            ChannelNumbers(i,j)=NaN;
            %ChannelNumbers_{i,j}=NaN;
        else
            ChannelNumbers(i,j)=ChannelNumbers_{i,j};
        end
    end
end

ConvertMatFlag = cell2mat(raw(:,9));

%last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
lastrow = find(1*(~isnan(ConvertMatFlag)),1,'last');
    raw(lastrow+1:end,:)=[];
    ChannelNumbers(lastrow+1:end,:)=[];
    ConvertMatFlag(lastrow+1:end,:)=[];

Filenames = raw(:,2:8);

% read spreadsheet (Options worksheet)
[~,~,options] = xlsread(AMasterSpreadsheet,2,'C3:C11');
F_samp = options{1};
ExportDataDirectory = char(options(2)); 
if ExportDataDirectory(end)~=particle
    ExportDataDirectory=[ExportDataDirectory particle];
end
% SpikeTimePositionBy = options{3}; % I don't think this is used, remove from spreadsheet?
settings.scoredarousalsinwake = options{3};


settings.EverySecondEEGisRef = cell2mat(raw(:,end))==1;   %NaN (empty) taken as zeros, negative  
settings.RunInReverse = logical(options{5});
settings.EverySecondEKGisRef = logical(options{6}); 
settings.FsLeak = double(options{7}); 
settings.skiparousalintensity = options{8};
settings.AnalyzeSnore = options{9};

if ~settings.RunInReverse
    runrange=(1:size(Filenames,1));
else
    runrange=(size(Filenames,1):-1:1);
end

%% Loop
for n=runrange
    errorlist=[];
    try
        if ConvertMatFlag(n)
            %tic
            disp(' '); % add row space for visual clarity in command window
            displaytext=['Starting Convert to .mat for: ' Filenames{n}(1:end-4)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            [DataEventHypnog_Mat,ChannelsList,ChannelsFs,Evts,WakeSleepInfo] = GeneralImport(n,Filenames,ChannelNumbers);
            try
                fname=[Filenames{n}(1:end-4) '_XHz' '.mat'];
                displaytext=['Saving data to: ' fname];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                save([ExportDataDirectory, fname],'DataEventHypnog_Mat','ChannelsList','ChannelsFs','Evts','WakeSleepInfo','-v7.3');
                displaytext='Finished saving data';
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
            catch me
                disp(me.message); 
                % disp(getReport(me));
            end
            %toc
        else
            fname=[Filenames{n}(1:end-4) '_XHz' '.mat'];
            displaytext=['User skipped: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
        end
    catch me
        %me.message
        fname=[Filenames{n}(1:end-4) '_XHz' '.mat'];
        displaytext=['Error converting: ' fname];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        errorlist = [errorlist n];
    end
end

delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['Convert process complete. Total time: ', char(D), ' (hh:mm:ss)']; 
disp(displaytext); set(handletext,'String',displaytext); drawnow;

if size(errorlist)>0
    displaytext=['Errors: ' num2str(errorlist)];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

diary off
%/end ConvertToMat



%% GeneralImport: Currently works for Minerva, Profusion, Spike
function [DataEventHypnog_Mat,ChannelsList,ChannelsFs,Evts,WakeSleepInfo] = GeneralImport(n,Filenames,ChannelNumbers)
global handletext F_samp settings particle

system = char(Filenames(n,7));

try
    directory = char(Filenames(n,4));
%     % option to append \ for consistency, but would need to change some code below
%     if directory(end)~='\'
%         directory=[directory '\'];
%     end
    fname=char(Filenames(n,1));
   
%     if strcmp(system,'ProfusionXML')|strcmp(system,'Minerva')|strcmp(system,'ProfusionXML_MrOs')
%         SignalFormat='EDF';
%     elseif strcmp(system,'Spike')
%         SignalFormat='Spike';
%     end
    
    switch system
        case {'ProfusionXML','Minerva','ProfusionXML_MrOs','BrainRT','RemLogic','RemLogicTextFr','Alice','Deltamed'}
            SignalFormat='EDF';
            if isnan(ChannelNumbers(n,1))
                disp('Flow is channel number NaN, skipping')
                WakeSleepInfo=NaN;
                ColumnHeads=NaN;
                DataEventHypnog_Mat=NaN;
                ChannelsList=NaN;
                ColumnHeadsList=NaN;
                return
            end
        case 'Spike'
            SignalFormat='Spike';
        case 'NoxMat'    
            SignalFormat='NoxMat';
    end
    
    %% Signals
    %This section imports respiratory flow and other polysomnography signal data from the EDF files
    displaytext='Loading signals';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
    %Default Signals List
    ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12','Position','Pmask','EKG','EKG2','LOC','ROC'};
    %SS: first channel renamed to Flow from Pnasal, otherwise later the function will not function
    
    if strcmp(SignalFormat,'EDF')
        %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
        Channels_Fs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        EDFfilenamedir = [directory fname];
        
        for i=1:length(ChannelsList)
            displaytext=['Collecting channel:' ChannelsList{i}];
            try
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                
                eval(['[' ChannelsList{i} ',Channels_Fs(i),~,~,LabelTemp,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
                %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory particle fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
                displaytext=['Found channel labelled:' LabelTemp ' at ' num2str(Channels_Fs(i)), ' Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                ChannelFound(i)=1;
            catch me
                displaytext=['No channel:' ChannelsList{i} ' ' me.message];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                %eval([ChannelsList{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
            end
        end
    elseif strcmp(SignalFormat,'Spike')
        %%
        %Overwrite default ChannelsList
        ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','LOC','ROC',...
            'Position','Pmask','EKG','alphaFlow','kFlow','Pes','Edi','Epi','GGpmax','Pnasal'...
            };
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        MultiplyPositionChannelby=5; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat','Sao2ear','Sao2fing','SaO2fing'};
        channelnameoptions.Position={'Position','position','Pos','pos','NoxPos'};
        %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal', 'CO2'};
        %channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
        channelnameoptions.EEG1={'EEG1','EEG_C3_A2_clean','EEG_C3_A2','C3_A2','C3-A2'};
        channelnameoptions.EEG2={'EEG2','EEG_C4_A1_clean','EEG_C4_A1','C4_A1','C4-A1'};
        channelnameoptions.EEG3={'EEG3','EEG_F3_A2_clean','EEG_F3_A2','F3_A2','F3-A2'};
        channelnameoptions.EEG4={'EEG4','EEG_O2_A1_clean','EEG_O2_A1','O2_A1','O2-A1'};
        channelnameoptions.LOC={'LOC'};
        channelnameoptions.ROC={'ROC'};
        channelnameoptions.EKG={'EKG','ECG'};
        channelnameoptions.Thorax={'ThNoxRIP','Thorax','RC','Chest','CHEST','Belt2'};
        channelnameoptions.Abdomen={'AbNoxRIP','Abdomen','ABD','Abdom','ABDM','Belt1'};
        channelnameoptions.alphaFlow={'alphaFlow'};
        channelnameoptions.kFlow={'kFlow'};
        channelnameoptions.Pes={'Pes','pes','Pes_clean'};
        channelnameoptions.Edi={'Edi','edi'}; %'EMGdi'
        channelnameoptions.Epi={'Epi','epi','Pepi','pepi'};
        %channelnameoptions.FlowEdi={'FlowEdi'};
        %channelnameoptions.FlowPes={'FlowPes'};
        channelnameoptions.GGpmax={'GGpmax','GGPmax'};
        channelnameoptions.Flow={'Vflow','Flow','flow','Pnasal','PNasal','Pmask','PMask'}; %'Vflow'
        channelnameoptions.Pnasal={'Pnasal','PNasal'}; %
        channelnameoptions.Pmask={'CPAP','Pmask','PMask'};
        channelnameoptions.Snore={'Snore','snore'};
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
        
        Channels_Fs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);

        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    Channels_Fs(i)=1./eval([ChannelsList{i} '.interval']);
                    
                    % this IF was added by DLM. For some reason, alphaFlow
                    % channel in 1161 had an interval of 0.0013, however 
                    % the data length was the same as other channels, 
                    % implying that the interval should have been the same 
                    % as the other channels (i.e. 0.008).
                    if 0 % this is not required with new spike to matlab 
                         % exports, and should be removed now.
                    if (Channels_Fs(i) > 740) && (Channels_Fs(i) < 741)
                        str = ['---- CAUTION ----   Manual set Fs for ', ChannelsList{i}]; disp(str);
                        Channels_Fs(i) = 125;
                        pause;
                    end
                    end   
                catch me
                    Channels_Fs(i)=F_samp; %assumed
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
        %%
        %Overwrite default ChannelsList
        ChannelsList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
            'X_Axis','Y_Axis','Z_Axis','Pmask','EKG','Pes','TcPCO2','AudioVolume','Pleth'...
            };
        filehandle = matfile([directory fname]);
        w = whos('-file',[directory fname]);
        MultiplyPositionChannelby=1; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
        
        clear channelnameoptions
        channelnameoptions.SpO2={'SpO2'};
        channelnameoptions.X_Axis={'X_Axis'};
        channelnameoptions.Y_Axis={'Y_Axis'};
        channelnameoptions.Z_Axis={'Z_Axis'};
        %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal'};
        %channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
        
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
        
        channelnameoptions.EKG={'EKG','ECG','ecg'};
        channelnameoptions.Thorax={'Inductance_Thorax'};
        channelnameoptions.Abdomen={'Inductance_Abdomen'};
        
        channelnameoptions.Pes={'Millar_Pes','channel_3','Channel_3','Channel_5'};

        channelnameoptions.Flow={'Channel_7','Pneumotach_Flow'}; %'Vflow'

        channelnameoptions.Pmask={'Mask_Pressure'};
        
        channelnameoptions.TcPCO2={'TcCO2'};
        channelnameoptions.AudioVolume={'Audio_Volume'};
        channelnameoptions.Pleth={'Pleth'};
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
        
        Channels_Fs = NaN*ones(length(ChannelsList),1);
        ChannelFound = zeros(length(ChannelsList),1);
        %ChannelStartTime = NaN*ones(length(ChannelsList),1);
        Signals = filehandle.Signals;
        
        TimeInfo = filehandle.TimeInfo;
        
        %% - 0.0165 
        
            
        if 0
            deltaC1 = 0%-0.0165;%-0.0165 Penn; -0.0216/-0.0235 Iceland;
            deltaC1sperhr = deltaC1/250*3600
            
        deltaA1 = 0.0132;
            X=4
        figure(1); %clf(1);
        
        if 0
        plota = EEG1;
        plota.data(length(EEG2.data)+1:end) = [];
        plota.data = plota.data - EEG2.data
        hold('on')
        plota.fsTrue = plota.fsTrue + deltaA1;
        temptime = (plota.StartTime:1./plota.fsTrue:plota.StartTime+(length(plota.data)-1)*1./plota.fsTrue)';
        ax1(1)=subplot(X,1,1); plot(temptime,plota.data);
        end
       
        if 0
        plotb = Pes;
        plotb.fsTrue = plotb.fsTrue + deltaC1;
        temptime = (plotb.StartTime:1./plotb.fsTrue:plotb.StartTime+(length(plotb.data)-1)*1./plotb.fsTrue)';
        ax1(2)=subplot(X,1,2); plot(temptime,plotb.data);
        linkaxes(ax1,'x');
        hold('on')
        end
        
        if 1
        plotb = Flow;
        plotb.fsTrue = plotb.fsTrue + deltaC1;
        %plotb.fsTrue = 250.02547
        temptime = (plotb.StartTime:1./plotb.fsTrue:plotb.StartTime+(length(plotb.data)-1)*1./plotb.fsTrue)';
        ax1(3)=subplot(X,1,3); plot(temptime,plotb.data);
        linkaxes(ax1,'x');
        end
        
        if 1
            hold('on')
        plotb = Pmask;
        plotb.fsTrue = plotb.fsTrue + deltaA1;
        %plotb.fsTrue = 200.008127
        temptime = (plotb.StartTime:1./plotb.fsTrue:plotb.StartTime+(length(plotb.data)-1)*1./plotb.fsTrue)';
        ax1(3)=subplot(X,1,3); plot(temptime,-plotb.data,'k');
        linkaxes(ax1,'x');
        end
        
        if 0
            hold('on')
        plotb = Abdomen;
        plotb.fsTrue = plotb.fsTrue + deltaA1;
        temptime = (plotb.StartTime:1./plotb.fsTrue:plotb.StartTime+(length(plotb.data)-1)*1./plotb.fsTrue)';
        ax1(3)=subplot(X,1,4); plot(temptime,plotb.data);
        linkaxes(ax1,'x');
        end
        
        if n==10
            xlim([89790 89880]);
            %pause
            %xlim([81000 81100]);
        end
        %xlim([103970    103980]);
%to zoom in and see asynchrony:
%xlim([91048 91064]);
%signals are more in sync earlier in night, suggesting fs not start times are the problem
%xlim([81048 81064]); 
        end
        
        %
        if 0
        if exist('Pes')
            Pes.fsTrue = Pes.fsTrue + deltaC1; 
        end
        if exist('Flow')
            Flow.fsTrue = Flow.fsTrue + deltaC1; 
        end   
        end
        %%
%         clear tempnewlabelsNoSpace
%         for i=1:length(Signals.labelsNoSpace)
%             temp = Signals.labelsNoSpace(i,:);
%             itemp = find(temp==' ');
%             temp(itemp)=[];
%             tempnewlabelsNoSpace{i,1}=temp;
%         end
%         Signals.labelsNoSpace = tempnewlabelsNoSpace;
        
%get sampling rates
        for i=1:length(ChannelsList)
            if exist(ChannelsList{i},'var')
                ChannelFound(i)=1;
                try
                    %I=find(strcmp(channelnamestemp,ChannelsList{i}));
                    %I2=find(strcmp(Signals.labelsNoSpace,ChannelTitlesTemp{I}));
                    %Channels_Fs(i)=eval([ChannelsList{i} '.fs']);
                    Channels_Fs(i)=eval([ChannelsList{i} '.fsTrue']);
                catch me
                    Channels_Fs(i)=F_samp; %assumed
                    disp(['missing sampling rate assuming F_samp']);
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
        %StartTime = min(ChannelStartTime);
        %correct start times, not all signals start simultaneously
%         ChannelStartTimeDelay = ChannelStartTime-min(ChannelStartTime);
%         ChannelStartTimeDelayi = round(ChannelStartTimeDelay.*Channels_Fs);
%         ChannelStartTime2 = ChannelStartTime-ChannelStartTimeDelayi./Channels_Fs;
%         for i=1:length(ChannelsList)
%             if ChannelFound(i)==1
%                 FirstValue(i) = eval([ChannelsList{i} '(1)']); %assumes first value is not NaN
%                 temp = FirstValue(i)*ones(ChannelStartTimeDelayi(i),1);
%                 eval([ChannelsList{i} '=[temp; ' ChannelsList{i} '];']); %e.g. Flow=[temp; Flow];
%             else
%                 FirstValue(i) = NaN;
%             end 
%         end

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
        
        I = find(strcmp(ChannelsList,'X_Axis')==1)
        %Fs_temp = Channels_Fs(I);
        ChannelFound(I+2)=0; %assumes X,Y,Zaxis are together
        ChannelsList{I}='Position';
        ChannelsList{I+1}='Incline';
        
    end %signal import
    
    %%
    ChannelsList(ChannelFound==0)=[];
    Channels_Fs(ChannelFound==0)=[];
    
    
    %%
    %all EEGs must be same length
    %Make channels the same length (sometimes these are off by a few samples)
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
            clear EEG2;
            Itemp=find(strcmp(ChannelsList,'EEG2')==1);
            ChannelsList(Itemp)=[];
            Channels_Fs(Itemp)=[];
        end
        if exist('EEG3','var')&&exist('EEG4','var')
            EEG3=EEG3-EEG4;
            clear EEG4;
            Itemp=find(strcmp(ChannelsList,'EEG4')==1);
            ChannelsList(Itemp)=[];
            Channels_Fs(Itemp)=[];
        end
        if exist('EEG5','var')&&exist('EEG6','var')
            EEG5=EEG5-EEG6;
            clear EEG6;
            Itemp=find(strcmp(ChannelsList,'EEG6')==1);
            ChannelsList(Itemp)=[];
            Channels_Fs(Itemp)=[];
        end
        if exist('EEG7','var')&&exist('EEG8','var')
            EEG7=EEG7-EEG8;
            clear EEG8;
            Itemp=find(strcmp(ChannelsList,'EEG8')==1);
            ChannelsList(Itemp)=[];
            Channels_Fs(Itemp)=[];
        end
        if exist('EEG9','var')&&exist('EEG10','var')
            EEG9=EEG9-EEG10;
            clear EEG10;
            Itemp=find(strcmp(ChannelsList,'EEG10')==1);
            ChannelsList(Itemp)=[];
            Channels_Fs(Itemp)=[];
        end
        if exist('EEG11','var')&&exist('EEG12','var')
            EEG11=EEG11-EEG12;
            clear EEG12;
            Itemp=find(strcmp(ChannelsList,'EEG12')==1);
            ChannelsList(Itemp)=[];
            Channels_Fs(Itemp)=[];
        end
    end
    
    if settings.EverySecondEKGisRef
        if exist('EKG')&&exist('EKG2')
            EKG=EKG-EKG2;
            clear EKG2;
            Itemp=find(strcmp(ChannelsList,'EKG2')==1);
                ChannelsList(Itemp)=[];
                Channels_Fs(Itemp)=[];
        end
    end
    
    %base the study duration on Flow signal (presumably present)
    N_flow=length(Flow);
    Fs_Flow=Channels_Fs(find(strcmp(ChannelsList,'Flow')==1));
    N_timeXHz = round((N_flow/Fs_Flow*F_samp));
    
    displaytext=['Get recording start time'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    if strcmp(SignalFormat,'EDF')
        fid = fopen([directory fname],'r');
        fseek(fid,176,-1);
        StartTimeText = char(fread(fid,8,'char')');
        fclose(fid); % Close file
        StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
        if StartTime<43200, StartTime=StartTime+86400; end
    elseif strcmp(SignalFormat,'Spike')
        
        StartTime=0;
        try
           eval(['StartTime=filehandle.StarttimeSpike']);
        catch me
           %add CED function here in try-catch
           %if StartTime<43200, StartTime=StartTime+86400; end
           %
           % Will need to load the matson library to use the CED functions
        end
    elseif strcmp(SignalFormat,'NoxMat')
        StartTime = TimeInfo.StartTime;
    end
    
    %Time array
    Time=(StartTime:(1/F_samp):StartTime+(N_timeXHz-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    ChannelsList = ['Time',ChannelsList];
    Channels_Fs = [F_samp;Channels_Fs];
    
  
    %% Epochs and events
    %default start time for scoring is start time of recording e.g. Time(1)
    starttimediff=0;
    
    switch system
        case 'ProfusionXML_MrOs' 
            %ProfusionXml Hypnogram
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
            
            S=xml2struct([directory fname]);
            
            %Xml events
            displaytext=['Get Events data:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            if ~isfield(S.PSGAnnotation.ScoredEvents,'ScoredEvent')
                displaytext=['Warning: No scored events'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
            end
            
            Nevents = length(S.PSGAnnotation.ScoredEvents.ScoredEvent);
            
            clear EventCodesList EventStart EventDuration EventCodesCategory
            for i=1:Nevents
                EventCodesList{i} = S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.EventConcept.Text;
                EventCodesCategory{i} = S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.EventType.Text;
                EventStart(i) = StartTime + str2num(S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.Start.Text);
                EventDuration(i) = str2num(S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.Duration.Text);
            end %W=0,N1=1,N2=2,N3=3,R=5
            EventCodesListUnique = unique(EventCodesList);
            EventCodesCategoryUnique = unique(EventCodesCategory);
            EventCategoriesLabels= {...
                'ASDA arousal|Arousal (ASDA)',...
                'ASDA',...
                'Central apnea|Central Apnea',...
                'Mixed apnea|Mixed Apnea',...
                'Obstructive apnea|Obstructive Apnea',...
                'Hypopnea|Hypopnea',...
                'Unsure|Unsure',...
                };
            EventCategoriesCodes=[1 1 3 5 2 4 4];
            EventExact1orStartsWith2 = [1 2 1 1 1 1 1];
            
            displaytext=['Get Hypnogram data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            IsStage=NaN*zeros(Nevents,1);
            for i=1:Nevents
                IsStage(i) = strcmp(EventCodesCategory{i},'Stages|Stages');
            end %W=0,N1=1,N2=2,N3=3,R=5
            NepochChanges = sum(IsStage);
            EpochCodes=NaN*zeros(NepochChanges,1);
            Ilist=find(IsStage==1);
            for i=1:NepochChanges
                EpochCodes(i) = str2num(EventCodesList{Ilist(i)}(end));
            end
            EpochDuration=EventDuration(Ilist)/30;
            Nepochs=sum(EpochDuration);
            Hypnogram=NaN*zeros(Nepochs,1);
            iEpochEnd=cumsum(EpochDuration);
            iEpochStart=[0 iEpochEnd(1:end-1)]+1;
            for i=1:length(EpochDuration)
                Hypnogram(iEpochStart(i):iEpochEnd(i))=EpochCodes(i);
            end
            %Convert to Terrill code
            Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
            Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
            Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
            
            %end ProfusionXML_MrOs version
            
        case 'ProfusionXML' %code is untested
            %elseif strcmp('ProfusionXML',system) %code is untested
            %ProfusionXml Hypnogram
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
            
            S=xml2struct([directory fname]);
            
            %Xml events
            displaytext=['Get Events data:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            if ~isfield(S.CMPStudyConfig.ScoredEvents,'ScoredEvent')
                displaytext=['Warning: No scored events'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
            end
            
            Nevents = length(S.CMPStudyConfig.ScoredEvents.ScoredEvent);
            
            clear EventCodesList EventStart EventDuration
            for i=1:Nevents
                EventCodesList{i} = S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Name.Text;
                EventStart(i) = StartTime + str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Start.Text);
                EventDuration(i) = str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Duration.Text);
            end %W=0,N1=1,N2=2,N3=3,R=5
            EventCodesListUnique = unique(EventCodesList);
            EventCategoriesLabels= {...
                'Arousal',...
                'Central Apnea',...
                'Mixed Apnea',...
                'Obstructive Apnea',...
                'Hypopnea',... % unconvincing hypop (for MESA, MrOS)
                'Unsure',... % Special case for large cohort studies with >50% reduction in flow (e.g. MESA)
                };
            EventCategoriesCodes=[1 3 5 2 4 4];
            EventExact1orStartsWith2 = [2 1 1 1 1 1];
            
            
            displaytext=['Get Hypnogram data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            Nepochs = length(S.CMPStudyConfig.SleepStages.SleepStage);
            Hypnogram=NaN*zeros(Nepochs,1);
            for i=1:Nepochs
                Hypnogram(i) = str2num(S.CMPStudyConfig.SleepStages.SleepStage{1,i}.Text);
            end %W=0,N1=1,N2=2,N3=3,R=5
            %Convert to Terrill code
            Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
            Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
            Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
            
            %end ProfusionXML
            
        case 'BrainRT' %code is untested
            
            %elseif strcmp('ProfusionXML',system) %code is untested
            %ProfusionXml Hypnogram
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
            
            %Xml events
            displaytext=['Get Events data:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            tic
            S=xml2struct([directory fname]); %slow
            toc
            
            if ~isfield(S.BrainRTResult.ResultElements.EventsResultElement.Events,'Event')
                displaytext=['Warning: No scored events'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
            end
            
            %Keep only these sections of the xml struct, and clear S to save memory
            Sevents = S.BrainRTResult.ResultElements.EventsResultElement.Events.Event;
            Shypevents = S.BrainRTResult.ResultElements.HypnogramResultElement.ManualScoring.Events.Event;
            SLightsout = S.BrainRTResult.ResultElements.HypnogramResultElement.LightsOutEvents.Events.Event;
            
            clear S
            
            %make cell array table of events
            Nevents = length(Sevents);
            
            EventCellArray={''};
            formatIn = 'HH:MM:SS.FFF';
            clear EventCodesList EventStart EventDuration
            for i=1:Nevents
                S1 = Sevents{i};
                EventCellArray{i,1}=S1.Type.Text;
                EventCellArray{i,2}=S1.SubType.Text;
                EventCellArray{i,3}=S1.Start.Text;
                EventCellArray{i,4}=S1.End.Text;
                try
                    EventCellArray{i,5}=S1.CH1.Text;
                    EventCellArray{i,6}=S1.CH2.Text;
                catch me
                end
                EventCellArray{i,7}=S1.Validated.Text;
                
                EventCodesList{i} = [EventCellArray{i,1} '_' EventCellArray{i,2}];
                %duration
                temp = EventCellArray{i,3}(12:end-3);
                temp = mod(datenum(temp,formatIn),1)*86400;
                if temp<86400/2,temp=temp+86400; end
                EventStart(i)=temp;
                temp2 = EventCellArray{i,4}(12:end-3);
                temp2 = mod(datenum(temp2,formatIn),1)*86400;
                if temp2<86400/2,temp2=temp2+86400; end
                EventDuration(i) = temp2-temp;
                EventCellArray{i,8}=num2str(temp2-temp);
            end
            %
            %Possible event channel options, debug only:
            UniqueEventChannels = unique(EventCellArray(cellfun(@isempty,EventCellArray(:,5))==0,5));

            EventCodesListUnique = unique(EventCodesList);
            EventCategoriesLabels= {...
%                 '133_6',... %arousal
%                 '133_2',... %arousal maybe
%                 '133_3',... %arousal maybe
                '135_1',... %arousal maybe
                '129_3',... %event OA
                '129_2',... %event CA
                '129_4',... %event MA
                '129_6',... %event HypC
                '129_7',... %event HypO
                '129_5',... %event Hyp, assumed obstructive
                };
            EventCategoriesCodes=[1 2 3 5 6 4 4];
            EventExact1orStartsWith2 = [1 1 1 1 1 1 1];
            
            
            %Xml hypnogram (manual)        
            %directory = char(Filenames(n,6));
            %fname = char(Filenames(n,3));
            
            displaytext=['Get Hypnogram data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;

            %Shyp=xml2struct([directory particle fname]); %slow %only if need second separate xml for sleep staging. not needed now.
            
            Nhypevents = length(Shypevents);
            
            HypCellArray={''};
            clear HypCodesList HypStart HypDuration
            for i=1:Nhypevents
                S1 = Shypevents{i};
                HypCellArray{i,1}=S1.Type.Text;
                HypCellArray{i,2}=S1.SubType.Text;
                HypCellArray{i,3}=S1.Start.Text;
                HypCellArray{i,4}=S1.End.Text;
                try
                    HypCellArray{i,5}=S1.CH1.Text;
                    HypCellArray{i,6}=S1.CH2.Text;
                catch me
                end
                
                EventCodesList{i} = [HypCellArray{i,1} '_' HypCellArray{i,2}];
                %duration
                temp = HypCellArray{i,3}(12:end-3);
                temp = mod(datenum(temp,formatIn),1)*86400;
                if temp<86400/2,temp=temp+86400; end
                HypStart(i,1)=temp;
                temp2 = HypCellArray{i,4}(12:end-3);
                temp2 = mod(datenum(temp2,formatIn),1)*86400;
                if temp2<86400/2,temp2=temp2+86400; end
                HypDuration(i,1) = temp2-temp;
                HypCellArray{i,7}=num2str(temp2-temp);
            end
            
            clear IsStage
            for i=1:Nhypevents
                IsStage(i,1) = strcmp(HypCellArray{i,1},'128');
            end %W=0,N1=1,N2=2,N3=3,R=5
            
            NepochChanges = sum(IsStage);
            EpochCodes=NaN*zeros(NepochChanges,1);
            Ilist=find(IsStage==1);
            for i=1:NepochChanges
                EpochCodes(i) = str2num(HypCellArray{Ilist(i),2});
            end
            iEpochDuration=(round(HypDuration(Ilist)/30));
            
            % check time difference between lgihts on and end of EDF
            EndTime = mod(datenum(datestr(Time(end)/86400,'HH.MM.SS'),'HH.SS.MM'),1);
            EndTimeLOn = mod(datenum(SLightsout.End.Text(12:end-3),formatIn),1);
            % if some substantial difference bewteen EDF end time and lights on (i.e. 10 mins)
            % generate a new time vector relative to lights on
            if (EndTimeLOn - EndTime)*86400 > 600 
                EndTime = mod(Time(end)/86400,1);
                StartTimeNew = StartTime/86400 + EndTimeLOn-EndTime;
                Time = (StartTimeNew*86400:(1/F_samp):StartTimeNew*86400+((N_timeXHz-1)*(1/F_samp)))';
            end

%             % Change Time to lights off moment
%             if 0 %maybe all you need, but how to know automaticlaly when to do this?
%                 temp_LightsOut = mod(datenum(SLightsout.Start.Text(12:end-3),formatIn),1)*86400;
%                 StartTime = temp_LightsOut;
%                 Time=(StartTime:(1/F_samp):StartTime+(N_timeXHz-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
%             end
%             
%             
%             if 0
%             temp_LightsOut = mod(datenum(SLightsout.Start.Text(12:end-3),formatIn),1)*86400;
%             temp_LightsOn = mod(datenum(SLightsout.End.Text(12:end-3),formatIn),1)*86400;
%             N_LightsOff = (1-(datenum(SLightsout.Start.Text(12:end-3),formatIn)-datenum(SLightsout.End.Text(12:end-3),formatIn)))*86400;
%             Time_LO =(temp_LightsOut:(1/F_samp):temp_LightsOut+N_LightsOff)';
%             Time_diff = temp_LightsOut-Time(1);
%             Time_diffN = length(Time)-length(Time_LO);
%             Time = Time_LO; % Remove if starttime should not equal lights off time
%             Starttime = temp_LightsOut;
%             end
%             
%             EndTimeText = datestr(Time(end)/86400,'HH.MM.SS');
%             EndTime = mod(datenum(EndTimeText,'HH.MM.SS'),1)*86400;
%             temp_LightsOn = datenum(SLightsout.End.Text(12:end-3),formatIn);
%             
%             if 0 % if lights off is way after end time
%                 temp_LightsOut = datenum(SLightsout.Start.Text(12:end-3),formatIn);
%                 N_LightsOff = (1-(temp_LightsOut-temp_LightsOn))*86400;
%                 Time_ =(temp_LightsOut:(1/F_samp):temp_LightsOut+N_LightsOff)';
%             end
            
            iEpochStart=(1+round((HypStart(Ilist)-Time(1))/30)); %HypStart(Ilist(1))
            iEpochEnd=iEpochStart+iEpochDuration-1;
            checktemp=[iEpochStart iEpochDuration iEpochEnd EpochCodes];
            starttimediff = HypStart(1)-Time(1);
            if iEpochStart(1)<1
                %scoring starts before recording
                starttimediff = starttimediff + (1-iEpochStart(1))*30;
            end
            
            %starttimediffN=round((starttimediff/30))
            displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            Nepochs=ceil((Time(end)-Time(1))/30); %make this based on Time(end)-Time(1)...sum(iEpochDuration)
            HypnogramX=NaN*zeros(Nepochs,1);
            Hypnogram=NaN*zeros(Nepochs,1);
%             iEpochEnd=cumsum(EpochDuration);
%             iEpochStart=[0 iEpochEnd(1:end-1)]+1;
            for i=1:length(iEpochDuration)
%                 if EpochCodes(i)==2 %not sure what this was.
%                     continue
%                 end                
                li=iEpochStart(i);
                ri=iEpochEnd(i);
                if li<1
                    li=1;
                end
                if ri>length(HypnogramX)
                    ri=length(HypnogramX);
                end
                HypnogramX(li:ri)=EpochCodes(i);
            end
            
            %Convert to Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
            Hypnogram(HypnogramX==2)=4;
            Hypnogram(HypnogramX==201)=3;
            Hypnogram(HypnogramX==301)=2;
            Hypnogram(HypnogramX==302)=1;
            Hypnogram(HypnogramX==303)=0;
            Hypnogram(HypnogramX==304)=0;    
            
%             clear S
            %end BrainRT
        case 'Minerva'
            %elseif strcmp('Minerva',system)
            [~,~,eventstxt] = xlsread([directory particle char(Filenames(n,2))],1,'A2:N5001'); %max number of events = 5k
            %delete excess
            temp = cell2mat(eventstxt(:,1));
            temp = find(isnan(temp),1);
            eventstxt(temp:end,:)=[];
            EventCodesList = eventstxt(:,6);
            EventStart = cell2mat(eventstxt(:,2)); %warning all data are nearest second
            EventStart=EventStart*86400;
            EventStart(EventStart<43200)=EventStart(EventStart<43200)+86400;
            EventDuration = cell2mat(eventstxt(:,3));
            EventCategoriesLabels= {...
                'ASDA Arousal' ,...
                'Obs. Apnea',...
                'Cnt. Apnea',...
                'Cnt. Appea',...
                'Obs. Flow <70%',...
                'Obs. Flow <50%',...
                'Cnt. Flow <50%',...
                'Cnt. Flow <70%',...
                'Sustained FL'};
            EventCategoriesCodes=[1 2 3 3 4 4 6 6 11];
            EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1];
            Hypnogram = xlsread([directory char(Filenames(n,3))],1,'B:B'); %[W=0, 1 2 3, R=5]
            %Convert to Terrill code
            Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
            Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
            Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
            
        case 'Alice'
            %elseif strcmp(system,'Alice') %code is untested
            directory = char(Filenames(n,6));
            fname=char(Filenames(n,3));
            displaytext=['Get Hypnogram data:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            [num,~,~]=xlsread([directory fname],1,'A2:A86401');
            Hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15
            
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
            displaytext=['Get Events data:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            Events_fname = [directory fname];
            
            [temp,~,~] = xlsread([Events_fname],'c1:c100000'); % Read the column defining the event start time
            Nlines = length(temp);
            
            [~,~,xlseventdataall] = xlsread([Events_fname],['a2:f' num2str(Nlines+1)]); % Read the event type column
            EventCodesList = xlseventdataall(:,1);
            
            %deal with blank (NaN)
            for i=1:length(EventCodesList)
                if isnan(EventCodesList{i})
                    EventCodesList{i}='';
                end
            end
%             
%             [~,~,EventCodesList] = xlsread([Events_fname],'a:a'); % Read the event type column
%             Nlines=length(EventCodesList);
%             EventCodesList(1)=[]; % remove the header line of the data
%           
            EventStartList = xlseventdataall(:,3);
            EventStart = NaN*(ones(length(EventStartList),1));
            for i=1:length(EventStartList)
                EventStart(i) = EventStartList{i};
            end
            %[EventStart,~,~] = xlsread([Events_fname],'c:c'); % Read the column defining the event start time
            EventStart(EventStart<0.5)=EventStart(EventStart<0.5)+1;
            EventStart=EventStart*86400;
            
            EventDurationText = xlseventdataall(:,6);
            %[~,~,EventDurationText] = xlsread([Events_fname],['f2:f' num2str(Nlines)]); % Read the column defining the event duration data
            
            EventDuration=NaN*EventStart;
            for m=1:length(EventDurationText)
                if length(EventDurationText{m})>1
                    Event_Duration{m}=sscanf((EventDurationText{m}),'%f'); %finds first number in text column
                else
                    Event_Duration{m}=EventDurationText{m};
                end
                if isnumeric(Event_Duration{m})
                    EventDuration(m)=Event_Duration{m};
                end
            end
            
%             Evts.codes(2:2:end,:)=[];
%             Evts.codes(:,2:4)=[];
%             Evts.endtimes=Evts.times(2:2:end);
%             Evts.starttimes=Evts.times(1:2:end);
%             Evts.durations=Evts.endtimes-Evts.starttimes;
%             Evts.text(2:2:end,:)=[];
%             Evts.text(:,7:end)=[];
%             EventDuration = Evts.durations;
            EventCategoriesLabels={'Central apnea','Obstructive apnea','Hypopnea','�-arousal','Mixed apnea'};
            EventCategoriesCodes=[3 2 4 1 5];
            EventExact1orStartsWith2 = [1 1 1 1 1];
%             for i=1:length(EventStart)
%                 EventCodesList{i} = deblank(Evts.text(i,:));
%             end
            
        case 'Spike'
            %elseif strcmp(system,'Spike')
            %
            %load events
            echannelnameoptions.Evts={'Evts','New_Evts'};
            echannelnameoptions.Evts2={'Evts2','New_Evts2'};
            echannelnameoptions.Epochs={'Epochs'};
            
            channelnamestemp=fieldnames(echannelnameoptions);
            
            for i=1:length(channelnamestemp)
                temp=eval(['echannelnameoptions.' char(channelnamestemp(i))]);
                foundamatch=0;
                for nn=1:length(temp)
                    %Does it exist?
                    for j=1:length(w)
                        if strcmp(w(j).name,char(temp(nn)))
                            eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                            foundamatch=1;
                            break
                        end
                    end
                    if foundamatch
                        break
                    end
                end
            end
            
            %process events
            %add extra row of zeros with time = inf just in case there is an event without an end in the first of two concatenated sets of events.
            Evts.codes(size(Evts.codes,1)+1,:)=[0 0 0 0];
            Evts.times(size(Evts.times,1)+1,:)=Inf;
            Evts.text=[Evts.text;['                                                ']];
            
            %delete all start rows without codes at start of list
            while 1
                if Evts.codes(1,1)==0
                    Evts.codes(1,:)=[];
                    Evts.times(1)=[];
                    Evts.length=Evts.length-1;
                    Evts.text(1,:)=[];
                else
                    break
                end
            end
            %repeat for Evts2 separately
            if exist('Evts2','var')
                while 1
                    if Evts2.codes(1,1)==0
                        Evts2.codes(1,:)=[];
                        Evts2.times(1)=[];
                        Evts2.length=Evts.length-1;
                        Evts2.text(1,:)=[];
                    else
                        break
                    end
                end
            end
            
            %combine events lists
            if exist('Evts2','var')&&exist('Evts','var')
                
                Evts.times=[Evts.times;Evts2.times];
                Evts.codes=[Evts.codes;Evts2.codes];
                Evts.text=[Evts.text;Evts2.text];
                Evts.length=Evts.length+Evts2.length;
            end
            
            %delete second row of zeros if there are two concurrent rows of zeros
            for i=length(Evts.codes):-1:2
                if Evts.codes(i)==0&&Evts.codes(i-1)==0
                    Evts.codes(i,:)=[];
                    Evts.times(i)=[];
                    Evts.text(i,:)=[];
                    Evts.length=Evts.length-1;
                    disp('delete_row_of_zeros');
                end
            end
            
            %add a row of zeros if there is a missing row of zeros. Assume end of event is the start of the next event (or end of file/Inf).
            for i=length(Evts.codes):-1:1
                if i==length(Evts.codes)
                    if Evts.codes(i)>0
                        Evts.codes=[Evts.codes;[0 0 0 0]];
                        Evts.times=[Evts.times;Inf];
                        Evts.text=[Evts.text;['                                            ']];
                        Evts.length=Evts.length+1;
                        disp('added_row_of_zeros to the end');
                    end
                end
                if i>1
                    if Evts.codes(i)>0&&Evts.codes(i-1)>0
                        Evts.codes=[Evts.codes(1:(i-1),:);[0 0 0 0];Evts.codes(i:end,:)];
                        Evts.times=[Evts.times(1:i);Evts.times(i:end)];
                        Evts.text=[Evts.text(1:(i-1),:);['                                            '];Evts.text(i:end,:)];
                        Evts.length=Evts.length+1;
                        disp('added_row_of_zeros');
                    end
                end
            end
            
            %remove zero rows:
            Evts.codes(2:2:end,:)=[];
            Evts.codes(:,2:4)=[];
            Evts.endtimes=Evts.times(2:2:end);
            Evts.starttimes=Evts.times(1:2:end);
            Evts.durations=Evts.endtimes-Evts.starttimes;
            Evts.text(2:2:end,:)=[];
            Evts.text(:,7:end)=[];
            EventStart=Evts.starttimes + StartTime;
            EventDuration = Evts.durations;
            EventCategoriesLabels= {...
                'AR' ,...
                'ApO',...
                'ApC',...
                'HypO',...
                'M',...
                'HypC',...
                'HypO2%',...
                'HypC2%',...
                'HypOx',...
                'HypCx'};
            EventCategoriesCodes=[1 2 3 4 5 6 7 8 9 10];
            EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1];
            for i=1:length(EventStart)
                EventCodesList{i} = deblank(Evts.text(i,:));
            end
            
            Hypnogram=double(Epochs.codes(:,1));
            
            Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
            Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
            Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
            Hypnogram(Hypnogram==-3)=NaN;
         case 'NoxMat' %code is untested   
             
            echannelnameoptions.AllEvents={'AllEvents'};
            echannelnameoptions.HypnogramInfo={'HypnogramInfo'};
            
            channelnamestemp=fieldnames(echannelnameoptions);
            
            for i=1:length(channelnamestemp)
                temp=eval(['echannelnameoptions.' char(channelnamestemp(i))]);
                foundamatch=0;
                for nn=1:length(temp)
                    %Does it exist?
                    for j=1:length(w)
                        if strcmp(w(j).name,char(temp(nn)))
                            eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                            foundamatch=1;
                            break
                        end
                    end
                    if foundamatch
                        break
                    end
                end
            end
        
        
        %fix tags    
        clear tempnewlabelsNoSpace
        for i=1:size(AllEvents.Tags,1)
            temp = AllEvents.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        AllEvents.Tags = tempnewlabelsNoSpace;
        
        EventCodesList = AllEvents.Tags;
            EventStart = AllEvents.StartEndDuration(1,:); %warning all data are nearest second
            EventStart=EventStart(:);
%             EventStart=EventStart*86400;
%             EventStart(EventStart<43200)=EventStart(EventStart<43200)+86400;
            EventDuration = AllEvents.StartEndDuration(3,:); 
            EventDuration=EventDuration(:);
            EventCategoriesLabels= {...
                'arousal' ,...
                'arousal-spontaneous', ...
                'arousal-limbmovement', ...
                'arousal-plm', ...
                'arousal-respiratory', ...
                'apnea-obstructive',...
                'apnea-central',...
                'hypopnea',...
                'hypopnea-obstructive',...
                'hypopnea-central',...
                'apnea-mixed',...
                };
            EventCategoriesCodes=[1 1 1 1 1 2 3 4 4 6 5];
            EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1 1 1];
            
        clear tempnewlabelsNoSpace
        for i=1:size(HypnogramInfo.Tags,1)
            temp = HypnogramInfo.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        HypnogramInfo.Tags = tempnewlabelsNoSpace;
        
        %
%         %Arousals.StartEndDuration = Arousals.StartEndDuration';
%         Events.StartEndDuration = Events.StartEndDuration';
        HypnogramInfo.StartEndDuration = HypnogramInfo.StartEndDuration';
        
        
        
            clear tempnewlabelsNoSpace
            
           % EventStart=[Arousals.StartEndDuration(:,1) ; Events.StartEndDuration(:,1)];
           % EventDuration=[Arousals.StartEndDuration(:,3) ; Events.StartEndDuration(:,3)];
            
            
           % EventCodesList = [Arousals.Tags ; Events.Tags];
            
            
            iEpochDuration=(round(HypnogramInfo.StartEndDuration(:,3)/30));
            iEpochStart=(1+round((HypnogramInfo.StartEndDuration(:,1)-Time(1))/30)); %HypStart(Ilist(1))
            iEpochEnd=iEpochStart+iEpochDuration-1;
            EpochCodes = HypnogramInfo.StartEndDuration(:,5);
            %checktemp=[iEpochStart iEpochDuration iEpochEnd EpochCodes];
            starttimediff = HypnogramInfo.StartEndDuration(1,1)-Time(1) - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding 
            %starttimediffN=round((starttimediff/30))
            displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            Nepochs=ceil((Time(end)-Time(1))/30); %make this based on Time(end)-Time(1)...sum(iEpochDuration)
            Hypnogram=NaN*zeros(Nepochs,1);
            for i=1:length(iEpochDuration)
%                 if EpochCodes(i)==2 
%                     continue
%                 end
                Hypnogram(iEpochStart(i):iEpochEnd(i))=EpochCodes(i);
            end
    case 'RemLogic' %code is untested
            %elseif strcmp('ProfusionXML',system) %code is untested
            %ProfusionXml Hypnogram
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
      
            %Xml events
            displaytext=['Get XML Events data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            %importing
            S=xml2struct([directory fname]);
            
            %S.EventExport.Events.Event{1,i}.Type.Text
            
            if ~isfield(S.EventExport.Events,'Event') %not checked
                displaytext=['Warning: No scored events'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
            end
            
            Nevents = length(S.EventExport.Events.Event);
            
            clear EventCodesList EventStart EventStop EventDuration
            formatIn = 'HH:MM:SS.FFF';
            for i=1:Nevents
                EventCodesList{i,1} = S.EventExport.Events.Event{1,i}.Type.Text;
                
                temp = S.EventExport.Events.Event{1,i}.StartTime.Text(12:end);
                temp = mod(datenum(temp,formatIn),1)*86400;
                
                test = [floor(temp/86400*24) floor(mod(temp/86400*24,1)*60) mod(mod(temp/86400*24,1)*60,1)*60];
                if temp<86400/2,temp=temp+86400; end
                EventStart(i,1) = temp;
                   %x  y x2
                temp = S.EventExport.Events.Event{1,i}.StopTime.Text(12:end);
                temp = mod(datenum(temp,formatIn),1)*86400;
                if temp<86400/2,temp=temp+86400; end
                EventStop(i,1) = temp;
            end
            EventDuration = EventStop - EventStart;
            
            SleepEventFlag = startsWith(string(EventCodesList),'SLEEP-'); %or should be "SLEEP-" ? as per other version
            
            EventTable = table(EventCodesList,EventStart,EventStop,EventDuration,SleepEventFlag);  
            SleepTable = EventTable(SleepEventFlag==1,:);
            
            %W=0,N1=1,N2=2,N3=3,R=5
            EventCodesListUnique = unique(EventCodesList);
            EventCategoriesLabels= {...
                'AROUSAL',...
                'APNEA-CENTRAL',...
                'APNEA-MIXED',...
                'APNEA-OBSTRUCTIVE',...
                'HYPOPNEA',...
                };
            EventCategoriesCodes=[1 3 5 2 4];
            EventExact1orStartsWith2 = [2 1 1 1 1];
            
            displaytext=['Get Hypnogram data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            SleepCategoriesLabels= {...
                'SLEEP-S0',...
                'SLEEP-REM',...
                'SLEEP-S1',...
                'SLEEP-S2',...
                'SLEEP-S3',...
                };
            SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
            
            %get codes
            temp = string(SleepTable{:,1})==SleepCategoriesLabels;
            temp = temp.*SleepCodes;
            temp = max(temp')';
            SleepTable.codes = temp;
            clear temp;
            
            NepochChanges = size(SleepTable,1);
            iEpochDuration=(round(SleepTable.EventDuration/30));
            iEpochStart=(1+round((SleepTable.EventStart-Time(1))/30)); %HypStart(Ilist(1))
            iEpochEnd=iEpochStart+iEpochDuration-1;
            checktemp=[iEpochStart iEpochDuration iEpochEnd SleepTable.codes];
            
            %first epochs filled with NaNs; starttimediff is just residual after rounding 
            starttimediff = SleepTable.EventStart(1)-Time(1) - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding 
            
            %starttimediffN=round((starttimediff/30))
            displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            Nepochs=ceil((Time(end)-Time(1))/30); %make this based on Time(end)-Time(1)...sum(iEpochDuration)
            Hypnogram=NaN*zeros(Nepochs,1);
            for i=1:length(iEpochDuration)
                Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
            end
            %end ProfusionXML
            
    case 'RemLogicTextFr' %code is untested
            
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
            
            clear eventstemp
            fid = fopen([directory '\' fname]);
            i=1;
            while 1
                eventstemp{i,:} = fgetl(fid);  %read in the data
                if eventstemp{i,:}==-1
                    eventstemp(i,:) = [];
                    break
                end
                i=i+1;
            end
            fclose(fid);   %close the file

            HeaderTextList = {'Ev�nement','Dur�e[s]'};
            I=[];
            for ii=1:length(HeaderTextList)
                I=[I,contains(eventstemp,HeaderTextList{ii})];
            end
            I=sum(I,2)/length(HeaderTextList);
            NPreHeaderlines=find(I==1)-1;
            HeaderTextExact = eventstemp{I==1};
            
            % Import Table (now that we know where the data starts)
            EventTable_temp = readtable([directory fname],'Delimiter','\t','ReadVariableNames',1,'HeaderLines',NPreHeaderlines);

%             clear timetemp
%             for i=1:size(EventTable_temp,1)
%                 eventtemp=EventTable_temp.Heure_hh_mm_ss_(i}; 
%                 %Time is in sec since the previous day's midnight.
%                 tempstr=eventtemp;%eventtemp(I(1)+1:I(2)-1);
%                 I2=find(tempstr==':');
%                 timetemp(i,1)=86400*(tempstr(I2(2)+4)=='A')+43200*(tempstr(I2(2)+4)=='P')+3600*mod(str2num(tempstr(1:I2(1)-1)),12)+60*str2num(tempstr(I2(1)+1:I2(1)+2)) + 1*str2num(tempstr(I2(2)+1:I2(2)+2)) + str2num(tempstr(end-3:end));
%             end
%             timetemp(timetemp<43200)=timetemp(timetemp<43200)+86400;
%             
%             clear timechangei
%             timechangei=find([0; diff(timetemp)]<0);
%             
%             while ~isempty(timechangei)
%                 timetemp(timechangei(1):end)=timetemp(timechangei(1):end)+86400;
%                 timechangei=find([0; diff(timetemp)]<0);
%             end
            
            % This code is replacing commented section above
            StartTimeTemp = seconds(EventTable_temp.Heure_hh_mm_ss_); % seconds from midnight
            [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
            StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference

            EventDuration = EventTable_temp.Dur_e_s_;
            EventStart = StartT;
            EventStop=EventStart+EventDuration;
            EventCodesList=EventTable_temp.Ev_nement;
                   
            EventTable = table(EventCodesList,EventStart,EventStop,EventDuration); 
            SleepTable = EventTable;
            
            %W=0,N1=1,N2=2,N3=3,R=5
            EventCodesListUnique = unique(EventCodesList);  
            EventCategoriesLabels= {...
                'AROUSAL',...
                'APNEA-CENTRAL',...
                'APNEA-MIXED',...
                'APNEA-OBSTRUCTIVE',...
                'HYPOPNEA',...
                };
            EventCategoriesCodes=[1 3 5 2 4];
            EventExact1orStartsWith2 = [2 1 1 1 2];
            
            displaytext=['Get Hypnogram data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            SleepCategoriesLabels= {...
                'SLEEP-S0',...
                'SLEEP-REM',...
                'SLEEP-S1',...
                'SLEEP-S2',...
                'SLEEP-S3',...
                };
            SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
            
            %get codes
            I = string(SleepTable.EventCodesList)==SleepCategoriesLabels;
            I2 = sum(I,2)>0;
            temp = I.*SleepCodes;
            temp = max(temp')';
            SleepTable.codes = temp;
            
            SleepTable = SleepTable(I2,:);
            temp = SleepTable.EventStart(2:end) - SleepTable.EventStart(1:end-1);
            templast = Time(end)-SleepTable.EventStart(end);
            SleepTable.EventDuration(1:end-1) =temp;
            SleepTable.EventDuration(end)=templast;
            SleepTable.EventStop=SleepTable.EventStart+SleepTable.EventDuration;
            clear temp;
            
            NepochChanges = size(SleepTable,1);
            iEpochDuration=(round(SleepTable.EventDuration/30));
            iEpochStart=(1+round((SleepTable.EventStart-Time(1))/30)); %HypStart(Ilist(1))
            iEpochEnd=iEpochStart+iEpochDuration-1;
            checktemp=[iEpochStart iEpochDuration iEpochEnd SleepTable.codes];
            
            %first epochs filled with NaNs; starttimediff is just residual after rounding 
            starttimediff = SleepTable.EventStart(1)-Time(1) - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding 
            
            %starttimediffN=round((starttimediff/30))
            displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            Nepochs=ceil((Time(end)-Time(1))/30); %make this based on Time(end)-Time(1)...sum(iEpochDuration)
            Hypnogram=NaN*zeros(Nepochs,1);
            for i=1:length(iEpochDuration)
                Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
            end
            %end ProfusionXML
        case 'Deltamed' %code is untested
            
            directory = char(Filenames(n,5));
            fname = char(Filenames(n,2));
            
            clear eventstemp
            fid = fopen([directory '\' fname]);
            i=1;
            while 1
                eventstemp{i,:} = fgetl(fid);  %read in the data
                if eventstemp{i,:}==-1
                    eventstemp(i,:) = [];
                    break
                end
                i=i+1;
            end
            fclose(fid);   %close the file

            HeaderTextList = {'Dur�e','Annotation'};
            I=[];
            for ii=1:length(HeaderTextList)
                I=[I,contains(eventstemp,HeaderTextList{ii})];
            end
            I=sum(I,2)/length(HeaderTextList);
            NPreHeaderlines=find(I==1)-1;
            HeaderTextExact = eventstemp{I==1};
            
            % Use regex to format the data
            eventstempNew = eventstemp;
            for rownum = NPreHeaderlines+1:length(eventstemp)
                temprow = eventstempNew{rownum,1};
                % insert duration into spots with missing data
                newtemprow = regexprep(temprow,'\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s','     00:00:00  ');

                % Make 5 space gaps into 6 space gaps
                newtemprow = regexprep(newtemprow,'([0-9])(\s\s\s\s\s)([0-9])','$1      $3');

                % Make 4 space gaps into 6 space gaps
                newtemprow = regexprep(newtemprow,'([0-9])(\s\s\s\s)([0-9])','$1      $3');

                % Make 4 space gaps into 6 space gaps
                newtemprow = regexprep(newtemprow,'([0-9])(\s\s\s)([0-9])','$1      $3');

%                 % Remove + signs
%                 newtemprow = regexprep(newtemprow,' \+','');
% 
                % Specific for one of the datasets
                newtemprow = regexprep(newtemprow,'(JMB1)(\s\s)(JMB2)','$1$3');

                % Replace ALL single spaced gaps -> replaces two commented
                % out secions above
                % finds spaces with NO space behind and NO space
                % ahead
                newtemprow = regexprep(newtemprow,'(?<!\s)(\s)(?!\s)','');

                % Replace h, m with colons
                newtemprow = regexprep(newtemprow,'([0-9])([hm])([0-9])','$1:$3');

                % Remove s
                newtemprow = regexprep(newtemprow,'([0-9])([s])','$1');

                % add back to events temp
                eventstempNew{rownum,1} = newtemprow;
            end

            % save as text
            fid2 = fopen([directory,fname(1:end-4),'_new.txt'],'w');
            for row = 1:length(eventstempNew)
                fprintf(fid2,'%s\n',eventstempNew{row,:});
            end
            fclose(fid2);

            % Import Table (now that we know where the data starts)
            EventTable_temp = readtable([directory,fname(1:end-4),'_new.txt'],'Delimiter',' ',...
                'ReadVariableNames',1,'HeaderLines',NPreHeaderlines,...
                'MultipleDelimsAsOne',1);
            
            % This code is replacing commented section above
            StartTimeTemp = seconds(duration(EventTable_temp.Heurer_elle)); % seconds from midnight
            [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
            StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference

            EventDuration = seconds(duration(EventTable_temp.Dur_e));
            EventStart = StartT;
            EventStop=EventStart+EventDuration;
            EventCodesList=EventTable_temp.Annotation;
                   
            EventTable = table(EventCodesList,EventStart,EventStop,EventDuration); 
            SleepTable = EventTable;
            
            %W=0,N1=1,N2=2,N3=3,R=5
            EventCodesListUnique = unique(EventCodesList);  
            EventCategoriesLabels= {...
                'Arousal',...
                'Apn�eCentrale',...
                'Apn�eMixte',...
                'Apn�eObstructive',...
                'Hypopn�e',...
                'hypopn�eCentrale',...
                };
            EventCategoriesCodes=[1 3 5 2 4 4];
            EventExact1orStartsWith2 = [2 2 2 2 2 2];
            
            displaytext=['Get Hypnogram data: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            % Import sleep stage table - alway starts at row 6
            SleepStageTable_temp = readtable([directory,fname(1:end-11),'.txt'],'Delimiter','\t',...
                'ReadVariableNames',0,'HeaderLines',5);
            
            % Get start time of sleep stages seconds
            StartTimeTemp = seconds(SleepStageTable_temp.Var1);
            [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
            StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference
            StageStartTemp = StartT;
            
            SleepStageList = SleepStageTable_temp.Var2;
            StageDurationTemp = nan(size(StageStartTemp));
            SleepStageTable = table(SleepStageList,StageStartTemp,StageDurationTemp);
            SleepStageTable.Properties.VariableNames{'StageDurationTemp'} = 'EventDuration';
            SleepStageTable.Properties.VariableNames{'StageStartTemp'} = 'EventStart';
            
            SleepCategoriesLabels= {...
                'Veille',...
                'Stade 1',...
                'Stade 2',...
                'Stade 3',...
                };
            SleepCodes = [4 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
            
            %get codes
            I = string(SleepStageTable.SleepStageList)==SleepCategoriesLabels;
            I2 = sum(I,2)>0;
            temp = I.*SleepCodes;
            temp = max(temp')';
            SleepStageTable.codes = temp;
            
            SleepStageTable = SleepStageTable(I2,:);
            temp = SleepStageTable.EventStart(2:end) - SleepStageTable.EventStart(1:end-1);
            templast = Time(end)-SleepStageTable.EventStart(end);
            SleepStageTable.EventDuration(1:end-1) =temp;
            SleepStageTable.EventDuration(end)=templast;
            SleepStageTable.EventStop=SleepStageTable.EventStart+SleepStageTable.EventDuration;
            clear temp;
            
            NepochChanges = size(SleepStageTable,1);
            iEpochDuration=(round(SleepStageTable.EventDuration/30));
            iEpochStart=(1+round((SleepStageTable.EventStart-Time(1))/30)); %HypStart(Ilist(1))
            iEpochEnd=iEpochStart+iEpochDuration-1;
            checktemp=[iEpochStart iEpochDuration iEpochEnd SleepStageTable.codes];
            
            %first epochs filled with NaNs; starttimediff is just residual after rounding 
            starttimediff = SleepStageTable.EventStart(1)-Time(1) - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding 
            
            %starttimediffN=round((starttimediff/30))
            displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            Nepochs=ceil((Time(end)-Time(1))/30); %make this based on Time(end)-Time(1)...sum(iEpochDuration)
            Hypnogram=NaN*zeros(Nepochs,1);
            for i=1:length(iEpochDuration)
                Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepStageTable.codes(i);
            end
            %end ProfusionXML        
    end   
    %% Plot Hypnogram
    
    if 0 %Plot Hypnogram
        xtime=(Time(1):30:(Time(1)+30*(length(Hypnogram)-1)));
        figure(); stairs(xtime/86400,Hypnogram);
        datetick('x','HH:MMPM');
    end
    
    %%  General Code: Warning
    if sum(Hypnogram==4)==length(Hypnogram)
        displaytext=['Warning: Entire hypnogram is wake'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    %% General Code: Make Events, Arousals, Epochs signals
    %EventCodes=NaN*EventStart;
    
    EventStart=EventStart(:);
    EventDuration=EventDuration(:);
    
    clear Evts
    EventCodesList=EventCodesList(:);
    Evts.EventCodesList = EventCodesList;
    Evts.EventCategoriesLabels = EventCategoriesLabels;
    Evts.Table1 = table(EventCodesList,EventStart,EventDuration);
    
    clear match
    for m=1:length(EventCategoriesLabels)% each column is an event type
        if EventExact1orStartsWith2(m)==1 %exact
            match(:,m)=string(EventCategoriesLabels(m))==string(EventCodesList');
        elseif EventExact1orStartsWith2(m)==2 %startswith
            match(:,m)=startsWith(string(EventCodesList'),string(EventCategoriesLabels(m)));
        end
    end
    temp = match.*EventCategoriesCodes;
    temp(temp==0)=NaN;
    EventCodes=min(temp'); %just in case two criteria are detected, will leave the lowest code value in list
    
    Evts.Table1.EventCodes=EventCodes(:);
    
    %check for events that should be in the list
    if ~isempty(EventCodesList(isnan(EventCodes)))
        disp(['Events found but unused: ' strjoin(unique(EventCodesList(isnan(EventCodes))),', ')]);
    else
        disp(['Events OK']);
    end

    %     for m=1:length(EventStart) %
    %         temp = strcmp(EventCategoriesLabels,EventCodesList(m));
    %         if sum(temp)>0
    %             EventCodes(m)=EventCategoriesCodes(temp);
    %         end
    %     end
    I=isnan(EventCodes);
    EventCodes(I)=[]; 
    EventStart(I)=[]; 
    EventDuration(I)=[];
    
    %remove very short resp events
    I=EventDuration(:)<2&EventCodes(:)>1;
    EventCodes(I)=[]; 
    EventStart(I)=[]; 
    EventDuration(I)=[];  
    
    if sum(EventCodes==1)==0
        displaytext=['Warning: No scored arousals'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    if std(EventDuration(EventCodes==1))<1
        displaytext=['Warning: Arousals are all similar durations'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    %initialize
    EventsResp=0*Time;
    EventsAr=0*Time;
    for m=1:length(EventStart) %Arousals
        lefti=round((EventStart(m)-StartTime)*F_samp)+1;
        righti=lefti+round((EventDuration(m))*F_samp);
        if lefti<1, lefti=1; end
        if righti>length(Time), righti=length(Time); end
        if EventCodes(m)>1
            EventsResp(lefti:righti)=EventCodes(m);
        elseif EventCodes(m)==1
            EventsAr(lefti:righti)=1;
        end
    end
    
    if 1 %set NaN, others, to 8
        Hypnogram(isnan(Hypnogram))=8; %no stage code
        Hypnogram(Hypnogram<0)=8; %no stage code
    end
    Hypnogram_t = StartTime+starttimediff+(0:30:(30*(length(Hypnogram)-1)))';
    Epochs = interp1(Hypnogram_t+15,Hypnogram,Time,'nearest','extrap');
    
    Evts.Hypnogram = Hypnogram;
    Evts.Hypnogram_t = Hypnogram_t;
%     if 1
%         Epochs(isnan(Epochs))=8;
%     end
    %figure(); plot(Time,Epochs);
    
    ChannelsList = [ChannelsList,{'Epochs','EventsAr','EventsResp'}];
    Channels_Fs = [Channels_Fs; F_samp;F_samp;F_samp];

    %% Sampling Frequency Errors, Time Drift
    %settings.FsLeak=-0.00008; %lag of EDF time vs actual in proportion
    if settings.FsLeak~=0
        disp('Warning: fixing time leak');
        Channels_Fs2 = Channels_Fs;
        Channels_Fs_fix = Channels_Fs*0 + 1;
        I=find(sum(string(ChannelsList)==string({'Epochs','EventsAr','EventsResp'})'));
        Channels_Fs_fix(I)=0;
        Channels_Fs2(Channels_Fs_fix==1) = Channels_Fs2(Channels_Fs_fix==1)*(1-settings.FsLeak);
        Channels_Fs = Channels_Fs2;
    end
    
    %% Snore analysis
    if settings.AnalyzeSnore
        fnameSnore = [fname(1:end-4),'Snore.mat'];
        filehandleSnore = matfile([directory fnameSnore]);
        Snore = filehandleSnore.snore; % need to change this to be flexible with naming of snore variable
        SnoreFs = 1/Snore.interval;
        % process snore data
        SnoreStruct = ProcessSnore(Snore);
        WinVar = SnoreStruct.WinVar; WinVar50 = SnoreStruct.WinVar50;
        frmnFreq1 = SnoreStruct.frmnFreq1; frmnFreq2 = SnoreStruct.frmnFreq2;
        frmnFreq3 = SnoreStruct.frmnFreq3; p70to250  = SnoreStruct.p70to250;
        p250to500 = SnoreStruct.p250to500; p500to750 = SnoreStruct.p500to750;
        p750to1k = SnoreStruct.p750to1k; p1kto1p5k = SnoreStruct.p1kto1p5k; 
        p1p5kto2k = SnoreStruct.p1p5kto2k; centroid = SnoreStruct.centroid;
        SnoreChannelsList = {'WinVar','WinVar50','p70to250','p250to500','p500to750',...
            'p750to1k','p1kto1p5k','p1p5kto2k','centroid','frmnFreq1','frmnFreq2','frmnFreq3'};
        SnoreChannels_Fs = ones(size(SnoreChannelsList))*F_samp;
        ChannelsList = [ChannelsList, SnoreChannelsList];
        Channels_Fs = [Channels_Fs; SnoreChannels_Fs'];
    end 
    
    %% EEG power analysis, EKG
    displaytext=['EEG processing'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
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
        Fs_EEG(i)=Channels_Fs(find(strcmp(ChannelsList,EEGsignals{i})==1));
    end
    Fs_EEG_ = mode(Fs_EEG); %average EEG sampling rate
    
    % build in option to resample or discard EEG signals that have a
    % different sampling rate to the average
    % Resample channels
    
    TimeEEG=(StartTime:(1/Fs_EEG_):StartTime+(length(eval(EEGsignals{1}))-1)*(1/Fs_EEG_))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(EEGsignals)
        if Fs_EEG(i)~=Fs_EEG_
            displaytext=['Resampling: ' EEGsignals{i} ' from ' num2str(Fs_EEG(i)) ' to ' num2str(Fs_EEG_) 'Hz'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            %eval([EEGsignals{i} ' = resample(' EEGsignals{i} ',Fs_EEG_,Channels_Fs(i));']); %only works if Fsamp / Fs are integer multiples of each other
           TimeTemp=(StartTime:(1/Fs_EEG(i)):StartTime+(length(eval(EEGsignals{i}))-1)*(1/Fs_EEG(i)))'; % This is the time vector associated with the 100Hz Flow data.
           Temp = interp1(TimeTemp,eval(EEGsignals{i}),TimeEEG,'linear');
           eval([EEGsignals{i} '= Temp;']);
           Fs_EEG(i)=Fs_EEG_;
           Channels_Fs(find(strcmp(EEGsignals{i},ChannelsList)))=Fs_EEG_;
        end
    end
    
    TimeEEG=(StartTime:(1/Fs_EEG_):StartTime+(length(eval(EEGsignals{1}))-1)*(1/Fs_EEG_))'; % This is the time vector associated with the 100Hz Flow data.
    
    fftlengthoptions = 2.^[6:9];
    [~,bestoptioni] = min(abs(Fs_EEG_ - fftlengthoptions));
    
    
    
    if exist('EKG','var')
        %EKG_original=EKG; %not needed
        Fs_EKG = Channels_Fs(find(strcmp(ChannelsList,'EKG')==1));
        if Fs_EKG~=Fs_EEG_
            displaytext=['Resampling EKG from ' num2str(Fs_EKG) ' to ' num2str(Fs_EEG_) 'Hz'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
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
    skipEKG=0;
    if exist('EKG','var')&&(1-skipEKG)
        displaytext=['EKG processing'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        [ECG_peak_i,HR] = EKGpeakdetection(EKGeeg,TimeEEG,1./Fs_EEG_,0,settings.EKGplot);
        ChannelsList = [ChannelsList,'HR'];
        Channels_Fs = [Channels_Fs;Fs_EEG_];
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
        
        disp(WakeSleepInfo.AUC_M);
        %Clean EEG signals
        for i=1:length(EEGsignals)
            eval([EEGsignals{i} '=EEGsignalsOut{i};'])
        end
    
        if 1
            %best EEG
            [~,besti]=max(WakeSleepInfo.AUC_M);
            EEG = eval([EEGsignals{besti}]);
            ChannelsList = [ChannelsList,{'EEG'}];
            Channels_Fs = [Channels_Fs;Fs_EEG_];
        end
    
            %WakeSleep
    %         temp = WakeSleep.values;
    %         clear WakeSleep;
    %         WakeSleep = temp;
        ChannelsList = [ChannelsList,'WakeSleep',Powersignals];
        Channels_Fs = [Channels_Fs ; F_samp ;F_samp ;F_samp ;F_samp ;F_samp ;F_samp];
        %Powersignals
        for i=1:length(Powersignals)
            eval([Powersignals{i} '=PowersignalsOut{i};'])
        end
        
        if settings.allpowers
        Powersignalsalllist=[];
        for j=1:length(EEGsignals)
            for i=1:length(PowersignalsOutAll{j})
                eval([Powersignals{i} num2str(j) '=PowersignalsOutAll{j}{i};']);
                Powersignalsalllist{length(Powersignalsalllist)+1} = [Powersignals{i} num2str(j)];
            end
        end
        PowersignalsalllistFs = F_samp * ones(length(Powersignalsalllist),1);
        
        ChannelsList = [ChannelsList,Powersignalsalllist];
        Channels_Fs = [Channels_Fs;PowersignalsalllistFs];
        
        end
        
        %%AA arousal intensity 
        if ~settings.skiparousalintensity
        try 
            [ArIntensityOrig,ArIntensity]=ArousalIntensityRun(EEGsignals,EventsAr,Epochs,Time,TimeEEG,WakeSleepInfo.AUC_M);
            ChannelsList = [ChannelsList,'ArIntensityOrig','ArIntensity'];
            Channels_Fs = [Channels_Fs ; F_samp ;F_samp ];
        catch meAI
            disp('Failed Arousal Intensity')
        end
        end
        
        clear EEGsignalsOut PowersignalsOut
    
    catch me 
        disp(me.message);
        EEG = eval([EEGsignals{1}]);
        ChannelsList = [ChannelsList,{'EEG'}];
        Channels_Fs = [Channels_Fs;Fs_EEG_];
        WakeSleepInfo.AUC_M=NaN;
    end
    
    else %No EEG present: Empty list of EEG channels
        disp('Warning: No EEG');
        WakeSleepInfo.AUC_M=NaN;
        if exist('EKG','var')
            displaytext=['EKG processing'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            try
            Fs_EKG = Channels_Fs(find(strcmp(ChannelsList,'EKG')==1));
            TimeTempEKG=(StartTime:(1/Fs_EKG):StartTime+(length(EKG)-1)*(1/Fs_EKG))'; % This is the time vector associated with the 100Hz Flow data.
            [ECG_peak_i,HR] = EKGpeakdetection(EKG,TimeTempEKG,1./Fs_EKG,0,settings.EKGplot);
            ChannelsList = [ChannelsList,'HR'];
            Channels_Fs = [Channels_Fs;Fs_EKG];
            catch me
                disp('failed EKG analysis');
            end
        end
    end
    clear TimeTempEKG

    

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
            Fs_EKG = Channels_Fs(find(strcmp(ChannelsList,'EKG')==1));
            if Fs_EKG~=F_samp
                'resampling EKG to match Resp'
                temp = resample(EKG,F_samp,Fs_EKG);
                deltalength = length(temp)-length(Time);
                if deltalength>0
                    'removing samples in EKG to match Resp'
                    temp(end-deltalength+1:end)=[];
                elseif deltalength<0
                    'adding samples in EKG to match Resp'
                    temp(end:end-deltalength)=temp(end);
                end
                Pes_clean = PesRemoveEKG(Pes,temp,1/F_samp,ECG_peak_i);
            else
                Pes_clean = PesRemoveEKG(Pes,EKG,1/F_samp,ECG_peak_i);
            end
            Pes = Pes_clean;
        end
    end
    catch me
        disp(me.message)
        'failed Pes decontaminate'
    end
    
    %% Remove from channel list list
    ChannelsRemoveList={};
    if 0
    ChannelsRemoveList={'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12'};
    end
    
    for i=length(ChannelsList):-1:1
        for j=length(ChannelsRemoveList):-1:1
            if strcmp(ChannelsRemoveList{j},ChannelsList{i})==1
                eval(['clear ' ChannelsList{i}]);
                ChannelsList(i)=[];
                Channels_Fs(i)=[];
            end
        end
    end
    
    %% Position channel fix using text
    % Fix Pos data if necessary
    % assumes times are in time since start recording in sec
    textfilename=[directory fname(1:end-4) '_pos.txt'];
    if exist(textfilename,'file')==2
        if ~(exist('Position','var')==1)
            Position = Time*0 + NaN;
        end
        displaytext = 'Correcting position data';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        %Position.values_original=Position.values;
        %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
        % textread was not recommended, upgraded to dlmread
        [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
        for i=1:(length(col1)-1)
            Position(Time-Time(1)>=col1(i)&Time-Time(1)<col1(i+1))=col2(i);
        end
        Position(Time-Time(1)>=col1(end))=col2(end);
        if sum(strcmp(ChannelsList,'Position'))==0
            ChannelsList = [ChannelsList,'Position'];
            Channels_Fs = [Channels_Fs;F_samp];
        end
    end

    
    %% Remove artifact in all signals using text files
    if 1
        for j=1:length(ChannelsList)
            textfilename=[directory fname(1:end-4) '_' ChannelsList{j} '_art.txt'];
            dt = 1/Channels_Fs(j);
            if exist(textfilename,'file')==2
                displaytext = ['Removing artifact from ' ChannelsList{j}];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                %[col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
                % textread was not recommended, upgraded to dlmread
                [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
                for i=1:length(col1)
                    %dttemp=eval([channels{j} '.interval']);
                    lefti=round((col1(i))/dt+1);
                    if lefti<1, lefti=1; end
                    righti=round((col2(i))/dt+1);
                    if righti>length(eval(ChannelsList{j})), righti=length(eval(ChannelsList{j})); end
                    eval([ChannelsList{j} '(lefti:righti)=NaN;']);
                end
                Percent_removed = sum(isnan(eval(ChannelsList{j})))/length(eval(ChannelsList{j}))*100;
                displaytext = [num2str(Percent_removed), '% of ', ChannelsList{j}, ' removed as artifact']; 
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                %figure(); plot(Time,eval([ChannelsList{j}]));
            end
        end
    end
    
    %% SpO2 artifact (set to zeros, NaN makes future artifact detection fail) 
    if 1
        try
            I = SpO2>101|SpO2<40;
            Ileft = [0;I(1:end-1)];
            Iright = [I(2:end);0];
            SpO2(I==1)=0;
            SpO2(Ileft==1)=0;
            SpO2(Iright==1)=0;
            %         figure()
            %         plot(SpO2)
        catch me
        end
    end
    
    
       
    %% Resample channels
    %SpO2_ = SpO2;
    ChannelsFs = Channels_Fs;
    for i=2:length(ChannelsList) %skip "Time"
        if Channels_Fs(i)~=F_samp
            displaytext=['Resampling: ' ChannelsList{i} ' from ' num2str(Channels_Fs(i)) ' to ' num2str(F_samp) 'Hz'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            if strcmp(ChannelsList(i),'SpO2')||strcmp(ChannelsList(i),'Position') %use nearest rather than interp
                TimeTemp=(StartTime:(1/Channels_Fs(i)):StartTime+(length(eval(ChannelsList{i}))-1)*(1/Channels_Fs(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(ChannelsList{i}),Time,'nearest');
                eval([ChannelsList{i} '= Temp;']);
                tol(i)=NaN;
            elseif 1%(F_samp/Channels_Fs(i))>2
                TimeTemp=(StartTime:(1/Channels_Fs(i)):StartTime+(length(eval(ChannelsList{i}))-1)*(1/Channels_Fs(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(ChannelsList{i}),Time,'linear');
                eval([ChannelsList{i} '= Temp;']);
                tol(i)=NaN;
            else
                tol(i)=0.00000001;
                while 1
                    [Q,P] = rat(Channels_Fs(i)/F_samp,tol(i));
                    if (Q*P)>2^31
                        tol(i)=tol(i)*10;
                    else
                        break
                    end
                end
                %eval([ChannelsList{i} ' = resample(' ChannelsList{i} ',round(F_samp),round(Channels_Fs(i)));']); %only works if Fsamp / Fs are integer multiples of each other
                eval([ChannelsList{i} ' = resample(' ChannelsList{i} ',P,Q);']); %only works if Fsamp / Fs are integer multiples of each other
            end
            displaytext=['Resampling: ' ChannelsList{i} ' from ' num2str(Channels_Fs(i)) ' to ' num2str(F_samp)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            Channels_Fs(i)=F_samp;
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
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
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
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
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
    fclose all;
    %pack;
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ------  Unused functions to merge with GeneralImport  -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



%% Alice
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Alice(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    %EDF Start time
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get recording start time from EDF: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    fid = fopen([directory '\' fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the EDF files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    %F_samp=100;
    
    EDFfilenamedir = [directory '\' fname];
    
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time;
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
    end

    %Alice Hypnogram
    
    directory = char(Filenames(n,6));
    fname=char(Filenames(n,3));
    displaytext=['Get Hypnogram data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    [num,~,~]=xlsread([directory '\' fname],1,'A:A');
    hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Alice Events
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    Events_fname = [directory '\' fname];
    
    [~,Event_Type,~] = xlsread([Events_fname],'a:a'); % Read the event type column
    Nlines=length(Event_Type);
    Event_Type(1)=[]; % remove the header line of the data
    
    [Event_StartTime,~,~] = xlsread([Events_fname],'c:c'); % Read the column defining the event start time
    Event_StartTime(Event_StartTime<0.5)=Event_StartTime(Event_StartTime<0.5)+1;
    Event_StartTime=Event_StartTime*86400;
    
    [~,~,EventDurationText] = xlsread([Events_fname],['f2:f' num2str(Nlines)]); % Read the column defining the event duration data
    for m=1:length(EventDurationText)
        if length(EventDurationText{m})>1
            Event_Duration{m}=sscanf((EventDurationText{m}),'%f'); %finds first number in text column
        else
            Event_Duration{m}=EventDurationText{m};
        end
    end
    
    %Add Events data in additional colums of the matrix:
    %**************************************************************************
    EventCategories={'Central apnea','Obstructive apnea','Hypopnea','Desaturation','�-arousal','Mixed apnea'};
    EventCols=[10 11 12 13 14 15];
    DataEventHypnog_Mat(:,10:15)=0; %initialize
    
    for m=1:length(Event_Type)
        for i=1:length(EventCategories)
            if strcmp(Event_Type{m},EventCategories{i})==1
                lefti=round((Event_StartTime(m)-StartTime)*F_samp)+1;
                righti=lefti+round(Event_Duration{m}*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
            end
        end
    end
    
    if 0 %check plot
        figure(1)
        plotcols=[2,11,12,3,9,10];
        downsamplefactor=20;
        hold('off')
        count=1;
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
        end
        linkaxes(ax,'x');
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','�-arousal','Mixed apnea'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                   3                     2        9           10        11                    12         13       14         4           5       6       7          8      15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive(/Mixed) 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Xml files from Profusion
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads,ColumnHeadsList,WakeSleepInfo] = ProfusionXML(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    %% EDF Start time
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get recording start time from EDF: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    fid = fopen([directory '\' fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
    if StartTime<43200, StartTime=StartTime+86400; end
    
    %% Signals
    %This section imports respiratory flow and other polysomnography signal data from the EDF files
    ChannelsList={'Pnasal','Thorax','Abdomen','SpO2','EEG1','EEG2','EEG3','EEG4','Position','Pmask','EKG'};
    Channels_Fs = NaN*ones(length(ChannelsList),1);
    %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
    ChannelFound = zeros(1,length(ChannelsList));
    EDFfilenamedir = [directory '\' fname];
    
    for i=1:length(ChannelsList)
        displaytext=['Collecting channel:' ChannelsList{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' ChannelsList{i} ',Channels_Fs(i),~,~,LabelTemp,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            displaytext=['Found channel labelled:' LabelTemp ' at ' num2str(Channels_Fs(i)), ' Hz'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            ChannelFound(i)=1;
        catch me
            displaytext=['No channel:' ChannelsList{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            %eval([ChannelsList{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    
    ChannelsList(ChannelFound==0)=[];
    Channels_Fs(ChannelFound==0)=[];
    %base the study duration on Pnasal signal (presumably present)
    N_flow=length(Pnasal);
    Fs_Flow=Channels_Fs(find(strcmp(ChannelsList,'Pnasal')==1));
    N_timeXHz = round((N_flow/Fs_Flow*F_samp));
    
    %Time array
    Time=(StartTime:(1/F_samp):StartTime+(N_timeXHz-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    ChannelsList = ['Time',ChannelsList];
    Channels_Fs = [F_samp Channels_Fs];
    
    %% Epochs and events
    
    %Xml Hypnogram
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get Hypnogram data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    S=xml2struct([directory '\' fname]);
    
    Nepochs = length(S.CMPStudyConfig.SleepStages.SleepStage);
    clear hypnogram
    hypnogram=NaN*zeros(Nepochs,1);
    for i=1:Nepochs
        hypnogram(i) = str2num(S.CMPStudyConfig.SleepStages.SleepStage{1,i}.Text);
    end %W=0,N1=1,N2=2,N3=3,R=5
    
    if sum(hypnogram==0)==length(hypnogram)
        displaytext=['Warning: Entire hypnogram is wake'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    hypnogram(hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
    hypnogram(hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
    hypnogram=3-hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    Epochs = interp1(hypnogram_t+15,hypnogram,Time,'nearest','extrap');
    
    %     clear Pnasal_Time; %save memory
    %
    %     DataEventHypnog_Mat(:,9)=EpochsXHz;
    %     clear EpochsXHz;
    %
    %     %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    %
    %Xml events
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    if ~isfield(S.CMPStudyConfig.ScoredEvents,'ScoredEvent')
        displaytext=['Warning: No scored events'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    Nevents = length(S.CMPStudyConfig.ScoredEvents.ScoredEvent);
    
    clear Events
    for i=1:Nevents
        Events{i}.name = S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Name.Text;
        Events{i}.start = StartTime + str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Start.Text);
        Events{i}.duration = str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Duration.Text);
    end %W=0,N1=1,N2=2,N3=3,R=5
    
    if 0
        for i=1:Nevents
            [num2str(i) ':' Events{i}.name]
        end %
    end
    
    %Add Events data in additional columns of the matrix:
    %**************************************************************************
    EventCategories={'Central Apnea','Obstructive Apnea','Mixed Apnea','Hypopnea','Unsure'} %Unsure is hypopnea in MESA
    
    %EventCols=[10 11 15 12 14 12]; %Mixed apneas are no longer treated as Obstructive apneas here. Unsure is hypopnea.
    EventCodes=[3 2 5 4 4];
    %initialize
    EventsResp=0*Time;
    EventsAr=0*Time;
    for m=1:Nevents %Arousals
        for i=1:length(EventCategories)
            strlength=length(EventCategories{i});
            if length(Events{m}.name)>=strlength&&(strcmp(Events{m}.name(1:strlength),EventCategories{i})==1)
                lefti=round((Events{m}.start-StartTime)*F_samp)+1;
                righti=lefti+round((Events{m}.duration)*F_samp);
                EventsResp(lefti:righti)=EventCodes(i);
            end
        end
    end
    for m=1:Nevents %Resp Events
        %for i=1:length(EventCategories)
        strlength=length('Arousal');
        if length(Events{m}.name)>=strlength&&(strcmp(Events{m}.name(1:strlength),'Arousal')==1)
            lefti=round((Events{m}.start-StartTime)*F_samp)+1;
            righti=lefti+round((Events{m}.duration)*F_samp);
            EventsAr(lefti:righti)=1;
        end
        %end
    end
    
    ChannelsList = [ChannelsList,{'Epochs','EventsAr','EventsResp'}];
    Channels_Fs = [Channels_Fs F_samp F_samp F_samp];
    
    %% EEG power analysis, EKG
    EEGsignals = {'EEG1','EEG2','EEG3','EEG4'};
    temp=[];
    for i=1:length(EEGsignals)
        if exist(EEGsignals{i},'var')
            temp{length(temp)+1}=EEGsignals{i};
        end
    end
    EEGsignals=temp; %ovewrite
    %Sampling rate of each EEG
    Fs_EEG=[];
    for i=1:length(temp)
        Fs_EEG(i)=Channels_Fs(find(strcmp(ChannelsList,EEGsignals{i})==1));
    end
    Fs_EEG_ = mode(Fs_EEG); %average EEG sampling rate
    % build in option to resample or discard EEG signals that have a
    % different sampling rate to the average
    
    TimeEEG=(StartTime:(1/Fs_EEG_):StartTime+(length(eval(EEGsignals{1}))-1)*(1/Fs_EEG_))'; % This is the time vector associated with the 100Hz Flow data.
    
    
    fftlengthoptions = 2.^[6:9];
    [~,bestoptioni] = min(abs(Fs_EEG_ - fftlengthoptions));
    
    if exist('EKG')
        %EKG_original=EKG; %not needed
        Fs_EKG = Channels_Fs(find(strcmp(ChannelsList,'EKG')==1));
        if Fs_EKG~=Fs_EEG_
            'resampling EKG to match EEG'
            EKG = resample(EKG,Fs_EEG_,Fs_EKG); %untested
        end
    end
    
    BetaPowerSettings.polyorder=1;
    BetaPowerSettings.suf1 = ''; %'.values'
    BetaPowerSettings.timeeegstr = 'TimeEEG';
    BetaPowerSettings.timestr = 'Time';
    BetaPowerSettings.flowstr = 'Pnasal';
    BetaPowerSettings.spo2str = 'SpO2';
    BetaPowerSettings.arstr = 'EventsAr';
    BetaPowerSettings.epochsstr = 'Epochs';
    BetaPowerSettings.fft_length = fftlengthoptions(bestoptioni); %256 %Best: 1710,1723:256, 1343:512
    BetaPowerSettings.scoredarousalsinwake = 0; %1 for Lauren's scoring.
    
    
    %EKG analysis
    
    [ECG_peak_i,HR] = EKGpeakdetection(EKG,TimeEEG,1./Fs_EEG_,0,0);
    
    [WakeSleep,WakeSleepInfo,EEGsignals,EEGsignalsOut,Powersignals,PowersignalsOut] = BetaPowerRun(EEGsignals,1/Fs_EEG_,BetaPowerSettings,ECG_peak_i);
    
    
    
    %Clean EEG signals
    for i=1:length(EEGsignals)
        eval([EEGsignals{i} '=EEGsignalsOut{i};'])
    end
    
    if 1
        %best EEG
        [~,besti]=max(WakeSleepInfo.AUC_M);
        EEG = eval([EEGsignals{besti}]);
        ChannelsList = [ChannelsList,{'EEG'}];
        Channels_Fs = [Channels_Fs,Fs_EEG_];
    end
    
    %WakeSleep
    %         temp = WakeSleep.values;
    %         clear WakeSleep;
    %         WakeSleep = temp;
    ChannelsList = [ChannelsList,'WakeSleep',Powersignals,'HR'];
    Channels_Fs = [Channels_Fs F_samp F_samp F_samp F_samp F_samp F_samp Fs_EEG_];
    %Powersignals
    for i=1:length(Powersignals)
        eval([Powersignals{i} '=PowersignalsOut{i};'])
    end
    
    clear TimeEEG EEGsignalsOut PowersignalsOut
    %% Remove from channel list list
    ChannelsRemoveList={'EEG1','EEG2','EEG3','EEG4','EKG'};
    
    for i=length(ChannelsList):-1:1
        for j=length(ChannelsRemoveList):-1:1
            if strcmp(ChannelsRemoveList{j},ChannelsList{i})==1
                eval(['clear ' ChannelsList{i}]);
                ChannelsList(i)=[];
                Channels_Fs(i)=[];
            end
        end
    end
    
    
    %% Resample channels
    Channels_Fs_original = Channels_Fs;
    for i=1:length(ChannelsList)
        if Channels_Fs(i)~=F_samp
            
            if strcmp(ChannelsList(i),'SpO2')||strcmp(ChannelsList(i),'Position') %use nearest rather than interp
                TimeTemp=(StartTime:(1/Channels_Fs(i)):StartTime+(length(eval(ChannelsList{i}))-1)*(1/Channels_Fs(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(ChannelsList{i}),Time,'nearest');
                eval([ChannelsList{i} '= Temp;']);
                tol(i)=NaN;
            else %(F_samp/Channels_Fs(i))>2
                TimeTemp=(StartTime:(1/Channels_Fs(i)):StartTime+(length(eval(ChannelsList{i}))-1)*(1/Channels_Fs(i)))'; % This is the time vector associated with the 100Hz Flow data.
                Temp = interp1(TimeTemp,eval(ChannelsList{i}),Time,'linear');
                eval([ChannelsList{i} '= Temp;']);
                tol(i)=NaN;
%             else
%                 tol(i)=0.00000001;
%                 while 1
%                     [Q,P] = rat(Channels_Fs(i)/F_samp,tol(i));
%                     if (Q*P)>2^31
%                         tol(i)=tol(i)*10;
%                     else
%                         break
%                     end
%                 end
%                 %eval([ChannelsList{i} ' = resample(' ChannelsList{i} ',round(F_samp),round(Channels_Fs(i)));']); %only works if Fsamp / Fs are integer multiples of each other
%                 eval([ChannelsList{i} ' = resample(' ChannelsList{i} ',P,Q);']); %only works if Fsamp / Fs are integer multiples of each other
            end
            displaytext=['Resampling: ' ChannelsList{i} ' from ' num2str(Channels_Fs(i)) ' to ' num2str(F_samp) 'Hz using factor ' num2str(tol(i))];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            Channels_Fs(i)=F_samp;
        end
    end
    
    %Make channels the same length (sometimes these are off by a few samples)
    for i=1:length(ChannelsList)
        if ~isempty(ChannelsList{i})
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
    
    
    %These PSG signals are expected (hardcoded) in this order:
    %If signals are not available, we need to either (1) make a fake channel (above) or (2) edit the Analysis program to handle the absence
    
    ColumnHeadsList = {'Time','Thorax','Pnasal','Epochs','EventsResp','EventsAr','Abdomen','SpO2','EEG','Position'};
    
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
    
    %DataEventHypnog_Mat=[Flow_Aux_Time.' RIP_Thorax RIP_Abdo Flow_Aux SpO2 C3A2 Position CPAP];
    %ColumnHeads=[1       2            4       8        9             10         3           5       6           7          11      ];
    %           [1=Time  2=RIP_Thorax 3=Vflow 4=Hypnog 5=EventsResp  6=Arousals 7=RIP_Abdo 8=SpO2 9=EEG 10=Position 11=CPAP
    %e.g. Arousal=DataEventHypnog_Mat(:,ColumnHeads(6))=DataEventHypnog_Mat(:,10)
    %Position=DataEventHypnog_Mat(:,ColumnHeads(10))=DataEventHypnog_Mat(:,7)
    
    %ColumnHeads=[1      2            4      8        9         10            11         12       13         3           5       6       7];
    %            [1=Time 2=RIP_Thorax 3=Flow 4=Hypnog 5=Central 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]
    
    %ChannelsList = ChannelsList;%{'Time';'Thorax';'Abdomen';'Pnasal';'SaO2';'EEG';'Position';'Epochs';'EventsResp';'EventsAr';'CPAP';'Pes';'Edi';'GGpmax';'FlowPes';'FlowEdi';'alphaFlow';'kFlow'};
    
    %Channels={'Pnasal','Thorax','Abdomen','SaO2','EEG','Position','CPAP','Pes','Edi'};
    
    displaytext='Combining data into large matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat = zeros(length(Time),length(ChannelsList));
    for x=1:length(ChannelsList)
        eval(['DataEventHypnog_Mat(:,x)=' ChannelsList{x} ';']);
    end
    %ChannelsList(find(cellfun(@isempty,ChannelsList)))=[];
    
    if 0
        plotchannels=[ColumnHeads([3 4 5 6 8])];
        figure(1);
        for i=1:length(plotchannels)
            ax(i)=subplot(length(plotchannels),1,i);
            plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,plotchannels(i)))
        end
        linkaxes(ax,'x');
    end
    
    close all
    fclose all
    pack
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Spike data
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Spike(n,Filenames)
global handletext F_samp
try
    
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get data from Spike in Matlab format: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    filehandle = matfile([directory '\' fname]);
    w = whos('-file',[directory '\' fname]);
    
    %load([directory '\' fname]);
    
    %% Options hardcoded
    MultiplyPositionChannelby=5; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
    
    clear channelnameoptions
    channelnameoptions.Evts={'Evts','New_Evts'};
    channelnameoptions.Evts2={'Evts2','New_Evts2'};
    channelnameoptions.Epochs={'Epochs'};
    channelnameoptions.SaO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat','Sao2fing'};
    channelnameoptions.Position={'Position','Pos','pos','position'};
    %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal'};
    %channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
    channelnameoptions.EEG={'WakeSleep','EEG','EEG_C3_A2','C3_A2'};
    
    channelnameoptions.EEG1={'EEG1','EEG_C3_A2_clean','EEG_C3_A2','C3_A2'};
    channelnameoptions.EEG2={'EEG2','EEG_C4_A1_clean','EEG_C4_A1','C4_A1'};
    channelnameoptions.EEG3={'EEG3','EEG_F3_A2_clean','EEG_F3_A2','F3_A2'};
    channelnameoptions.EEG4={'EEG4','EEG_O2_A1_clean','EEG_O2_A1','O2_A1'};
    channelnameoptions.EKG={'EKG','ECG'};
    
    channelnameoptions.Thorax={'ThNoxRIP','Thorax','RC','Chest','CHEST','Belt2'};
    channelnameoptions.Abdomen={'AbNoxRIP','Abdomen','ABD','Abdom','ABDM','Belt1'};
    
    channelnameoptions.alphaFlow={'alphaFlow'};
    channelnameoptions.kFlow={'kFlow'};
    
    channelnameoptions.Pes={'Pes_clean','Pes','pes','Pepi','pepi'};
    channelnameoptions.Edi={'Edi'};
    channelnameoptions.FlowEdi={'FlowEdi'};
    channelnameoptions.FlowPes={'FlowPes'};
    channelnameoptions.GGpmax={'GGpmax','GGPmax'};
    %channelnameoptions.WakeSleep={'WakeSleep'};
    
    %channelnameoptions.EKG={'EKG','ECG'};
    if 0 %use pnasal signal
        channelnameoptions.Pnasal={'Pnasal','PNasal','Pmask','PMask'}; %'Vflow'
    else %use pneumotach flow
        channelnameoptions.Pnasal={'Vflow','Flow','PNasal','Pmask','PMask'}; %'Vflow'
    end
    channelnameoptions.CPAP={'CPAP','Pmask','PMask'};
    channelnamestemp=fieldnames(channelnameoptions);
    
    for i=1:length(channelnamestemp)
        temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
        foundamatch=0;
        for nn=1:length(temp)
            %Does it exist?
            for j=1:length(w)
                if strcmp(w(j).name,char(temp(nn)))
                    eval([channelnamestemp{i} '=filehandle.' char(temp(nn))]);
                    foundamatch=1;
                    break
                end
            end
            if foundamatch
                break
            end
        end
    end
    
    Channels={'Pnasal','Thorax','Abdomen','SaO2','EEG','Position','CPAP','Pes','Edi','GGpmax','FlowEdi','FlowPes','EEG1','EEG2','EEG3','EEG4','EKG','alphaFlow','kFlow'};
    optionalchannelsstartat=8;
    dt=1/F_samp;
    
    %% Resample and adjust [move EKG and EEG analysis to before the resampling here]
    for i=1:length(Channels)
        displaytext=['Collecting channel: ' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval(['Fs=(1/' Channels{i} '.interval);']);
            if round(Fs)~=F_samp
                displaytext=['Resampling: ' Channels{i} ' from ' num2str(round(Fs)) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} '.values = resample(' Channels{i} '.values,F_samp,round(Fs));']); %only works if Fsamp / Fs are integer multiples of each other
                eval([Channels{i} '.interval = 1/F_samp;']); %only works if Fsamp / Fs are integer multiples of each other
                eval([Channels{i} '.length = length(' Channels{i} '.values);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel: ' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            if i<optionalchannelsstartat %make empty channel with values = -1;
                eval([Channels{i} '.values = -1+ 0*Pnasal.values;']); %assumes we at least have a Pnasal channel
                eval([Channels{i} '.interval = 1/F_samp;']);
                eval([Channels{i} '.length = length(' Channels{i} '.values);']);
            else
                Channels{i}=[]; %make cell empty, later removed...
            end
        end
    end
    
    Channels(find(cellfun(@isempty,Channels)))=[];
    
    %% Get timing, epoch and event data
    
    %get time data
    timeXHz=(Pnasal.start:(1/F_samp):(Pnasal.start+(length(Pnasal.values)/F_samp)-1/F_samp))';
    N=length(timeXHz);
    
    %reserve memory for the large matrix
    DataEventHypnog_Mat=zeros(N,12);
    DataEventHypnog_Mat(:,1)=timeXHz;
    %clear timeXHz;
    
    %% New: fix channel lengths if not = N (not the same as Pnasal) -- incorporate this into Alice and others
    for i=1:length(Channels)
        if isempty(Channels{i})
            continue
        end
        eval(['ChannelN(i)=length(' Channels{i} '.values);']);
        if ChannelN(i)>N
            displaytext=['Length of ' Channels{i} ' channel is being altered from ' num2str(ChannelN(i)) ' to ' num2str(N)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values(N+1:end)=[];']);
            eval([Channels{i} '.length=N;']);
        end
        if ChannelN(i)<N
            displaytext=['Length of ' Channels{i} ' channel is being altered from ' num2str(ChannelN(i)) ' to ' num2str(N)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values((ChannelN(i)+1):N)=' Channels{i} '.values(end);']);
            eval([Channels{i} '.length=N;']);
        end
    end
    
    %% Epoch info
    if exist('Epochs','var')
        displaytext=['Get Hypnogram data:' fname];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
        tempdiff=round(double(Epochs.times));
        tempdiff2=diff(tempdiff);
        CountEpochsNot30s=sum(tempdiff2(1:end-1)~=30);
        if CountEpochsNot30s==0
            epochsXHz.values=interp1((Epochs.times+15)',double(Epochs.codes(:,1)'),DataEventHypnog_Mat(:,1),'nearest','extrap');
        else
            epochsXHz.values=zeros(1,N);
            for i=1:length(Epochs.times)
                if i<length(Epochs.times)
                    epochsXHz.values(DataEventHypnog_Mat(:,1)>=Epochs.times(i)&DataEventHypnog_Mat(:,1)<Epochs.times(i+1))=double(Epochs.codes(i,1));%(Epochs.codes(i,1));
                elseif i==length(Epochs.times)
                    %epochsXHz.values(timeXHz>=Epochs.times(i))=double(Epochs.codes(i));
                end
            end
        end
        
        recode=1
        if recode
            %Spike code: W=0,N1=1,N2=2,N3=3,R=5,?=8
            %recode Epochs to match Phil's naming: N3=0,N2=1,N1=2,R=3,W=4
            Epochs.values=3-epochsXHz.values;
            Epochs.values(epochsXHz.values==0)=4;
            Epochs.values(epochsXHz.values==5)=3;
            Epochs.values(epochsXHz.values==8)=NaN;
        end
        clear epochsXHz;
    end
    
    
    %% Events info
    
    %%%to do: disp error if no arousals / handle error with first row of zeros %%%
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    if exist('Evts','var'), New_Evts=Evts; end
    if exist('Evts2','var'), NewEvts2=Evts2; end
    if exist('New_Evts2','var'), NewEvts2=New_Evts2; end
    
    while New_Evts.codes(1,1)==0
        New_Evts.codes(1,:)=[];
        New_Evts.times(1)=[];
        New_Evts.text(1,:)=[];
        New_Evts.length=New_Evts.length-1;
        'removed_one'
    end
    if exist('NewEvts2','var')
        while NewEvts2.codes(1,1)==0
            NewEvts2.codes(1,:)=[];
            NewEvts2.times(1)=[];
            NewEvts2.text(1,:)=[];
            NewEvts2.length=NewEvts2.length-1;
            'removed_one 2'
        end
        if exist('NewEvts2','var')&&exist('New_Evts','var')
            %combine events lists
            New_Evts_backup1=New_Evts;
            New_Evts.times=[New_Evts.times;NewEvts2.times];
            New_Evts.codes=[New_Evts.codes;NewEvts2.codes];
            New_Evts.text=[New_Evts.text;NewEvts2.text];
            New_Evts.length=New_Evts.length+NewEvts2.length;
        end
    end
    
    ii=1; %loop cycles multiple times, removing any double rows of zeros.
    while ii<length(New_Evts.codes)
        ii=1;
        for i=2:length(New_Evts.times)
            evennumber=1-mod(i,2);
            oddnumber=1-evennumber;
            if 1
                if oddnumber&&New_Evts.codes(i)==0&&New_Evts.codes(i-1)==0
                    New_Evts.codes(i,:)=[];
                    New_Evts.times(i)=[];
                    New_Evts.text(i,:)=[];
                    New_Evts.length=New_Evts.length-1;
                    %'removed_one double row of zeros'
                    break
                end
            end
            %     if evennumber&New_Evts.codes(i)>0
            %         'saysomething'
            %         break
            %     end
            ii=ii+1;
        end
    end
    
    
    %add a row of zeros if there is a missing row of zeros. Assume end of event is the start of the next event (or end of file/Inf).
    for i=length(New_Evts.codes):-1:1
        if i==length(New_Evts.codes)
            if New_Evts.codes(i)>0
                New_Evts.codes=[New_Evts.codes;[0 0 0 0]];
                New_Evts.times=[New_Evts.times;Inf];
                New_Evts.text=[New_Evts.text;['                                            ']];
                New_Evts.length=New_Evts.length+1;
                disp('added_row_of_zeros to the end');
            end
        end
        if i>1
            if New_Evts.codes(i)>0&&New_Evts.codes(i-1)>0
                New_Evts.codes=[New_Evts.codes(1:(i-1),:);[0 0 0 0];New_Evts.codes(i:end,:)];
                New_Evts.times=[New_Evts.times(1:i);New_Evts.times(i:end)];
                New_Evts.text=[New_Evts.text(1:(i-1),:);['                                            '];New_Evts.text(i:end,:)];
                New_Evts.length=New_Evts.length+1;
                disp('added_row_of_zeros');
            end
        end
    end
    
    %% Make list of event types
    
    clear EventTypeList_text EventTypeList_code tempi
    EventTypeList_text=[];
    for i=1:2:size(New_Evts.codes,1)
        codematch=0;
        for j=1:size(EventTypeList_text,1)
            if strcmp(char(EventTypeList_text(j,:)),New_Evts.text(i,:))
                %'found code match'
                codematch=1;
                break
            end
        end
        if codematch==0
            tempi=size(EventTypeList_text,1);
            EventTypeList_code(tempi+1)=New_Evts.codes(i);
            EventTypeList_text(tempi+1,:)=New_Evts.text(i,:);
        end
    end
    EventTypeList_code=EventTypeList_code'
    char(EventTypeList_text)
    
    %Make list fit the following categories:
    %1 AR
    %2 ApO | Ap-O
    %3 ApC | Ap-C
    %4 HypO | H (Y)
    %5 M (Y)
    %6 HypC
    %7 HypO2% (Y)
    %8 HypC2%
    %9 HypOx (Y)
    %10 HypCx
    
    New_Evts_backup2=New_Evts;
    
    for i=1:size(New_Evts.codes,1)
        if strncmp(New_Evts.text(i,:),'AR',2)
            New_Evts.codes(i,1)=1;
        elseif strncmp(New_Evts.text(i,:),'ApO',3)
            New_Evts.codes(i,1)=2;
        elseif strncmp(New_Evts.text(i,:),'Ap-O',4)
            New_Evts.codes(i,1)=2;
        elseif strncmp(New_Evts.text(i,:),'ap-O',4)
            New_Evts.codes(i,1)=2;
        elseif strncmp(New_Evts.text(i,:),'ApC',3)
            New_Evts.codes(i,1)=3;
        elseif strncmp(New_Evts.text(i,:),'Ap-C',4)
            New_Evts.codes(i,1)=3;
        elseif strncmp(New_Evts.text(i,:),'M',1)
            New_Evts.codes(i,1)=5;
        elseif strncmp(New_Evts.text(i,:),'HypO2',5)
            New_Evts.codes(i,1)=7;
        elseif strncmp(New_Evts.text(i,:),'HypC2',5)
            New_Evts.codes(i,1)=8;
        elseif strncmp(New_Evts.text(i,:),'HypOx',5)
            New_Evts.codes(i,1)=9;
        elseif strncmp(New_Evts.text(i,:),'HypCx',5)
            New_Evts.codes(i,1)=10;
        elseif strncmp(New_Evts.text(i,:),'HypO',4)
            New_Evts.codes(i,1)=4;
        elseif strncmp(New_Evts.text(i,:),'HypC',4)
            New_Evts.codes(i,1)=6;
        elseif strncmp(New_Evts.text(i,:),'Hyp-C',5)
            New_Evts.codes(i,1)=6;
        elseif strncmp(New_Evts.text(i,:),'H',1)
            New_Evts.codes(i,1)=4;
        end
    end
    
    %%%Add this in%%%
    % for x=1:2:length(New_Evts.times)
    %     New_Evts.duration((x-1)/2+1)=New_Evts.times(x+1)-New_Evts.times(x);
    % end
    
    %% Make events in continuous time
    % if exist('New_Evts','var')
    
    %channel_list = {... 'EventsCentralApnea';'EventsObstrApnea';'EventsObstrHyp';'Desats';'EventsAr';'EventsMixedApnea';'EventsCentralHyp'};
    EventsAr.values=zeros(N,1); EventsAr.interval=dt;
    EventsResp.values=zeros(N,1); EventsResp.interval=dt;
    %     EventsObstrApnea.values=zeros(N,1); EventsObstrApnea.interval=dt;
    %     EventsCentralApnea.values=zeros(N,1); EventsCentralApnea.interval=dt;
    %     EventsObstrHyp.values=zeros(N,1); EventsObstrHyp.interval=dt;
    %     EventsMixedApnea.values=zeros(N,1); EventsMixedApnea.interval=dt;
    %     EventsCentralHyp.values=zeros(N,1); EventsCentralHyp.interval=dt;
    %
    for x=1:2:length(New_Evts.times)
        lefti=round((New_Evts.times(x)-DataEventHypnog_Mat(1,1))*F_samp)+1;
        righti=lefti+round((New_Evts.times(x+1)-New_Evts.times(x))*F_samp);
        if righti>N, righti=N; end
        if New_Evts.codes(x)==1
            EventsAr.values(lefti:righti)=New_Evts.codes(x);
        end
        if New_Evts.codes(x)>1
            EventsResp.values(lefti:righti)=New_Evts.codes(x);
        end
    end
    
    % Plot events
    if 0
        figure(1);
        plot(DataEventHypnog_Mat(:,1),EventsAr.values)%,New_Evts.times,New_Evts.codes(:,1),'.'); ylabel('EventsIndex'); end
        set(gca,'FontName','Arial Narrow','FontSize',10,'box','off','YTickLabel',{'AR','OA','CA','OH','M','CH','OH2','CHx','OHx','CHx'},'YTick',[1 2 3 4 5 6 7 8 9 10],'XTick',[],'Xcolor',[1 1 1]);
        ylabel('Events/Arousal');
    end
    
    %% Position data at XHz.
    
    clear Mexception
    try %if channel doesn't exist: program will continue
        channel_name = 'Position';
        Position.values=Position.values*MultiplyPositionChannelby;
        displaytext='Scaled position data';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    catch Mexception
        displaytext='No position data, looking for file';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    % Fix Pos data if necessary
    directorytemp = directory;     if directorytemp(end)=='\',directorytemp(end)=[];end
    textfilename=[directorytemp '\' fname(1:(length(fname)-4)) '_pos.txt'];
    if exist(textfilename,'file')==2
        displaytext='Found position file';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
        %PosXHz_original=Position.values;
        [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
        for i=1:(length(col1)-1)
            Position.values(DataEventHypnog_Mat(:,1)>=col1(i)&DataEventHypnog_Mat(:,1)<col1(i+1))=col2(i);
        end
        Position.values(DataEventHypnog_Mat(:,1)>=col1(end))=col2(end);
    end
    %channel_exist=sum(abs(Position.values))>0;
    
    if 0
        figure(1), plot(DataEventHypnog_Mat(:,1),Position.values);
        set(gca,'FontName','Arial Narrow','FontSize',10,'box','off','YTickLabel',{'Up','Left','Supine','Prone','Right'},'YTick',[0 1 2 3 4],'XTick',[],'Xcolor',[1 1 1],'Ylim',[0 5]);
        ylabel('Position');
    end
    %clear PosXHz_original
    
    %% Remove artifact in all signals using text files
    for j=1:length(Channels)
        if isempty(Channels{j})
            continue
        end
        directorytemp = directory;     if directorytemp(end)=='\',directorytemp(end)=[];end
        textfilename=[directorytemp '\' fname(1:end-4) '_' Channels{j} '_art.txt'];
        if exist(textfilename,'file')==2
            ['found ' Channels{j} ' artifact file and removed artifact']
            %Position.values_original=Position.values;
            [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
            for i=1:length(col1)
                dttemp=eval([Channels{j} '.interval']);
                try
                    lefti=round((col1(i)-eval([Channels{j} '.start']))/dttemp +1);
                catch me
                    lefti=round((col1(i)-0)/dttemp +1);
                end
                if lefti<1, lefti=1; end
                try
                    righti=round((col2(i)-eval([Channels{j} '.start']))/dttemp +1);
                catch me
                    righti=round((col2(i)-0)/dttemp +1);
                end
                if righti>length(eval([Channels{j} '.values'])), righti=length(eval([Channels{j} '.values'])); end
                eval([Channels{j} '.values(lefti:righti)=NaN;']);
            end
            Percent_removed(j) = sum(isnan(eval([Channels{j} '.values'])))/length(eval([Channels{j} '.values']))*100
            %figure(5); plot(Time,eval([channelnamestemp{j} '.values']));
        end
    end
    
    %% BetaPower--Find best EEG
    EEGsignals = {'EEG1','EEG2','EEG3','EEG4'};
    temp=[];
    for i=1:length(EEGsignals)
        if exist(EEGsignals{i},'var')
            temp{length(temp)+1}=EEGsignals{i};
        end
    end
    EEGsignals=temp;
    
    [ECG_peak_i,HR] = EKGpeakdetection(EKG.values,timeXHz,dt);
    
    
    BetaPowerSettings.polyorder=1;
    BetaPowerSettings.suf1 = '.values';
    BetaPowerSettings.timestr = 'DataEventHypnog_Mat(:,1)';
    BetaPowerSettings.timeeegstr = 'DataEventHypnog_Mat(:,1)';
    BetaPowerSettings.flowstr = 'Pnasal';
    BetaPowerSettings.arstr = 'EventsAr';
    BetaPowerSettings.epochsstr = 'Epochs';
    BetaPowerSettings.spo2str = 'SpO2';
    
    
    fftlengthoptions = 2.^[6:9];
    [~,bestoptioni] = min(abs(1/dt - fftlengthoptions));
    
    
    BetaPowerSettings.fft_length = fftlengthoptions(bestoptioni); %256 %Best: 1710,1723:256, 1343:512
    BetaPowerSettings.scoredarousalsinwake = 1; %1 for Lauren's scoring.
    %[WakeSleep,WakeSleepInfo,EEGsignals,EEGsignalsOut,Powersignals,PowersignalsOut] = BetaPowerRun(EEGsignals,dt,BetaPowerSettings);
    [WakeSleep,WakeSleepInfo,EEGsignals,EEGsignalsOut,Powersignals,PowersignalsOut] = BetaPowerRun(EEGsignals,dt,BetaPowerSettings,ECG_peak_i);
    
    
    
    
    %Clean EEG signals
    for i=1:length(EEGsignals)
        eval([EEGsignals{i} '.values=EEGsignalsOut{i};'])
    end
    
    %best EEG
    [~,besti]=max(WakeSleepInfo.AUC_M);
    EEG.values = eval([EEGsignals{besti} '.values']);
    
    %WakeSleep
    temp = WakeSleep;
    clear WakeSleep;
    WakeSleep.values = temp;
    
    %Powersignals
    for i=1:length(Powersignals)
        eval([Powersignals{i} '.values=PowersignalsOut{i};'])
    end
    
    %% Channel data at XHz: Combining data into large matrix
    %DataEventHypnog_Mat=[Flow_Aux_Time.' RIP_Thorax RIP_Abdo Flow_Aux SpO2 C3A2 Position CPAP];
    ColumnHeads=[1       2            4       8        9             10         3           5       6           7          11      ];
    %           [1=Time  2=RIP_Thorax 3=Vflow 4=Hypnog 5=EventsResp  6=Arousals 7=RIP_Abdo 8=SpO2 9=EEG 10=Position 11=CPAP
    %Arousal=DataEventHypnog_Mat(:,ColumnHeads(6))=DataEventHypnog_Mat(:,10)
    %Position=DataEventHypnog_Mat(:,ColumnHeads(10))=DataEventHypnog_Mat(:,7)
    
    %ColumnHeads=[1      2            4      8        9         10            11         12       13         3           5       6       7];
    %            [1=Time 2=RIP_Thorax 3=Flow 4=Hypnog 5=Central 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]
    
    ChannelsList = {'Time';'Thorax';'Abdomen';'Pnasal';'SaO2';'EEG';'Position';'Epochs';'EventsResp';'EventsAr';'CPAP';'Pes';'Edi';'GGpmax';'FlowPes';'FlowEdi';'alphaFlow';'kFlow'};
    ChannelsList = [ChannelsList;'WakeSleep';Powersignals'];
    %Channels={'Pnasal','Thorax','Abdomen','SaO2','EEG','Position','CPAP','Pes','Edi'};
    
    displaytext='Combining data into large matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    optionalchannelsstartat=12;
    index=1;
    for x=2:length(ChannelsList)
        if x<optionalchannelsstartat
            index=x;
        end
        try %if channel doesn't exist: program will continue but an error message will show
            if x<optionalchannelsstartat
                eval(['DataEventHypnog_Mat(:,x)=' ChannelsList{x} '.values;']);
            else
                if exist(ChannelsList{x},'var')
                    index=index+1;
                    displaytext=['Optional channel found: ' ChannelsList{x}];
                    disp(displaytext); set(handletext,'String',displaytext); drawnow;
                    ColumnHeads(index)=index;
                    eval(['DataEventHypnog_Mat(:,index)=' ChannelsList{x} '.values;']);
                else
                    displaytext=['No optional channel found: ' ChannelsList{x}];
                    disp(displaytext); set(handletext,'String',displaytext); drawnow;
                    ChannelsList{x}=[];
                end
            end
        catch me
            displaytext=['No channel found: ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
        end
    end
    ChannelsList(find(cellfun(@isempty,ChannelsList)))=[];
    
    if 0
        plotchannels=[2 3 4 6 20]
        figure(1);
        for i=1:length(plotchannels)
            ax(i)=subplot(length(plotchannels),1,i);
            plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,plotchannels(i)))
        end
        linkaxes(ax,'x');
    end
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Alice Sleepware G3
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = AliceSleepwareG3(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get rml data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    S=xml2struct([directory '\' fname]);
    
    displaytext=['Get recording start time from rml: ' directory];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    StartText = S.PatientStudy.Acquisition.Sessions.Session.RecordingStart.Text(12:19);
    StartTime = mod(datenum(StartText,'HH:MM:SS'),1)*86400
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the .csv files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    NChannels = length(Channels);
    %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
    %EDF Start time -- have to import codes
    directory = char(Filenames(n,4));
    
    fname=char(Filenames(n,1));
    EDFfilenamedir = [directory '\' fname];
    
    %M = csvread(filename)
    if 1
        for i=1:length(Channels)
            displaytext=['Collecting channel:' Channels{i}];
            try
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                
                eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
                %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
                if Fs~=F_samp
                    displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                    disp(displaytext); set(handletext,'String',displaytext); drawnow;
                    eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
                end
            catch me
                displaytext=['No channel:' Channels{i} ' ' me.message];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
            end
        end
    end
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time;
    %clear Pnasal_Time; %save memory
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
        eval(['clear ' Channels{i}]); %save memory
    end
    
    %Xml Hypnogram
    displaytext=['Get Hypnogram data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    
    StageTextOptions = {'Wake','NonREM1','NonREM2','NonREM3','REM','NotScored'};
    StageTextCodes = [4,2,1,0,3,-1];
    Ntransitions = length(S.PatientStudy.ScoringData.StagingData.UserStaging.NeuroAdultAASMStaging.Stage);
    clear Transitions
    for i=1:Ntransitions
        Transitions{i}.Text = S.PatientStudy.ScoringData.StagingData.UserStaging.NeuroAdultAASMStaging.Stage{i}.Attributes.Type;
        Transitions{i}.Time = str2num(S.PatientStudy.ScoringData.StagingData.UserStaging.NeuroAdultAASMStaging.Stage{i}.Attributes.Start)/30 +1;
        Transitions{i}.Code = -1; %Default is unknown
    end
    for i=1:Ntransitions
        for m=1:length(StageTextOptions)
            if strcmp(StageTextOptions{m},Transitions{i}.Text)==1
                %[StageTextOptions{m} ' ' num2str(StageTextCodes(m))]
                Transitions{i}.Code = StageTextCodes(m);
            end
        end
    end
    
    RecordingDuration = str2num(S.PatientStudy.Acquisition.Sessions.Session.Duration.Text);
    Nepochs = ceil(RecordingDuration/30);
    
    clear hypnogram
    hypnogram=NaN*zeros(Nepochs,1);
    for i=1:Ntransitions
        indexi = Transitions{i}.Time;
        hypnogram(indexi:end) = Transitions{i}.Code;
    end
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz;
    clear EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Xml/rml events [tested and working]
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    
    TempEvents = S.PatientStudy.ScoringData.Events.Event;
    EventTypesOfInterest = {'ObstructiveApnea','Hypopnea','CentralApnea','MixedApnea','Arousal','RelativeDesaturation'};
    
    Nevents = length(S.PatientStudy.ScoringData.Events.Event);
    clear Events
    for i=1:Nevents
        Events{i}.name = S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Type;
        Events{i}.start = StartTime + str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Start);
        Events{i}.duration = str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Duration);
    end %W=0,N1=1,N2=2,N3=3,R=5
    
    %Add Events data in additional colums of the matrix:
    %**************************************************************************
    EventCategories={'CentralApnea','ObstructiveApnea','MixedApnea','Hypopnea','RelativeDesaturation','Arousal'}
    
    EventCols=[10 11 15 12 13 14]; %Mixed apneas are no longer treated as Obstructive apneas here.
    
    for m=1:Nevents
        for i=1:length(EventCategories)
            if strcmp(Events{m}.name,EventCategories{i})==1
                %Events{m}.name
                lefti=round((Events{m}.start-StartTime)*F_samp)+1;
                righti=lefti+round((Events{m}.duration)*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
            end
        end
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','�-arousal'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
    
    
    if 1 %check plot
        figure(1)
        plotcols=[2,11,12,3,9,7];
        downsamplefactor=4;
        hold('off')
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
            ylabel(ChannelsList{plotcols(i)});
        end
        linkaxes(ax,'x');
    end
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% GrassTwin [needs checking]
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = GrassTwin(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get xls data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    [~,xlstxt,~] = xlsread([directory '\' fname],1,'C:D')
    
    
    displaytext=['Get recording start time from xlstxt'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    StartText = xlstxt{1,1};
    StartTime = mod(datenum(StartText,'HH:MM:SS'),1)*86400
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the .csv files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    NChannels = length(Channels);
    %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
    %EDF Start time -- have to import codes
    directory = char(Filenames(n,4));
    
    fname=char(Filenames(n,1));
    EDFfilenamedir = [directory '\' fname];
    
    %M = csvread(filename)
    if 1
        for i=1:length(Channels)
            displaytext=['Collecting channel:' Channels{i}];
            try
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                
                eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
                %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
                if Fs~=F_samp
                    displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                    disp(displaytext); set(handletext,'String',displaytext); drawnow;
                    eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
                end
            catch me
                displaytext=['No channel:' Channels{i} ' ' me.message];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
            end
        end
    end
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time;
    %clear Pnasal_Time; %save memory
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
        eval(['clear ' Channels{i}]); %save memory
    end
    
    %Xls Hypnogram
    displaytext=['Get Hypnogram data from xlstxt'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    
    StageTextOptions = {'Stage - W','Stage - R','Stage - N1','Stage - N2','Stage - N3','Stage - No Stage'};
    StageTextCodes = [4,3,2,1,0,-1]; %Terrill code: W=4,R=3,N1=2,N2=1,N3=0 [unknown=-1]
    
    xlstxthyponly = xlstxt;
    temp = 'Stage - ';
    for i=size(xlstxthyponly,1):-1:1
        temphyp=xlstxthyponly{i,2};
        if length(temphyp)<8||~strcmp(temphyp(1:8),'Stage - ')
            xlstxthyponly(i,:)=[];
        end
    end
    Ntransitions = size(xlstxthyponly,1);
    LastEpochStartTime = mod(datenum(xlstxthyponly{Ntransitions,1},'HH:MM:SS'),1)*86400
    if LastEpochStartTime<43200, LastEpochStartTime=LastEpochStartTime+86400; end
    Nepochs = round((LastEpochStartTime-StartTime)/30+1);
    
    clear HypTimes
    for i=1:Ntransitions
        HypTimes(i,:) = mod(datenum(xlstxthyponly{i,1},'HH:MM:SS'),1)*86400
        if HypTimes(i)<43200, HypTimes(i)=HypTimes(i)+86400; end
    end
    
    clear HypCodes
    for i=1:Ntransitions
        for m=1:length(StageTextOptions)
            if strcmp(xlstxthyponly{i,2},StageTextOptions{m})
                HypCodes(i,:) = StageTextCodes(m);
            end
        end
    end
    
    HypEpochIndex = round((HypTimes-StartTime)/30+1);
    %RecordingDuration = ... **might be more robust to use the EDF data for this in future
    %Nepochs = ceil(RecordingDuration/30);
    
    clear hypnogram
    hypnogram=NaN*zeros(Nepochs,1);
    for i=1:Ntransitions %write from current epoch to end file with new stage code
        indexi = round(HypEpochIndex(i));
        hypnogram(indexi:end) = HypCodes(i);
    end
    
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz;
    clear EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Xls events
    displaytext=['Get Events data from xlstxt'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    
    %xlstxt
    EventTypesOfInterest = {'Central Apnea','Obstructive Apnea','Mixed Apnea','Hypopnea','RERA','Desaturation','Arousal'};
    EventCodes = [10 11 15 12 12 13 14]; %RERAs are treated as hypopneas here
    %Nevents = length(S.PatientStudy.ScoringData.Events.Event);
    index = 1;
    clear Events
    for i=1:size(xlstxt,1)
        for m=1:length(EventTypesOfInterest)
            if ~isempty(findstr(xlstxt{i,2},EventTypesOfInterest{m}))
                Events.name{index}=EventTypesOfInterest{m};
                Events.code(index)=EventCodes(m);
                Events.row(index)=i;
                timetemp = mod(datenum(xlstxt{i,1},'HH:MM:SS'),1)*86400;
                if timetemp<43200, timetemp=timetemp+86400; end
                Events.start(index)=timetemp;
                tempk1 = findstr(xlstxt{i,2},'Dur:') + 5;
                tempk2 = findstr(xlstxt{i,2},'sec') - 2;
                Events.duration(index) = str2num(xlstxt{i,2}(tempk1:tempk2));
                index=index+1;
            end
            continue
        end
    end
    Nevents = length(Events.code);
    
    %Add Events data in additional colums of the matrix:
    for i=1:Nevents
        lefti=round((Events.start(i)-StartTime)*F_samp)+1;
        righti=lefti+round((Events.duration(i))*F_samp);
        DataEventHypnog_Mat(lefti:righti,Events.code(i))=1;
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','�-arousal'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
    
    
    if 1 %check plot
        figure(1)
        plotcols=[2,11,12,3,9,7];
        downsamplefactor=4;
        hold('off')
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
            ylabel(ChannelsList{plotcols(i)});
        end
        linkaxes(ax,'x');
    end
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Remlogic1p1
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Remlogic1p1(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    %EDF Start time
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get recording start time from EDF: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    fid = fopen([directory '\' fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the EDF files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    %F_samp=100;
    
    EDFfilenamedir = [directory '\' fname];
    
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time;
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
    end
    
    %Hypnogram
    directory = char(Filenames(n,6));
    fname=char(Filenames(n,3));
    displaytext=['Get Hypnogram data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    clear fileinfo
    fid = fopen([directory '\' fname]);
    i=1;
    while 1
        fileinfo{i,:} = fgetl(fid);  %read in the data
        if fileinfo{i,:}==-1
            fileinfo(i,:) = [];
            break
        end
        i=i+1;
    end
    fclose(fid);   %close the file
    
    starttext = 'Sleep Stage';
    %%%up to here
    for i=1:length(fileinfo)
        if ~isempty(findstr(fileinfo{i,:},starttext))
            break
        end
    end
    startline = i+1;
    for i=1:length(fileinfo)
        if ~isempty(findstr(fileinfo{i,:},starttext))
            break
        end
    end
    sleeptable = fileinfo(startline:end);
    
    clear SleepCodesTxt
    for i=1:length(sleeptable);
        SleepCodesTxt{i} = sleeptable{i}(1:2);
    end
    
    StageTextOptions = {'W','R','N1','N2','N3'};
    StageTextCodes = [4,3,2,1,0];
    
    clear hypnogram
    for i=1:length(SleepCodesTxt)
        for m=1:length(StageTextOptions)
            if findstr(SleepCodesTxt{i},StageTextOptions{m})
                hypnogram(i)=StageTextCodes(m);
                break
            end
        end
    end
    
    HypStartTimeTxt = sleeptable{1}(11:18);
    
    timestrrangetemp=findstr(sleeptable{1},':');
    timestrrange = timestrrangetemp(1)-2:timestrrangetemp(2)+2;
    HypStartTimeTxt = sleeptable{1}(timestrrange);
    
    HypStartTime = mod(datenum(HypStartTimeTxt,'HH:MM:SS'),1)*86400;
    if HypStartTime<43200, HypStartTime=HypStartTime+86400; end
    
    hypnogram_t = HypStartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Events
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    Events_fname = [directory '\' fname];
    
    clear fileinfo
    fid = fopen(Events_fname);
    i=1;
    while 1
        fileinfo{i,:} = fgetl(fid);  %read in the data
        if fileinfo{i,:}==-1
            fileinfo(i,:) = [];
            break
        end
        i=i+1;
    end
    fclose(fid);   %close the file
    
    starttext = 'Sleep Stage';
    for i=1:length(fileinfo)
        if ~isempty(findstr(fileinfo{i,:},starttext))
            break
        end
    end
    startline = i+1;
    for i=1:length(fileinfo)
        if ~isempty(findstr(fileinfo{i,:},starttext))
            break
        end
    end
    EventsTable = fileinfo(startline:end);
    
    
    %up to here
    clear Event_StartTime Event_StartTimeTxt
    for i=1:length(EventsTable)
        timestrrangetemp=findstr(EventsTable{i},':');
        timestrrange = timestrrangetemp(1)-2:timestrrangetemp(2)+2;
        Event_StartTimeTxt{i} = EventsTable{i}(timestrrange);
        Event_StartTime(i,:) = mod(datenum(Event_StartTimeTxt{i},'HH:MM:SS'),1)*86400;
        if Event_StartTime(i)<43200, Event_StartTime(i)=Event_StartTime(i)+86400; end
    end
    
    clear Event_DurationText Event_Duration
    for i=1:length(EventsTable)
        Event_DurationText{i}=EventsTable{i}(end-1:end);
        Event_Duration(i,:)=str2num(Event_DurationText{i});
    end
    
    clear Event_StartTime Event_StartTimeTxt
    for i=1:length(EventsTable)
        timestrrangetemp=findstr(EventsTable{i},':');
        timestrrange = timestrrangetemp(1)-2:timestrrangetemp(2)+2;
        Event_StartTimeTxt{i} = EventsTable{i}(timestrrange);
        Event_StartTime(i,:) = mod(datenum(Event_StartTimeTxt{i},'HH:MM:SS'),1)*86400;
        if Event_StartTime(i)<43200, Event_StartTime(i)=Event_StartTime(i)+86400; end
    end
    %up to here
    EventCategories={'APNEA-CENTRAL','APNEA-OBSTRUCTIVE','APNEA-MIXED','HYPOPNEA','AROUSAL'}; %note there is no desat column but an error will occur if it is missing
    EventCols=[10 11 15 12 14];
    
    %initialize
    DataEventHypnog_Mat(1:size(DataEventHypnog_Mat,1),10:15)=0;
    
    clear Event_Type
    for m=1:length(EventsTable)
        for i=length(EventCategories):-1:1
            temp = findstr(EventsTable{m},EventCategories{i});
            if ~isempty(temp)
                Event_Type{m,:}=EventsTable{m}(temp:temp+length(EventCategories{i})-1);
                if EventCols(i)==14&&Event_Duration(m)<3
                    Event_Duration(m)==3; %force arousal duration to be at least 3 s
                end
                lefti=round((Event_StartTime(m)-StartTime)*F_samp)+1;
                righti=lefti+round(Event_Duration(m)*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
                break
            end
        end
    end
    
    if 0 %check plot
        figure(1)
        plotcols=[2,11,12,14,3,9];
        downsamplefactor=5;
        hold('off')
        count=1;
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
        end
        linkaxes(ax,'x');
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','�-arousal'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end
