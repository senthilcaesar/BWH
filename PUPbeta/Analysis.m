function Analysis(n1,settings1,AMasterSpreadsheet1)
%% global variables and settings
global ChannelsList ChannelsFs Evts settings AMasterSpreadsheet n Info

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

n=n1;
ChannelsList = []; % ChannelsList is new global from here, so start with clean slate

try
    %% load file
    settings.filename=[char(settings.patients(n,1))]; %seems to be unused.
    % DLM says, settings.filename is used in LGfromFlowBetaPart1 to
    % save flowdrive plot, and Part2 to save loop gain plot
%     directoryn=char(settings.patients(n,2));
%     if directoryn(end)~=filesep
%         directoryn=[directoryn filesep];
%     end
    MATfilename=[settings.ConvertedDirectory char(settings.patients(n,1))];
    
    Evts=struct();
    temp=[];
    
    disp(['Loading File:' MATfilename]);
    
    temp=load(MATfilename);
    if isfield(temp, 'DataEventHypnog_Mat') % convert DataEventHypnog_Mat to SigT
        SigT=array2table(temp.DataEventHypnog_Mat);
        SigT.Properties.VariableNames = temp.ChannelsList;
        temp = rmfield(temp, {'DataEventHypnog_Mat','ChannelsList'});
    else
        SigT=temp.SigT;
    end
    
    
    Evts=temp.Evts;
    %WakeSleepInfo=temp.WakeSleepInfo; % removed after modifying 'Info' in Convert
    ChannelsFs=temp.ChannelsFs;
    ChannelsList=SigT.Properties.VariableNames;
    try Info=temp.Info;
    catch
        Info=temp.WakeSleepInfo;
    end
    
    % using automated flow upright detector
    if (isfield(settings,'PnasalUprightAuto') && settings.PnasalUprightAuto==1) || ~isfield(settings,'Pnasaldownisinsp')
        %Attempt AutoUpright
        try
            successtemp=0;
            disp('Warning: Using Automated Flow Upright Detection from Converted File');
            if exist('Info') && isfield(Info,'FlowQ') && isfield(Info.FlowQ,'PrUpright')
                if isnan(Info.FlowQ.PrUpright)
                    disp('Info.FlowQ.PrUpright is not a number, failed');
                else
                    settings.Pnasaldownisinsp = Info.FlowQ.PrUpright<0.5;
                    disp(['Success: Info.FlowQ.PrUpright = ' num2str(Info.FlowQ.PrUpright)]);
                    successtemp=1;
                end
            else
                disp('Info.FlowQ.PrUpright not found');
            end
            
            if successtemp==0
                disp('Calculating Info.FlowQ.PrUpright')
                [Info.FlowQ.PrUpright,Info.FlowQ.FnoiseAll]=FlowInvertedDetectTool(SigT.Flow,SigT.Time);
                if isnan(Info.FlowQ.PrUpright)
                    disp('Info.FlowQ.PrUpright is not a number, failed')
                else
                    settings.Pnasaldownisinsp = Info.FlowQ.PrUpright<0.5;
                    disp(['Success: Info.FlowQ.PrUpright = ' num2str(Info.FlowQ.PrUpright)]);
                end
            end
        catch
            disp('Automated Flow Upright Detection failed unexpectedly');
        end
    else
        disp(['Pnasaldownisinsp: ', num2str(settings.Pnasaldownisinsp)]);
        disp(['invertflow: ', num2str(settings.invertflowlist(n1))]);
        
    end
    
    
    %% Set this save location up now, save time later
    % will end up saving within indiv pt folders within this location
    
    % DLM added AnalyzedDirectory, because it was called but not set up
    % it does see a little redundant, and a bit dangerous to have two
    % global variables initialized to same value...
    % DV: if you do NOT isfield, then in the next line it looks for the
    % field in settings and has an error, I removed just so I could run
    % my code.
    if 0
        if ~isfield(settings,'RunPUPA2Z') && ~settings.RunPUPA2Z % causing some weird issues in A to Z mode
            savedir=[settings.AnalyzedDirectory settings.savename];
            
            if ~(exist(savedir, 'dir') == 7)
                mkdir(savedir);
            end
            
        end
    end
    
    settings.OutputDataDirectory = settings.AnalyzedDirectory; %default
    if 1 % this is the original code
        savedir=[settings.OutputDataDirectory settings.savename];
        if ~(exist(savedir, 'dir') == 7)
            mkdir(savedir);
        end
    else % this is just applied here to keep things working, not ideal!
        settings.OutputDataDirectory = settings.AnalyzedDirectory;
    end
    
    %% Set position codes according to protocol
    settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
    try
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
    catch
        % above code failing if settings.protocol is char array.
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(n,:));
    end
    
    settings.supinepositioncode = settings.positioncodes(1);
    
    if settings.plotfigure
        figure(101); clf(101); fig101 = gcf;
        fig101.Color = [1 1 1]; % set background colour to white
        fig101.Units = 'pixels';
        MonitorPositions = get(0,'MonitorPositions');
        % just use a portion of the main display
        if ~isempty(MonitorPositions)
            win_x = MonitorPositions(1,3)*0.3;
            win_y = MonitorPositions(1,4)*0.6;
            offset_x = MonitorPositions(1,3)*0.65;
            offset_y = MonitorPositions(1,4)*0.4;
            fig101.Position = [offset_x offset_y win_x win_y];
        end
    end
    
    %% get CPAP data
    if ~settings.ignoreCPAPdata
        [CPAPoff,CPAP]=getCPAP(SigT,settings,ChannelsList);
    else %a priori knowledge that CPAP is not administered
        CPAPoff=ones(height(SigT),1);
        CPAP=zeros(height(SigT),1);
    end
    
    %% get AHI and events data
    if ~isfield(Evts,'AHIdata') && ...
            (~isfield(settings, 'SiteOfCollapse') || (isfield(settings, 'SiteOfCollapse') && settings.SiteOfCollapse == 0))
        disp('No AHI Data- Run getEvtsAddOn.m to get Evts Variable');
    end
    
    %% get SpO2 data
    if ~isfield(Evts,'SpO2') && ...
            (~isfield(settings, 'SiteOfCollapse') || (isfield(settings, 'SiteOfCollapse') && settings.SiteOfCollapse == 0))
        disp('No Spo2 Data- Run getEvtsAddOn.m to get Evts Variable');
        SpO2data=NaN;
    end
    
    %% Check FlowSignal - overwrite flow related settings if necessary
    % The following settings are only applicable to Pnasal derived flow
    %  (1) invertflowlist, (2) pnasaldownisinsp, (3) sqrt_scaling, and
    %  (4) scalingexponent.
    % So if the origin of the flow signal is pneumotach data, then the
    % settings do not apply, and we manually overwrite them below.
    % If the origin is pnasal, then the settings do apply.
    
    % if both 'Flow' and 'Pnasal' exist, set pneumotach/pnasal to Flow,
    % in constrast, if only 'Flow' exists then do not change Flow
    if ( sum(strcmp(ChannelsList,'Flow')==1) && ...
            sum(strcmp(ChannelsList,'Pnasal')==1) )
        disp('Detected both Flow and Pnasal signals');
        % if both of these exist, and the user wants to use pnasal,
        % then we need to overwrite the Flow variable, which
        % currently holds pneumotach data, with the pnasal data.
        % if they want pneumotach data, then it's already set as such
        if strcmp(settings.FlowSignal,'SecondPnasal')
            SigT.Flow = SigT.Pnasal;
            disp('Using a secondary Pnasal signal in this analysis');
        else
            disp('Using primary Flow signal rather than secondary Pnasal in this analysis');
        end
    else % both Flow and pnasal don't exist
        % if the user wants to use pnasal, set this as the "flow" signal
        %             if strcmp(settings.FlowSignal, 'Pnasal')
        %                 % we need to set up a "flow'
        %                 DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1)) = ...
        %                 DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pnasal')==1));
        %                 disp('Using Pnasal signal in this analysis');
        %             end
    end
    
    % The flow related settings are manually overwritten for real flow,
    % and are unchanged for pnasal flow.
    if 0&&strcmp(settings.FlowSignal,'Pneumotach')
        settings.invertflowlist = false(size(invertflowlist,1),1);
        settings.Pnasaldownisinsp = 0;
        settings.sqrt_scaling = 0;
        settings.scalingexponent = 1;
    end
    
    %% downsample flow
    if (settings.downsampledFs~=0)&&(settings.downsampledFs~=settings.Fs)
        % if downsampling flow, doing downsample here so that everything
        % uses the same signal. This is then upsampled, as Fs is used in
        % subsequent processing.
        % input: Time, Flow
        % output: Flow, downsampled, but restored to previous Fs
        
        if strcmp(settings.FlowSignal,'SecondPnasal')
            Flow = SigT.Pnasal;
            disp('Using Pnasal signal in this analysis');
        else
            Flow = SigT.Flow;
        end
        
        if 0
            figure(89)
            plot(SigT(:,1),Flow)
            
        end
        
        %FlowOrig = Flow;
        
        if 0 % this is the old resampling/downsampling code
            % but this fails if fs's are not integer multiples
            displaytext = ['WARNING: Resampling flow data'];
            disp(displaytext);
            %                 set(handletext,'String',displaytext); drawnow;
            % downsample
            [p,q] = rat(settings.downsampledFs/originalFs);
            Flow_ds = resample(Flow,p,q);
            % upsample
            %[p,~] = rat(originalFs/settings.downsampledFs);
            Flow_ds2 = resample(Flow_ds,q,p);
            
            Flow=Flow_ds2;
            
        else
            
            displaytext=['Resampling: Flow data from ' num2str(settings.Fs) ' to ' num2str(settings.downsampledFs) 'Hz...(and back)'];
            disp(displaytext);
            %                 set(handletext,'String',displaytext); drawnow;
            TotalTime = (height(SigT)-1)*(1/settings.Fs);
            Time_ = [0:(1/settings.Fs):TotalTime]';
            TimeDS=[0:(1/settings.downsampledFs):TotalTime]';
            FlowDS = interp1(Time_,Flow,TimeDS,'linear'); % downsample
            Flow = interp1(TimeDS,FlowDS,Time_,'linear');  % upsample back to original length.
            clear Time_ TimeDS FlowDS
            %figure();
            %plot(Flow,'r'); hold on;
            %plot(FlowUS, 'b');
        end
        
        % ensure correct length
        mismatch = length(Flow)-height(SigT);
        if mismatch>=0
            % truncate Flow
            Flow = Flow(1:end-mismatch);
        else
            % add spacer to Flow
            Flow = [Flow; zeros(abs(mismatch),1)];
        end
        %then write it back into DataEventHypnog_Mat
        if strcmp(settings.FlowSignal,'SecondPnasal')
            SigT.Pnasal = Flow;
        else
            SigT.Flow = Flow;
        end
        
    end
    
    
    %% Remove artefact in CCH-MC data
    % Some of the CCH-MC studies have a LOT of artefact.
    % Artefact removal per signal is usually performed at ConvertToMat,
    % and the times are determined by the rows in the text _art files,
    % with each row being start and end as seconds since start of study
    % Here, artefact applies only to the nasal pressure/flow signal,
    % but in contrast to doing this at convert to mat, it has no effect
    % on the converted file. It simply replaces the signal with NaN for
    % analysis only.
    % The artefact times in these text files are Time of Day (HH:MM:SS)
    % and require conversion to the time in DataEventHypnogMat
    
    % NOTE: artefact scoring should be moved to summary analysis.
    %
    
    if isfield(settings,'ArtefactPeriodsInCCHMCdata') && settings.ArtefactPeriodsInCCHMCdata==1
        % open text file containing artefact time data.
        % files in this dataset are named LG_0001_EDF2_XHz.mat, and
        % the corresponding artefact file would be named LG_0001_art.txt
        % each period to aretfact is a row, with two columns, as
        %  ArtefactStartTime ArtefactEndTime in HH:MM:SS (24 hour)
        textfilename=[settings.ConvertedDirectory settings.filename(1:end-9) '_art.txt'];
        dt = 1/settings.Fs;
        if exist(textfilename,'file')==2
            displaytext = ['Removing artefact from ' settings.filename(1:end-9)]; disp(displaytext);
            fid = fopen( textfilename );
            cac = textscan( fid, '%s%s', 'CollectOutput', true );
            sts = fclose( fid );
            try
                dataLength = size(SigT,1);
                
                % get the start and end times from DataEventHypnog as actual time
                if istable(SigT)
                    StartTime = datevec(datestr(SigT{1,1}/86400,'HH:MM:SS'));
                    EndTime = datevec(datestr(SigT{end,1}/86400,'HH:MM:SS'));
                else
                    StartTime = datevec(datestr(SigT(1,1)/86400,'HH:MM:SS'));
                    EndTime = datevec(datestr(SigT(end,1)/86400,'HH:MM:SS'));
                end
                
                for AP = 1:size(cac{1,1},1)
                    % read the row data as times for start and end of artefact
                    vec_start = datevec( char(cac{1,1}{AP,1}), 'HH:MM:SS' );
                    vec_end = datevec( char(cac{1,1}{AP,2}), 'HH:MM:SS' );
                    % if it's close to the start, just make it the start
                    if abs(datetime(vec_start)-datetime(StartTime))<minutes(1)
                        vec_start = StartTime;
                    end
                    if vec_start(4)<12
                        vec_start(3) = vec_start(3)+1;
                    end
                    % if it's close to the end, just make it the end
                    if abs(datetime(vec_end)-datetime(EndTime))<minutes(1)
                        vec_end = EndTime;
                    end
                    if vec_end(4)<12
                        vec_end(3) = vec_end(3)+1;
                    end
                    % find the time as seconds from start of study
                    StartSecs = seconds(datetime(vec_start)-datetime(StartTime));
                    EndSecs = seconds(datetime(vec_end)-datetime(StartTime));
                    
                    % convert these times to samples
                    lefti=round(StartSecs/dt+1);
                    righti=round(EndSecs/dt+1);
                    
                    % check these indices are within the data range
                    if lefti>dataLength; disp('Index of artefact start time greater than data length, adjusting...'); lefti=dataLength-1; end
                    if righti>dataLength; disp('Index of artefact end time greater than data length, adjusting...'); righti=dataLength; end
                    
                    % then NaN the flow data in this range
                    if istable(SigT)
                       SigT{lefti:righti,find(strcmp(ChannelsList,'Flow')==1)}=NaN;
                    else
                       SigT(lefti:righti,find(strcmp(ChannelsList,'Flow')==1))=NaN;
                    end

                    disp(['Artefact period from ', num2str(StartSecs),' to ', num2str(EndSecs),'(seconds)'])                        
                end
            catch
                disp('WARNING: CCH-MC artefact not processed');
                textfilename
                AP
            end
        end
    end
    
    %% run windows
    %clear AnalysisIndex LGplusinfo EventsInfo LG_QualityInfo DataOut ArousalDat fitQual MiscEpochData
    [WinT,BreathDataTable,BreathFLDataTable,BreathSnoreTable,LocalSignals,SigT2] ...
        = WindowSelectAndRun(SigT,CPAPoff,CPAP); %pi/16
    %EvtsData added as input for function to calculate position per window
    
    %% Label Site of Collapse
    if isfield(settings,'SiteOfCollapse') && settings.SiteOfCollapse == 1
        BreathDataTable = LabelBreaths(BreathDataTable,Evts);
    end
    
    %% set data to save
    datatosavelist={'WinT','AnalysisIndex','LGplusinfo','EventsInfo','LG_QualityInfo','DataOut','ArousalDat','fitQual','SleepData','CPAPData','StoNData','BreathDataTable','BreathFLDataTable','BreathSnoreTable','Evts','SpikesInfo'};
    if ~isempty(SigT2)
        datatosavelist = [datatosavelist,'SigT2'];
    end
    
    %% save method (append or one file per pt)
    if isfield(settings,'SaveAnalysisOff')&&settings.SaveAnalysisOff==1
        disp('Warning, not saving Analyzed data')
        %do nothing
    else
        %directory switching is needed at present to handle file paths with spaces.
        %localdir = cd; % backup current directory
        try
            displaytext=['Saving analysis data to ' settings.savename '_' num2str(n) '.mat'];
            disp(displaytext);
            savename = [settings.AnalyzedDirectory settings.savename '_' num2str(n) '.mat'];
            
            settings.FileName_In = MATfilename;
            settings.FileName_Out = savename;
            
            %what to save:
            clear datatosave
            datatosave.settingsAnalyzed = settings;
            datatosave.datatosavelist = datatosavelist;
            for i=1:length(datatosavelist)
                if exist(datatosavelist{i})
                    datatosave.(datatosavelist{i}) = eval(datatosavelist{i});
                end
            end
            %Just in case, later not needed, but will need similar code after opening old Analyzed files
            for i=1:length(datatosavelist)
                if exist(datatosavelist{i}) && iscell(eval(datatosavelist{i})) && size(eval(datatosavelist{i}),1)==1 ...
                        && size(eval(datatosavelist{i}),2)==1 && sum(size(eval([datatosavelist{i} '{1}'])))>2 %exist, be a 1x1 cell, contains a cell, that is larger than 1x1.
                    datatosave.(datatosavelist{i}) = eval([datatosavelist{i} '{1}']);
                end
            end
            
            % code to place in Summary Analysis and elsewhere to handle older Analyzed files
            % shifts all fields from from datatosave to current workspace:
            % cellfun(@(x,y) assignin('caller',x,y),fieldnames(datatosave),struct2cell(datatosave))
            
            if ~(isfield(settings,'NSRR') && settings.NSRR==1) %requires moving directly to local directory to save
                localdir = cd; % backup current directory
                try
                    cd(settings.AnalyzedDirectory);
                    save(savename, '-struct', 'datatosave','-v7.3');
                    cd(localdir);
                catch me
                    cd(localdir);
                end
            else
                save(savename, '-struct', 'datatosave','-v7.3')
            end
            displaytext=['Finished Saving Data to ' settings.savename '_' num2str(n) '.mat'];
            disp(displaytext);
            clear datatosave
            
            %if save local signals
            if ~strcmp(settings.GetLocalSignals,'None')
                try
                    %cd(settings.AnalyzedDirectory);
                    savename2 = [settings.AnalyzedDirectory settings.savename '_' num2str(n) '_localsignals.mat'];
                    LocalSignalsData.LocalSignals=LocalSignals;
                    save(savename2, '-struct', 'LocalSignalsData','-v7.3')
                    clear LocalSignals
                catch me
                end
            end
            %cd(localdir); % restore previous current directory
        catch me
            disp('failed saving');
        end
    end
catch me
    displaytext=me.message; disp(me.getReport);
    disp(displaytext);
    %         set(handletext,'String',displaytext); drawnow;
end


%% wrap it up
try close(figure(101)); catch; end
try close(figure(99)); catch; end
try close(figure(2)); catch; end
% NB, these figure closers may also appear elsewhere...

% progressbar(1,1,1); % close progress bar
% delta_t = etime(clock, t_start); % delta in seconds
% D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
% displaytext = ['Complete. Total Analysis time: ', char(D), ' (hh:mm:ss)'];
% disp(displaytext); set(handletext,'String',displaytext); drawnow;
% diary off
end

%% Window Selection
function [WinT,BreathDataTable,BreathFLDataTable,BreathSnoreTable,LocalSignals,SigT2] = ...
    WindowSelectAndRun(SigT,CPAPoff,CPAP)

%% A bunch of variables no longer get assigned in this function and this
% causes and error so I am assigned NaNs to them here
LGplusinfo=nan; EventsInfo=nan; LG_QualityInfo=nan; DataOut=nan; 
ArousalDat=nan; fitQual=nan;
%global settings handletext n winNum ChannelsList ChannelsFs
global settings n winNum ChannelsList ChannelsFs Evts Info
% options
timewav=SigT.Time;
respwav=SigT.Flow;
hypnog=SigT.Epochs;
% if settings.ignoreCPAPdata
%     CPAP=0*CPAP; %a priori knowledge that CPAP is absent.
% end
% filterCPAPdata=1;
% if filterCPAPdata
%     filter_HFcutoff_butter0 = 1/8;
%     filter_order0 = 4;
%     [B_butterHcut,A_butterHcut] = butter(filter_order0,[filter_HFcutoff_butter0]/(1*settings.Fs/2),'low');
%     cpap = filtfilt(B_butterHcut,A_butterHcut,cpap);
%     % time = DataEventHypnog_Mat(:,1);
%     % figure(5); plot(time,[cpap cpapfilt])
% end
duration=length(respwav);
%added a plus one in equation below:
% Dan Vena changed the way number of windows are calculated (original method commented) -
% important for my short DISE studies 07/22/2021
% numwind=floor((duration-settings.windowlength*60*settings.Fs)/(settings.WindowStep*settings.Fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
numwind=ceil((duration-settings.windowlength*60*settings.Fs)/(settings.WindowStep*settings.Fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
if numwind<=0 && duration-settings.windowlength*60*settings.Fs < 0
    numwind=1;
end

disp(['N possible windows: ' num2str(numwind)])
clear allsleep CPAPoff_
    

winnumrange=1:1:numwind;
if isfield(settings,'rerunspecificwindows') && ~isempty(settings.rerunspecificwindows)
    disp('Warning: rerun specific windows on');
    winnumrange=settings.rerunspecificwindows;
end

SleepData=NaN*zeros(numwind,7);
CPAPData=NaN*zeros(numwind,4);
%PositionData = NaN*zeros(numwind,1);

%% Load FlowDrive models
global FlowDriveModels
load('FinalModel_Table.mat');
load('FinalModel_Table10.mat');
load('FLcertModel.mat');
FlowDriveModels.FinalModelTable = FinalModelTable;
FlowDriveModels.FinalModelTable10 = FinalModelTable10;
FlowDriveModels.FLcert = FLcert;


%% spike detector
if settings.spikedetector
    respwav=SigT.Flow;
    if settings.plotfigure
        figure(34); clf(34);
        plot(timewav,respwav,'g');
        hold('on');
    end
    
    [respwav,lrespwavtmp,I]=Spikedetect(timewav,respwav,settings.Fs);
    if 1
        SigT.Flow=respwav;
    else
        'warning not saving spike-removed data'
    end
    ArtifactLengthFull = length(respwav)-lrespwavtmp; % length of nans ie; samples removed
    Artifact.Flow = I;
    %                 upperthresj(j)=upperthres;
    %             lowerthresj(j)=lowerthres;
    disp(['%spikes = ' num2str(100*ArtifactLengthFull/length(respwav))]);
    spikesper = (100*ArtifactLengthFull/length(respwav));
    SpikesInfo=spikesper;
    
    %% Block of old code cut from here, pasted at end of this file
    
else
    SpikesInfo=NaN;
end % spikedetector

%% signal to noise
warning('off'); % turn off, absent signal will cause warning
noisewav = FlowSignalToNoise(timewav,respwav,settings.plotfigure);
warning('on');

if settings.spikedetector
    if settings.plotfigure
        figure(101);
        subplot(4,1,4);
        hold('on');
        %plot(timewav,I*upperthres*0.2);
        hold('off');
        noisewav(I==1)=3; %set spikes as noise
    end
end

%% Run windows
% Disconnect periods (if present and enabled) were set to zeros during ConvertToMat.
% So that these are not incorrectly analyzed as central apneas or other erroneous type,
% here we exclude the effected windows from Analysis.
% First, we make a full length (at Fs) signal of 'disconnect', (0=normal, 1=disconnected)
% then, in the next block, we overwrite criteria to skip were disconnect=1.
if isfield(settings,'AddDisconnectPeriods') && settings.AddDisconnectPeriods==1
    DisconnectData = 0 .* SigT{:,1}; % set as all good to start with
    if isfield(Info, 'DisconnectPeriods')
        try
            cac = Info.DisconnectPeriods;
            for disconnectNum=1:size(cac{1,1},1)
                % get the disconnect time (from start of study)
                [~,~,~,theH,theM,theS]=datevec( cac{1,1}{disconnectNum,1}, 'HH:MM:SS.FFF' );
                % then work out how many samples into datamat the disconnect period starts at
                holdStartsAt = (theH*60*60*settings.Fs) + (theM*60*settings.Fs) + (floor(theS*settings.Fs));
                % get the disconnect duration
                holdTime = cac{1,2}(disconnectNum);  %cac{i,2} is the duration of disconnect
                holdSamples = floor(holdTime*settings.Fs); % how many samples in holdTime
                % inject hold time
                DisconnectData(holdStartsAt:(holdStartsAt+holdSamples))=1;
            end
        catch
            disp('WARNING: failed while setting the disconnect periods');
        end
    end
end



%% Run Window Loop - Set criteria for which windows to run
StoNData.Fnoise=zeros(length(winnumrange),3)*NaN;
for winNum=winnumrange
    
    % note that winNum starts from zero
    %disp(['Processing window: ', num2str(winNum)]);
    % Extract the relevant variables from the particular window.
    lefti=(winNum-1)*settings.WindowStep*settings.Fs+1;
    righti=(winNum-1)*settings.Fs*settings.WindowStep+settings.windowlength*60*settings.Fs;
    
    % Dan Vena added an condition to get end index of last partial window
    % 7/22/2021
    if righti > length(timewav)
        righti = length(timewav);
    end
    
    hypnogwind=hypnog(lefti:righti);
    CPAPwind=CPAP(lefti:righti);
    noisewind=noisewav(lefti:righti);
    %findest longest durations of wakefulness
    clear I I2 I3 I4
    I=hypnogwind==4; I2=[I(1)==1;diff(I)]; I3=find(I2==1); I4=find(I2==-1);
    if length(I4)<length(I3), I4(length(I3),1)=length(hypnogwind)+1; end
    lengthwake=(I4-I3)/settings.Fs;
    if isempty(lengthwake),lengthwake=0; end
    LongestWake=max(lengthwake);
    %Position data in window
    %     PositionData(winNum) = mode(Evtsdata{n}.PositionVariable(lefti:righti));
    
    % If it is N1, N2 or N3 sleep the whole window, calculate loop gain:
    allsleep(winNum)=(sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0);
    CPAPmedian=median(CPAPwind);
    CPAPstd=std(CPAPwind);
    CPAPabs95=max([prctile(CPAPwind,95),-prctile(CPAPwind,5)]);
    if 1
        CPAPoff_(winNum)=CPAPabs95<settings.minabsPmaskforCPAPoff;%CPAPabs95<settings.minabsPmaskforCPAPoff;
    else
        CPAPoff_(winNum)=prctile(CPAPoff(lefti:righti),1); %original method
    end
    %if round(CPAPmedian*100)/100==-1, CPAPoff=0;end
    CPAPData(winNum,:)=[CPAPoff_(winNum) CPAPmedian CPAPstd CPAPabs95];
    Fnoise3=sum(noisewind==3)/length(noisewind);
    Fnoise2(winNum)=sum(noisewind>=2)/length(noisewind);
    Fnoise1=sum(noisewind>=1)/length(noisewind);
    StoNData.Fnoise(winNum,:)=[Fnoise1 Fnoise2(winNum) Fnoise3];
    FNREM1=sum(hypnogwind==2)/length(hypnogwind);
    FNREM2=sum(hypnogwind==1)/length(hypnogwind);
    FNREM3=sum(hypnogwind==0)/length(hypnogwind);
    FREM=sum(hypnogwind==3)/length(hypnogwind);
    FNREM=FNREM1+FNREM2+FNREM3;
    FWAKE=sum(hypnogwind==4)/length(hypnogwind);
    SleepData(winNum,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
    ContainsUnknownSleep(winNum)=sum(hypnogwind==-1|hypnogwind==8|isnan(hypnogwind))>0;
    
    REMdetected(winNum)=(sum(hypnogwind==3)>0);
    
    %choose criteria; mostly just run everything here and tidy up in SummaryAnalysis
    criteria(winNum) = Fnoise2(winNum)<0.1;

    %criteria(winNum) = Fnoise2(winNum)<=1&CPAPoff_(winNum)==1; % for LOFT subject with high spikes

    
    if isfield(settings,'SiteOfCollapse') && settings.SiteOfCollapse == 1
        ContainsUnknownSleep(winNum)=false;
        %        criteria(winNum) =  Fnoise2(winNum)<0.1&CPAPoff_(winNum)==1;% don't care about sleep stages
        criteria(winNum) =  CPAPoff_(winNum)==1;% don't care about sleep stages. Had error with Fnoise so ignoring this too
        
    end
    
    if isfield(settings,'AddDisconnectPeriods') && settings.AddDisconnectPeriods==1
        DisconnectWind=DisconnectData(lefti:righti);
        if sum(DisconnectWind)>0
            criteria(winNum) = 0; % force OFF any window with a disconnect
        end
    end
end



%% Breath level data analysis - LGfromFlowBetaPart1
for winNum=winnumrange
    
    try
        progressbar([],winNum/(length(winnumrange)),[]); % update progress bar
    catch me
        me.message
    end
    
    lefti=(winNum-1)*settings.WindowStep*settings.Fs+1;
    righti=(winNum-1)*settings.Fs*settings.WindowStep+settings.windowlength*60*settings.Fs;
    
    % Dan Vena added an condition to get end index of last partial window
    % 7/22/2021
    if righti > size(SigT,1)
        righti = size(SigT,1);
    end
    
    if ContainsUnknownSleep(winNum)==1
        displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
            ', contains unknown sleep'];
        disp(displaytext);
    end
    if Fnoise2(winNum)>=0.1
        displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
            ', noisy or absent flow signal'];
        disp(displaytext);

    end
    if ~CPAPoff_(winNum)
        displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
            ', CPAP on'];
        disp(displaytext);
    end
    
    AnalysisIndex(winNum,:)=[(winNum-1)*settings.WindowStep*settings.Fs+1 (winNum-1)*settings.Fs*settings.WindowStep+settings.windowlength*60*settings.Fs];
        
    if criteria(winNum)
        displaytext=['Analyzing ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)];
        disp(displaytext);
        try
                [BreathDataTable{winNum},BreathFLDataTable{winNum},BreathSnoreTable{winNum},LocalSignals{winNum},Vflow_out{winNum}]= ...
                    LGfromFlowBetaPart1(SigT(lefti:righti,:)); %,settings,n,winNum,ChannelsList
        catch me
            disp(['error evaluating LGfromFlow or saving its data: ' me.message])
            me.getReport
            if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
                keyboard
            end
            BreathDataTable{winNum} = NaN;
            BreathFLDataTable{winNum} = NaN;
            BreathSnoreTable{winNum} = NaN;
            LocalSignals{winNum} = NaN;
            Vflow_out{winNum} = [];
        end
    else
        BreathDataTable{winNum} = NaN;
        BreathFLDataTable{winNum} = NaN;
        BreathSnoreTable{winNum} = NaN;
        LocalSignals{winNum} = NaN;
    end % if criteria for LGfromFlow analysis
end % window loop
progressbar([],1,[]); % update progress bar

%% Adding Analyzed Window time for each epochs
% chk if completed or not?
if 0
    [Evts]=EpochsTableA(Evts,DataEventHypnog_Mat,ChannelsList,BreathDataTable);
end

%% AutoScore Events
SigT2=[];
if isfield(settings,'AutoScoreRespEvents') &&  settings.AutoScoreRespEvents==1 %separate parts
    disp([num2str(n) ') ' 'warning: autoscoring events']);
    temp{1} = BreathDataTable;
    
    if ~isfield(settings,'minEventDuration')
        settings.minEventDuration=10;
    end    
   
    
    try
        [SigT,ChannelsList,ChannelsFs,Evts] = AutoScoreEvents(temp,SigT,ChannelsList,ChannelsFs,Evts,1,settings.minEventDuration);
    
        plotchanneltext = { ...
            {'Epochs'},{'Flow'},{'EventsResp'},{'EventsRespAuto'},{'VE','VEeupnea'},{'Thorax'},{'Abdomen'},{'SpO2'},{'EventsAr','WakeSleep'}...
            };
        PlotXHzData();
    
        settings.UseAutoScoredEventsForAHI=1;
        RespT = EventRespSigtoRespT(SigT.EventsRespAuto,SigT.Time); %Autoscored Resp Events
        
        BreathDataAvailable = AnalyzedWindowsXHzFromBDT(BreathDataTable,timewav);
        ROImask=(CPAPoff & BreathDataAvailable)*1; % analyzable period per breath analysis and CPAPoff
        
        try
            %running standard version (actual best EEG for arousal analysis based on available and good staging)
            if isfield(Evts,'EvtsAuto') && isfield(Evts.EvtsAuto,'ArT')
                %Evts.EvtsAuto.ArT;
                %do nothing
            else
                Evts.EvtsAuto.ArT = EventRespSigtoRespT(SigT.EventsArWS,SigT.Time);
            end
            Evts.EvtsAuto.RespT = RespT;
            
            %Time = (Evts.EpochT.EpochStart(1):Evts.EpochT.EpochStart(end)+30)';
            %             AnalyzedXHz = AnalyzedWindowsXHzFromBDT(BreathDataTable,timewav);
            %             ROImask=(CPAPoff & AnalyzedXHz)*1; % add analyzable period here 1 and 0's
            %
            
            
            [Evts.EvtsAuto]=getAHIAll(SigT,ROImask,Evts.EvtsAuto,ChannelsList,settings,'EventsRespAuto','EventsArWS');
        end
        
        try
            %running B version (estimated best EEG for arousal analysis; for when staging is not available or good)
            if isfield(Evts,'EvtsAutoB') && isfield(Evts.EvtsAutoB,'ArT')
                %Evts.EvtsAuto.ArT;
                %do nothing
            else
                Evts.EvtsAutoB.ArT = EventRespSigtoRespT(SigT.EventsArWSB,SigT.Time);
            end
            Evts.EvtsAutoB.RespT = RespT;
            [Evts.EvtsAutoB]=getAHIAll(SigT,ROImask,Evts.EvtsAutoB,ChannelsList,settings,'EventsRespAuto','EventsArWSB');
        end
        
        try
            %running version using manually scored arousals paired with auto scored respiratory events:
            if isfield(Evts,'EvtsAutoRespOnly') && isfield(Evts.EvtsAutoRespOnly,'ArT')
                %Evts.EvtsAuto.ArT;
                %do nothing
            else
                Evts.EvtsAutoRespOnly.ArT = EventRespSigtoRespT(SigT.EventsAr,SigT.Time);
            end
            Evts.EvtsAutoRespOnly.RespT = RespT;
            [Evts.EvtsAutoRespOnly]=getAHIAll(SigT,ROImask,Evts.EvtsAutoRespOnly,ChannelsList,settings,'EventsRespAuto','EventsAr');
            
        end
        
        SigT2 = SigT(:,{'Time','VE','VEeupnea','VEpeupnea','EventsRespAuto'});
    catch faultAHI
        disp(faultAHI.getReport);
    end
end

try
    AHI3paA = Evts.EvtsAutoRespOnly.AHIdata2{1}.AllSleepAllPahi(2)
    AHI3paB = Evts.EvtsAutoB.AHIdata2{1}.AllSleepAllPahi(2)
end

%% Loop n windows to get LG - LGfromFlowBetaPart2
if (isfield(settings,'LGfromFlowOriginal') &&  settings.LGfromFlowOriginal==1)
    %do nothing
else
    if settings.flowshapesonly~=1 %1 skips LGmodel section to focus on flowshape
        for winNum=winnumrange
            
            try
                progressbar([],[],winNum/(length(winnumrange))); % update progress bar
            catch me
                me.message
            end
            
            lefti=(winNum-1)*settings.WindowStep*settings.Fs+1;
            righti=(winNum-1)*settings.Fs*settings.WindowStep+settings.windowlength*60*settings.Fs;
            
            % Dan Vena added an condition to get end index of last partial window
            % 7/22/2021
            if righti > size(SigT,1)
                righti = size(SigT,1);
            end
            if ContainsUnknownSleep(winNum)==1
                displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
                    ', contains unknown sleep'];
                disp(displaytext);
            end
            if Fnoise2(winNum)>=0.1
                displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
                    ', noisy or absent flow signal'];
                disp(displaytext);
            end
            if ~CPAPoff_(winNum)
                displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
                    ', CPAP on'];
                disp(displaytext);
            end
            
            if criteria(winNum)
                displaytext=['Analyzing LG ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)];
                disp(displaytext);

                try
                    [WinInfo(winNum,:),BreathDataTable{winNum}]= ...
                        LGfromFlowBetaPart2(SigT(lefti:righti,:),BreathDataTable{winNum},Vflow_out{winNum});
                    
                catch me
                    disp(['error evaluating LGfromFlowPart2 or saving its data: ' me.message])
                    me.getReport
                    if isfield(settings,'HaltOnErrors') && settings.HaltOnErrors==1
                        keyboard
                        % should probably be setting outputs to NaNs, same as else case
                    end
                end
            else
                %disp(['ignore window ' num2str(n) '/' num2str(winNum) ', allsleep=' num2str(allsleep(winNum)) ', CPAPoff=' num2str(CPAPoff(winNum))]);
                WinInfo(winNum,1:22)=NaN;
                BreathDataTable{winNum} = NaN;
                BreathFLDataTable{winNum} = NaN;
                BreathSnoreTable{winNum} = NaN;
            end % if criteria for LGfromFlow analysis
            
        end % window loop
    else %Skipped LG
        WinInfo=NaN;
    end
    progressbar([],[],1); % update progress bar
end

%%
%%
%BreathDataTableLength
for i=1:length(BreathDataTable)
    BreathDataTableLength(i)=size(BreathDataTable{i},2);
end
BreathDataTableLength(BreathDataTableLength==0)=NaN;
BreathDataTableLength(BreathDataTableLength==1)=NaN;
Ntablesdifferentcols = sum(~isnan(BreathDataTableLength) & BreathDataTableLength~=mode(BreathDataTableLength)) %check BreathDataTable is consistent width
if Ntablesdifferentcols>0
    disp('warning: BreathDataTable has different widths for different windows');
end

%% convert everything to Table
WinT=table;

WinT.Time0 = SigT.Time(AnalysisIndex(:,1));

AnalysisIndexT = array2table(AnalysisIndex);
AnalysisIndexT.Properties.VariableNames = {'AnalyzeIndL','AnalyzeIndR'};
WinT = [WinT AnalysisIndexT];

SleepDataT = array2table(SleepData);
SleepDataT.Properties.VariableNames = {'FWake','FNREM','FNREM1','FNREM2','FNREM3','FREM','LongestWake'};
WinT = [WinT SleepDataT];
                
try
WinInfoT = array2table(WinInfo);
WinInfoT.Properties.VariableNames = {'TimeB1','LG0','tau','tau2','LGn','Tn','LG1','LG2','delay','VRA','VRA2','ArThres','MSE','TtotMean','TtotStd','TtotMedian','TtotIQR','Narousals','Nevents','meanNotE','OneMinusRsq','VIsd'};
WinT = [WinT WinInfoT];
catch
end

try
StoNDataT = array2table(StoNData.Fnoise);
StoNDataT.Properties.VariableNames = {'FNoise1','FNoise2','FNoise3'};
WinT = [WinT StoNDataT];
catch
end

try
CPAPDataT = array2table(CPAPData);
CPAPDataT.Properties.VariableNames = {'CPAPoff','CPAPmedian','CPAPstd','CPAPabs95'};
WinT = [WinT CPAPDataT];
catch
end

end % function


