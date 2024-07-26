
%% RUN START HERE FIRST
global settings AMasterSpreadsheet 
 

%% settings
settings.Fs=64;
settings.savefigure=1;

%% get info from spreadsheet
[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
ChannelNumbers_ = MasterWorksheet(:,17:end-1);
ChannelNumbers = NaN*ones(size(ChannelNumbers_));
for i = 1:size(ChannelNumbers_,1)
    for j = 1:size(ChannelNumbers_,2)
        if ischar(ChannelNumbers_{i,j})
            ChannelNumbers(i,j)=NaN;
        else
            ChannelNumbers(i,j)=ChannelNumbers_{i,j};
        end
    end
end
ConvertMatFlag = cell2mat(MasterWorksheet(:,9));
%last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
lastrow = find(1*(~isnan(ConvertMatFlag)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];

ChannelNumbers(lastrow+1:end,:)=[];
ConvertMatFlag(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,2:8);

%% start processing

for n=1:size(Filenames,1)

    %% read edf
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    disp(['Processing:' num2str(n) '=' fname]);
    
    ChannelsList={'Flow'};
    ChannelsFs = NaN*ones(length(ChannelsList),1);
    ChannelFound = zeros(length(ChannelsList),1);
    EDFfilenamedir = [directory fname]
    for i=1:length(ChannelsList)
        displaytext=['Collecting channel:' ChannelsList{i}];
        try
            eval(['[' ChannelsList{i} ',ChannelsFs(i),~,~,LabelTemp,~,Filt,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            displaytext=['Found channel labelled:' LabelTemp ' at ' num2str(ChannelsFs(i)), ' Hz'];
            ChannelFound(i)=1;
        catch me
            displaytext=['No channel:' ChannelsList{i}];
        end
    end
    
    %% Time information
    N_Flow=length(Flow);
    Fs_Flow=ChannelsFs(find(strcmp(ChannelsList,'Flow')==1));
    N_timeXHz = round((N_Flow/Fs_Flow*settings.Fs));
    fid = fopen([directory fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    try
        StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
        if StartTime<43200; StartTime=StartTime+86400; end
    catch me
        StartTime=0;
        disp('Failed to import EDF Start Time, set to 0 sec');
    end
    Time=(StartTime:(1/settings.Fs):StartTime+(N_timeXHz-1)*(1/settings.Fs))'; % This is the time vector associated with the _XHz Flow data.
    EndTime=Time(end);
    ChannelsList = ['Time',ChannelsList];
    ChannelsFs = [settings.Fs;ChannelsFs];
    
    %If needed for handling flow etc:
    TimeFlow=(StartTime:(1/Fs_Flow):StartTime+(N_Flow-1)*(1/Fs_Flow))'; % This is the time vector associated with the _XHz Flow data.

    
    %% flow QC    
     try
         settings.fname=fname;
         
        clear ReportTable NSRRFlowwarnings Ttemp FlowFilterDetect
        
        % 1. find if flow is inverted & overall noise level
        [PrUpright,FnoiseAll]=FlowInvertedDetectTool(Flow,TimeFlow);
        if FnoiseAll>=0.1
            disp('noisy or absent flow signal')
        end
        
        
        % 2. find if a high pass/low pass filter is applied to flow signal
        verbose=0;
        ploton=0;
        FlowRS = interp1(TimeFlow,Flow,Time,'linear');
        FlowFilterDetect = FlowFilterDetector(FlowRS,Time,ploton,verbose);
        clear FlowRS;
       
        LowPassPredict=FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
        HighPassPredict=FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
        
        FlowFilterDetectTemp=struct2table(FlowFilterDetect);
        
        if FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1)<8
            disp('Warning: Flow appears smoothed or downsampled');
        end
        
        if FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1)>0.011
            disp('Warning: Flow appears distorted by baseline adjustment (high pass)');
        end
        
        % 3. find if flow is clipped or not
        [~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[]);
        FclippedTotal=FclippedE+FclippedI;
        
        
        if FclippedTotal>0.005
            disp('Warning: Flow appears clipped');
        end
        disp(['   Clipping fraction: ' num2str(100*FclippedTotal,2) ' %']);
            
        SubId=string(extractBefore(settings.fname,'.edf'));
        
            if ~exist('FlowQcT')
                FlowQcT = table(SubId,Fs_Flow,PrUpright,FnoiseAll,LowPassPredict,HighPassPredict,FclippedTotal);
                FlowFilterDetectT=[table(SubId) FlowFilterDetectTemp];
            else
                temp=table(SubId,Fs_Flow,PrUpright,FnoiseAll,LowPassPredict,HighPassPredict,FclippedTotal);
                FlowQcT = [FlowQcT;temp];
                temp2=[table(SubId) FlowFilterDetectTemp];
                FlowFilterDetectT=[FlowFilterDetectT;temp2];
            end
    catch
    end
    
   

end
saveFiledir=[settings.workdir '\Summary\'];
if ~(exist(saveFiledir, 'dir') == 7)
    mkdir(saveFiledir);
end
saveFilename=[saveFiledir 'FlowQc.mat'];
save(saveFilename,'FlowQcT','FlowFilterDetectT');
disp('finished saving all files')






