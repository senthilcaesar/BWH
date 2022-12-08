%% Batch signal frequency analysis

% open directory (use converted files i.e. _XHz.mat s)
% for every file
% open/load file
% get flow signal
% do frequency analysis
% save plot
% next


%% start up
fh=findall(0,'type','figure');
close(fh)
%clear global settings
clear
clc
warning('off','all'); % the GUI fig can generate warnings, so best off

global settings

%% set path
[currentdir,~,~] = fileparts(mfilename('fullpath'));
cd(currentdir); addpath(genpath(pwd));

%% add paths to filter tool and clipping tool
addpath('M:\Dropbox\PUPbeta_git\PUPbeta\FlowFilterDetection');
addpath('M:\Dropbox\PUPbeta_git\PUPbeta\Clipping');

%% start processing
% open directory (use converted files i.e. _XHz.mat s)
selpath = uigetdir('Select directory with _XHz.mat data');
dirx = dir(fullfile(selpath,'*_XHz.mat')); %dirx(1:2)=[];

LPs = [];
HPs = [];
PrUprights = [];
FnoiseAlls = [];
Clippings = [];
savefigs = 1;
for i=1:length(dirx) % for every file
    
    disp('.'); % add line space in command window for ease of viewing
    filedir = [selpath filesep dirx(i).name];
    str = ['File: ', dirx(i).name, ', loading ...'];  disp(str);
    try
        Fdata = open(filedir); % open
    catch Fopen_fail
        str = ['File: ', dirx(i).name, ', could not open ...'];  disp(str);
        continue
    end
    
    %% get time and flow signals
    if 0 %old
        Time = Fdata.DataEventHypnog_Mat(:,find(strcmp(Fdata.ChannelsList,'Time')==1));
        Flow = Fdata.DataEventHypnog_Mat(:,find(strcmp(Fdata.ChannelsList,'Flow')==1));
    else
        Time = Fdata.SigT.Time;
        Flow = Fdata.SigT.Flow;
    end
    
    %%  Pnasal Polarity Probability Prediction Process
    settings.Fs = 1/(Time(2)-Time(1));
    info = [];
    [PrUpright,FnoiseAll,info]=FlowInvertedDetectTool(Flow,Time,info);
    PrUprights{i,1} = PrUpright;
    FnoiseAlls{i,1} = FnoiseAll;

    if savefigs
        fig = gcf; fig.Units = 'Inches'; fig.Position = [10    6   10    6];
        [filepath,nameOnly,ext] = fileparts(filedir);
        savedir = [selpath, filesep, 'SignalQuality', filesep];
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        saveas(fig, [savedir,nameOnly,'_ExampleFlowTraces'], 'png');
        close(fig);
    end
        
    %% frequency of flow analysis
    ploton=1;   verbose=1; plotdims = [10 6];
    FlowFilterDetect = FlowFilterDetector(Flow,Time,ploton,verbose, plotdims);
    LPs{i,1} = num2str(FlowFilterDetect.LowPassPredict_1stOrderEquivalent(1),'%.2g');
    HPs{i,1} = num2str(FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1),'%.2g');
    if savefigs
        fig = gcf; fig.Units = 'Inches'; fig.Position = [10    6   10    6];
        [filepath,nameOnly,ext] = fileparts(filedir);
        savedir = [selpath, filesep, 'SignalQuality', filesep];
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        saveas(fig, [savedir,nameOnly,'_FlowFrequencyAnalysis'], 'png');
        close(fig);
    end
    
    %% calculation of clipping concentration
    ClipThresholdFmax=0.90;
    ClipFthresholdFnormal=0.002; %higher value removes more (i.e. false) clipping (0.002)
    [~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[],ClipThresholdFmax,ClipFthresholdFnormal,1);
    FclippedTotal=FclippedE+FclippedI;
    Clippings{i,1} = 100*FclippedTotal; % should be <0.5
    if savefigs
        fig = gcf; fig.Units = 'Inches'; fig.Position = [10    6   10    6];
        [filepath,nameOnly,ext] = fileparts(filedir);
        savedir = [selpath, filesep, 'SignalQuality', filesep];
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        saveas(fig, [savedir,nameOnly,'_ClippingAnalysis'], 'png');
        close(fig);
    end

end

keyboard

FN =[];
for i=1:length(dirx) % for every file
    FN{i,1} = dirx(i).name;
end
ConvertedDataSummary = table(FN,LPs,HPs,PrUprights,FnoiseAlls,Clippings, 'VariableNames', {'File','LP','HP','Upright','Noise','Clipping'});

