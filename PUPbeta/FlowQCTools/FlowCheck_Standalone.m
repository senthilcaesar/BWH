close
clear
clc

%% setup
% addpath(genpath('M:\Dropbox\PUPbeta_git\PUPbeta\')); % change or remove as req'd

ploton=1;
verbose=1;
savefigs = 1;

%% file handling
% get file
[file, path] = uigetfile(pwd, '*_XHz.mat'); %
fileNoExt = file(1:length(file)-4);

% load file          
DataIn = load([path, file]);
% find and assign Time and Flow data
Time = DataIn.DataEventHypnog_Mat(:,find(strcmp(DataIn.ChannelsList,'Time')==1));
Flow = DataIn.DataEventHypnog_Mat(:,find(strcmp(DataIn.ChannelsList,'Flow')==1));

%%  Pnasal Polarity Probability Prediction Process
[PrUpright,FnoiseAll]=FlowInvertedDetectTool(Flow,Time);
if verbose
    PrUpright
    FnoiseAll
end

if savefigs
    fig = gcf;
    saveas(fig, [path, fileNoExt, '_FlowExampleTracesPlot'], 'png');
end

%% Flow Filter Frequency Formulation
if verbose
   disp('    ------Flow quality analysis------    ');
end
FlowFilterDetect = FlowFilterDetector(Flow,Time,ploton,verbose);

% add clip detection:
ClipThresholdFmax=0.90;
ClipFthresholdFnormal=0.002; %higher value removes more (i.e. false) clipping (0.002)
[~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[],ClipThresholdFmax,ClipFthresholdFnormal,1);

FlowFilterDetect.FclippedI=FclippedI;
FlowFilterDetect.FclippedE=FclippedE;
FlowFilterDetect.FclippedTotal=FclippedE+FclippedI;

if verbose
if FlowFilterDetect.FclippedTotal(1)>0.005
    disp('Warning: Flow appears clipped');
else
    disp('Checked: Flow appears free of clipping');
end
disp(['   Clipping fraction: ' num2str(100*FlowFilterDetect.FclippedTotal,2) ' %']);
end

if savefigs
    fig = figure(999);
    saveas(fig, [path, fileNoExt, '_FlowClippingPlot'], 'png');
    
    fig = figure(10);
    saveas(fig, [path, fileNoExt, '_FlowFrequencyAnalysisPlot'], 'png');
end





