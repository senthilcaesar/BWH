clc; clear; close all;
% below is the directory that contains the converted DISE .mat files
% please replace with the actual directory 
convertedFileDirectory = 'C:\Users\Desktop\PUPbeta\Converted';
%convertedFileDirectory = 'C:\Users\s4235757\Desktop\temp\PAH_PVT\PUPbeta-2021-12-15\PUPbeta\MouthBreathing\Converted';
if exist([convertedFileDirectory filesep 'MouthBreathingMetrics.mat'], 'file')==2
    delete([convertedFileDirectory filesep 'MouthBreathingMetrics.mat']);
end
matFname = dir(fullfile(convertedFileDirectory, '*.mat'));
patientCount = size(matFname,1);

MouthBreathingMetrics = zeros(patientCount,7) * NaN;
for i = 1 : patientCount
    %% Load converted file and get nasal flow
    fileName = matFname(i);
    fileName = [fileName.folder filesep fileName.name];
    text = sprintf('Loading file %s', fileName); disp(text);
    load(fileName);
    flow = SigT.Flow; % Please replace with whatever the nasal flow channel is called, e.g. Pnasal
    flow(isnan(flow)) = [];
    flowPositive = flow;
    flowNegative = -1 * flow;
    flowPositive(flowPositive < 0) = 0;
    flowNegative(flowNegative < 0) = 0;
    %% Inverse 2/3 power transform to simulate cannula
    flowPositive = flowPositive .^ (3/2);
    flowNegative = flowNegative .^ (3/2);
    flowOrig = flow;
    flow = flowPositive - flowNegative;
    %figure; hold on; plot(flow); plot(flowOrig); legend('transformed','original');
    %% Calculate mouth breathing related metrics 
    medianFlow = median(flow);
    meanFlow = mean(flow);
    stdDevFlow = std(flow);
    IQRFlow = iqr(flow);
    percentDifference = abs((medianFlow - meanFlow) / stdDevFlow) * 100;
    if (stdDevFlow == 0)
        percentDifference = 0;
    end
    MouthBreathingPDThreshold = 39.925; % used as threshold to detect mouth breathing
    mouthBreathing = 0; % no mouth breathing by default
    if (percentDifference >= MouthBreathingPDThreshold)
        if (meanFlow > medianFlow)
            % expiratory mouth breathing
            mouthBreathing = 1;
        else
            % inspiratory mouth breathing
            mouthBreathing = -1;
        end
    end
    MouthBreathingMetrics(i,:)=[i, meanFlow, medianFlow, stdDevFlow, IQRFlow, percentDifference, mouthBreathing];
end

try
    MouthBreathingMetrics = array2table(MouthBreathingMetrics);
    MouthBreathingMetrics.Properties.VariableNames = {'patientNum', 'meanFlow', 'medianFlow', 'stdDevFlow', 'IQRFlow', 'percentDifference', 'mouthBreathing'};
    save([convertedFileDirectory filesep 'MouthBreathingMetrics.mat'], 'MouthBreathingMetrics');
catch
end