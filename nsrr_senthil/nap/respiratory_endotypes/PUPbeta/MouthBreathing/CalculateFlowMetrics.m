% EAS added on 2022-01-04 to detect and handle mouth breathing
% Calculate mouth breathing related metrics for a given flow window
function [meanFlow, medianFlow, stdDevFlow, IQRFlow, meanAD, medianAD, percentDifference, mouthBreathing] = CalculateFlowMetrics(windowFlow)
    global settings
    medianFlow = median(windowFlow);
    meanFlow = mean(windowFlow);
    stdDevFlow = std(windowFlow);
    IQRFlow = iqr(windowFlow);
    meanAD = mad(windowFlow);
    medianAD = mad(windowFlow,1);
    percentDifference = abs((meanFlow - medianFlow) / meanAD) * 100;
    percentDifference(meanAD == 0) = 0;    
    mouthBreathing = 0; % no mouth breathing by default
    if (percentDifference >= settings.MouthBreathingPDThreshold)
        if (meanFlow > medianFlow)
            % expiratory mouth breathing
            mouthBreathing = 1;
        else
            % inspiratory mouth breathing
            mouthBreathing = -1;
        end
    end
end