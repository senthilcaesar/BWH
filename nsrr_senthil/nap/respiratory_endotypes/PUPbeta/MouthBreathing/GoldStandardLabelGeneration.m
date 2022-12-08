clc; clear; close all;
% below is the directory that contains the event analyzed DISE .mat files
% please replace with the actual directory 
eventAnalyzedFileDirectory = 'C:\Users\Desktop\PUPbeta\Analyzed\EventAnalyzed';
%eventAnalyzedFileDirectory = 'C:\Users\s4235757\Desktop\temp\PAH_PVT\PUPbeta-2021-12-15\PUPbeta\MouthBreathing\Analyzed\EventAnalyzed';
if exist([eventAnalyzedFileDirectory filesep 'GoldStandardLabels.mat'], 'file')==2
    delete([eventAnalyzedFileDirectory filesep 'GoldStandardLabels.mat']);
end
% assumes there are two directories called Oral and Nasal within the EventAnalyzed directory
% each directory contains either oral or nasal event analyzed .mat files
eventAnalyzedFileDirectoryOral = [eventAnalyzedFileDirectory filesep 'Oral'];
eventAnalyzedFileDirectoryNasal = [eventAnalyzedFileDirectory filesep 'Nasal'];
matFnameOral = dir(fullfile(eventAnalyzedFileDirectoryOral, '*.mat'));
matFnameNasal = dir(fullfile(eventAnalyzedFileDirectoryNasal, '*.mat'));
patientCount = size(matFnameOral,1);

GoldStandardLabels = zeros(patientCount,2) * NaN;
for i = 1 : patientCount
    %% Load .mat files and get oral and nasal breath counts
    fileNameOral = matFnameOral(i);
    fileNameOral = [fileNameOral.folder filesep fileNameOral.name];
    text = sprintf('Loading file %s', fileNameOral); disp(text);
    load(fileNameOral);
    breathTableOral = BreathDataTableFulls;
    breathCountOral = height(breathTableOral);

    fileNameNasal = matFnameNasal(i);
    fileNameNasal = [fileNameNasal.folder filesep fileNameNasal.name];
    text = sprintf('Loading file %s', fileNameNasal); disp(text);
    load(fileNameNasal);
    breathTableNasal = BreathDataTableFulls;
    breathCountNasal = height(breathTableNasal);

    %% Determine whether mouth breathing has occurred or not
    totalBreaths = breathCountOral + breathCountNasal;
    mouthBreathingPercent = breathCountOral / totalBreaths * 100;
    mouthBreathing = 0; % no mouth breathing by default
    if (mouthBreathingPercent >= 50)
        mouthBreathing = 1;
    end
    GoldStandardLabels(i,:)=[i, mouthBreathing];
end

try
    GoldStandardLabels = array2table(GoldStandardLabels);
    GoldStandardLabels.Properties.VariableNames = {'patientNum', 'mouthBreathing'};
    save([eventAnalyzedFileDirectory filesep 'GoldStandardLabels.mat'], 'GoldStandardLabels');
catch
end