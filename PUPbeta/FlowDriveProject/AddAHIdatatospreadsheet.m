% script to add AHI data to AnalyzeDataSpreadsheet

close all
clear global AnalyzeDataSpreadsheet
clear
clc

%% Local variables and settings
%AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_Pnasal.xlsx';

datasetsize = 54; % 54 for FlowDrive, 32 for OA

%% pre-processing 
% read spreadsheet (options worksheet)
[~,~,raw] = xlsread(AnalyzeDataSpreadsheet,2,'C3:C32');
settings.savename = char(raw{1});
settings.OutputDataDirectory = char(raw{18});
if settings.OutputDataDirectory(end)~='\'
    settings.OutputDataDirectory=[settings.OutputDataDirectory '\'];
end

filelist = dir([settings.OutputDataDirectory, '*.mat']);

%% get the AHI data from the individual files
AHIdataMat = NaN(datasetsize, 184); % 184 is the curent size of the AHI data 
D=struct('settings',settings);
for pt = 1:datasetsize
    try
        % load an analysis file
        filename = [D.settings.OutputDataDirectory, D.settings.savename, '_', num2str(pt), '.mat'];
        if exist(filename,'file') == 2
            displaytext=['Loading saved data: ', filename]; disp(displaytext);
            load(filename);
                       
            % add this pt AHI data to the growing table
            AHIdataMat(pt,:) = AHIdata{1,1};
                      
            % update settings variable (no need to do this more than once)
            if size(fieldnames(D.settings),1)<3
                D.settings=settings;
                %D.settings.savename = ['PnasalDrive_only25Hz_exp067'];
                if D.settings.OutputDataDirectory(end)~='\'
                    D.settings.OutputDataDirectory=[D.settings.OutputDataDirectory '\'];
                end
            end
            
        else
            displaytext=['No saved data for Pt: ' num2str(pt)]; disp(displaytext);
            % NaN already set for this pt
        end
        % clear vars in readiness for next file
        clearvars -except D pt datasetsize AHIdataMat AnalyzeDataSpreadsheet
    catch me
        disp(me.message);
    end
end

%% write the AHI data to the spreadsheet (no error checking - should close file outside Matlab beforehand)
xlswrite(AnalyzeDataSpreadsheet,AHIdataMat,1,'AA3');
            


