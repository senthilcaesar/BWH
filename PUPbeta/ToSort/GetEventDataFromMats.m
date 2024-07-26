function GetEventDataFromMats()
%% The function just looks through the event data from individaul patients
% It uses the same analysis worksheet previously used for the production of 
% the individual files. 

%% Inputs / Outputs
% Inputs
%   - file, the individual patient files (analysed .mat file)
%   - file, the analysis spreadsheet used to create these files
%
% Outputs
%   - Text to command window, showing Pt num and Event counts

close all
clear global AnalyzeDataSpreadsheet
clear
clc

%% Local variables and settings
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
datasetsize = 54; % 54 for FlowDrive, 32 for OA

%% pre-processing 
% read spreadsheet (options worksheet)
[~,~,raw] = xlsread(AnalyzeDataSpreadsheet,2,'C3:C32');
settings.savename = char(raw{1});
settings.OutputDataDirectory = char(raw{18});
if settings.OutputDataDirectory(end)~='\'
    settings.OutputDataDirectory=[settings.OutputDataDirectory '\'];
end
if 0 % overwrite settings
    settings.savename = ['FlowDrive_only25Hz'];
    settings.OutputDataDirectory =['C:\PSG_Data\FlowDrive\Analyzed\'];
end
%filelist = dir([settings.OutputDataDirectory, '*.mat']);
D=struct('settings',settings);
%% processing
for pt = 1:datasetsize
    try
        filename = [D.settings.OutputDataDirectory, D.settings.savename, '_', num2str(pt), '.mat'];
        if exist(filename,'file') == 2
            %displaytext=['Loading saved data: ', filename]; disp(displaytext);
            load(filename);
            
            Evt = Evtsdata{1,1}.codes;
            Unknown = find(Evt>10); %|Evt<0);
            ApO = find(Evt==2);
            HypO = find(Evt==4);
            Mixed = find(Evt==5);
            ApC = find(Evt==3);
            HypC = find(Evt==6);
            if ~isempty(Unknown); str=['Patient ', num2str(pt), ' - Unknown: ', num2str(length(ApC))]; disp(str); end
%             if ~isempty(ApC); str=['Patient ', num2str(pt), ' - ApC: ', num2str(length(ApC))]; disp(str); end
%             if ~isempty(HypC); str=['Patient ', num2str(pt), ' - HypC: ', num2str(length(HypC))]; disp(str); end
           
        else
            %displaytext=['No saved data for Pt: ' num2str(pt)]; disp(displaytext);
        end
        % clear vars in readiness for next file
        clearvars -except D pt datasetsize datatosave
    catch me
        disp(me.message);
    end
end

disp('done');
