function CombineIndivPtAnalysisMats()
%% The function combines individaul patient analysis files into one output
% It uses the same analysis worksheet previously used for the production of 
% the individual files. 
% It simply adds the data from each patient to a growing struct.cell array

%% Inputs / Outputs
% Inputs
%   - file, the individual patient files (analysed .mat file)
%   - file, the analysis spreadsheet used to create these files
%
% Outputs
%   - file, one (potentially very large) file with all patient data

close all
clear global AnalyzeDataSpreadsheet
clear
clc

%% Local variables and settings
if verLessThan('matlab', '9.2')   % 2017a is ver 9.2,
    % i.e. matlab 2016b or anything older
    if 1
        AnalyzeDataSpreadsheet  = 'AnalyzeDataSpreadsheet_MESA.xlsx';
        %AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_OA.xlsx';
        datasetsize = 100; % 54 for FlowDrive, 32 for OA, 100 for MESA
    else
        AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_Pnasal.xlsx';
        datasetsize = 54; % 54 for FlowDrive, 32 for OA
    end
else
    % i.e. matlab 2017a (only)
    AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_pneumotach.xlsx';
    datasetsize = 54; % 54 for FlowDrive, 32 for OA
end


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
  
filelist = dir([settings.OutputDataDirectory, '*.mat']);

%% processing
datatosave={'AHIdata','SpO2data','AnalysisIndex','LGplusinfo','EventsInfo',...
    'LG_QualityInfo','DataOut','ArousalDat','fitQual','MiscEpochData','SleepData',...
    'CPAPData','StoNData','BreathDataTable','BreathFLDataTable','LocalSignals','Evtsdata'};
D=struct('settings',settings);
for pt = 1:datasetsize
    try
        filename = [D.settings.OutputDataDirectory, D.settings.savename, '_', num2str(pt), '.mat'];
        if exist(filename,'file') == 2
            displaytext=['Loading saved data: ', filename]; disp(displaytext);
            load(filename);
            
            %str=['FileIn: ', settings.FileName_In]; disp(str);
            %str=['FileOut: ', settings.FileName_Out]; disp(str);
            
            % add this pt to the megafile 'D'
            for i=1:length(datatosave)
                eval(['D.',datatosave{i},'{1,',num2str(pt),'}=',datatosave{i},';']);  
            end
            
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
            % set NaN for this pt in the cell aray
            for i=1:length(datatosave)
                eval(['D.',datatosave{i},'{1,',num2str(pt),'}=NaN;']); 
            end
        end
        % clear vars in readiness for next file
        clearvars -except D pt datasetsize datatosave
    catch me
        disp(me.message);
    end
end

%% post-processing
% push each variable back to the workspace (i.e destruct the D struct) :)
settings = D.settings;
pts = length(D.AHIdata);
for i=1:length(datatosave)
    for pt=1:pts
        var_as_text = ['D.', datatosave{i},'{1,',num2str(pt),'}'];
        if eval(['iscell(',var_as_text,')'])
            eval([datatosave{i}, '{1,',num2str(pt),'} = D.',datatosave{i},'{1,',num2str(pt),'}{1,1};']);
        else
            eval([datatosave{i}, '{1,',num2str(pt),'} = NaN;']);
        end
    end
end

% tidy up
clear D pt i datasetsize datatosave

%% save output
disp(' ');
displaytext=['Saving combined data']; disp(displaytext);
if 0
save([settings.OutputDataDirectory, settings.savename, '_All'],'-v7');
else
    save([settings.OutputDataDirectory, settings.savename,'_v73'],'-v7.3');
end
% DLM comment
% -v7 wins !
% same data saved with options 'v7' and 'v7.3'. 
% Filesize for each was 300Mb and 1.3Gb respectively. 
% BUT, v7 doesn't support fgiles >2Gb in size

end

%% Description of analysis output variables
% each variable is a 1 x p cell for each of the p patients
% in any given patient, there may be w windows analyzed
% each window has b breaths, and e events, and t time
% 'SleepData' {1,p}(w,7)
% 'LG_QualityInfo' {1,p}(w,13)
% 'DataOut' {1,p}{1,w}(b,24) % not used anymore
% 'BreathDataTable' {1,p}{1,w}(b,22)
% 'BreathFLDataTable' {1,p}{1,w}(b,100)
% 'LocalSignals' {1,p}{1,w}(t,2) col 1 is Flow, col 2 is Edi
% Note that the BreathDataTable and BreathFLDataTable may be
% shorter than the others if no data in the end windows

% Variables available, but currently not loaded
% 'MiscEpochData' {1,p}{1,w}, each being a single value
% 'AnalysisIndex' {1,p} each being (w,2) pair of start stop values
% 'LGplusinfo' {1,p}(w,17)
% 'EventsInfo' {1,p}(w,10)
% 'ArousalDat' {1,p}{1,w}(e,2) where e is the number of arousals?
% 'AHIdata' {1,p}(1,128)
% 'CPAPData' {1,p}(w,4)