function WhichPatientWhatData()
%% The function compiles a list
% It uses the same analysis worksheet previously used
% It outputs a spreadhseet with filenames and what data each carry

%% Inputs / Outputs
% Inputs
%   - file, the individual patient files (converted datahypnog_mat)
%   - file, the analysis spreadsheet
%
% Outputs
%   - file, one file with list of patients and which channels they have

close all
clear
clc

%% Local variables and settings
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_OA.xlsx';

%% pre-processing 
% read spreadsheet (files worksheet)
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'B3:E150');
D = [];

analyzelist = logical(num(:,2));

%% processing
for pt=1:size(patients,1)
    if analyzelist(pt)==0
        displaytext=['Skipping: n=' num2str(pt) ', ' char(patients(pt,1))];
        disp(displaytext);
        continue
    end
    % load pt
    directoryn=char(patients(pt,2));
    if directoryn(end)~='\'
        directoryn=[directoryn '\'];
    end
    filename=[directoryn, char(patients(pt,1)), '.mat'];
    if exist(filename,'file') == 2
        displaytext=['Working on Pt: ' num2str(pt)]; disp(displaytext);
        load(filename);
        Fch = @(x) nnz(find(strcmp(ChannelsList,x)==1)); % anon fn for finding channel
        % use the channels list to confirm presence of Flow, Pnasal, Edi, Pes, Epi
        D = [D; [pt, Fch('Flow'), Fch('Pnasal'), Fch('Edi'), Fch('Pes'), Fch('Epi')] ];
    else
        displaytext=['No saved data for Pt: ' num2str(pt)]; disp(displaytext);
        % set NaN for this pt
        D = [D; [pt, 0, 0, 0, 0, 0] ];
    end 
end
displaytext=['End Pt processing. Performing post-processing']; disp(displaytext);

%% post-processing
% convert to table
VarNames = {'PT','Flow','Pnasal','Edi','Pes','Epi'};
T = array2table(D, 'VariableNames', VarNames);

% write the results sheet out
writetable(T, 'DataChannels.xlsx');

