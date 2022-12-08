% startup
close all
clear
clc

% path to data
OldData = 'C:\PSG_Data\QAO\FS_FD_9.mat'; 
NewData = 'C:\PSG_Data\FlowDrive\Analyzed\FlowDrive_wFlut_w25and125Hz_All.mat';

% selectively load variables of interest
Old = load(OldData,'BreathDataTable');
New = load(NewData,'BreathDataTable');

% match up pts in old and new versions
OldSpreadsheet = 'C:\Users\uqdmann\Dropbox\DOS_Scoring\DOS_Scoring.xlsx';
NewSpreadsheet = 'C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO\AnalyzeDataSpreadsheet.xlsx';
T_old = readtable(OldSpreadsheet,'Range','A1:C31');
T_new = readtable(NewSpreadsheet,'Range','A2:B54');

[f1, f2, f3] = intersect(T_old.Converted, T_new.MATFilename);
[NumOrder] = sortrows([f2 f3]);

DataOutOLDnNEW = cell(1,30);

%% patient processing
for i=1:29

    str = ['Processing old patient ', num2str(NumOrder(i,1)), ' and new patient ', num2str(NumOrder(i,2))]; disp(str);
    str = ['Old patient name: ', char(T_old.Converted(NumOrder(i,1)))]; disp(str);
    str = ['New patient name: ', char(T_new.MATFilename(NumOrder(i,2)))]; disp(str);

    % show the windows that we analysed
    numOLDWindows = length(Old.BreathDataTable{NumOrder(i,1)});
    DataOLDWindows = NaN(numOLDWindows,1);
    for w=1:numOLDWindows
        if (size(Old.BreathDataTable{NumOrder(i,1)}{w},1)==1&&isnan(Old.BreathDataTable{NumOrder(i,1)}{w}))||...
                isempty(Old.BreathDataTable{NumOrder(i,1)}{w})
            continue
        else
            DataOLDWindows(w) = Old.BreathDataTable{NumOrder(i,1)}{w}.Time0(1);
        end
    end
    
    % in the SourceMat .mat files for some studies, there is a variable 
    % "StarttimeSpike", which will need to be subtracted to line up times.
    timeoffset = 0;
    if ~ismember(NumOrder(i,1), [1 2 7 11 14 15 17 18 19 20 21])
        sourcedir = 'C:\PSG_Data\FlowDrive\SourceMat 20171123\';
        sourcename = char(T_old.Converted(NumOrder(i,1)));
        sourcename = sourcename(1:end-4);
        sourcename = [sourcename, '.mat'];
        t = load([sourcedir sourcename],'StarttimeSpike');
        try
            timeoffset = t.StarttimeSpike;
            disp('Adjusting start time in this study (new version)');
        catch me
            disp(me.message);
        end
    end
    
    % show the windows that we analysed
    numNEWWindows = length(New.BreathDataTable{NumOrder(i,2)});
    DataNEWWindows = NaN(numNEWWindows,1);
    for w=1:numNEWWindows
        if (size(New.BreathDataTable{NumOrder(i,2)}{w},1)==1&&isnan(New.BreathDataTable{NumOrder(i,2)}{w}))||...
                isempty(New.BreathDataTable{NumOrder(i,2)}{w})
            continue
        else 
            DataNEWWindows(w) = New.BreathDataTable{NumOrder(i,2)}{w}.Time0(1)-timeoffset;
        end
    end
        
    figure(NumOrder(i,1)); clf(figure(NumOrder(i,1)));
    plot([DataOLDWindows, DataOLDWindows+180], [1 1], 'k', 'LineWidth',10); hold on;
    plot([DataNEWWindows, DataNEWWindows+180], [2 2], 'k', 'LineWidth',10);
    fig = gcf; fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [20.5 2.5 8 2.5];
    box off;
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick=[1 2];
    ax.YLim = [0.5 2.5];
    ax.YTickLabels={'Old' 'New'};    
    ax.FontSize=12;
    tempname = char(T_new.MATFilename(NumOrder(i,2)));
    str = ['Patient: ', tempname(1:end-4), ', Old: ', num2str(NumOrder(i,1)), ', New: ', num2str(NumOrder(i,2))];
    title(str);

    str = ['C:\Users\uqdmann\Dropbox\QAO\Figures\', tempname(1:end-4)];
    %savefig(str);
    saveas(fig, str, 'png');
    
    if size(DataOLDWindows) == size(DataNEWWindows)
        t = [DataOLDWindows,DataNEWWindows];
    else
        mismatch = size(DataOLDWindows,1) - size(DataNEWWindows,1);
        if mismatch<0
            t = [[DataOLDWindows;NaN(abs(mismatch),1)],DataNEWWindows]; 
        else
            t = [DataOLDWindows,[DataNEWWindows;NaN(mismatch,1)]]; 
        end
    end
    DataOutOLDnNEW{NumOrder(i,1)} = t;
end
    
clearvars -except DataOutOLDnNEW

disp('end');    
        
        
                 
        