%% calculate AHI for each pt
global settings ChannelsList ColumnHeads
try
    
    AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
    
    %% read spreadsheet
    % files worksheet - used to get converted files
    [~,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'B3:C56');
    % options worksheet - used to get 'settings
    [~,~,raw] = xlsread(AnalyzeDataSpreadsheet,2,'C3:C33');
    settings.plotfigure=logical(raw{6});
    settings.Fs=raw{13};
    settings.ignoreCPAPdata=logical(raw{14});
    if ischar(raw{22})
        eval(['settings.supinepositioncode=[' raw{22} '];']); %e.g. [0 5 -1]
    else
        settings.supinepositioncode=raw{22}; %e.g. 5
    end
    settings.minabsPmaskforCPAPoff=raw{24};
    
    for n=1:size(patients,1)
        %% load file
        directoryn=char(patients(n,2));
        if directoryn(end)~='\'
            directoryn=[directoryn '\'];
        end
        
        MATfilename=[directoryn char(patients(n,1)) '.mat'];
        if exist(MATfilename, 'file') == 2
            load(MATfilename,'DataEventHypnog_Mat','ChannelsList','ColumnHeads');
            
            %% get CPAP data
            [CPAPoff,CPAP]=getCPAP(DataEventHypnog_Mat);
            
            %% get AHI and events data
            [AHIdata{n},Evtsdata{n}]=getAHI(DataEventHypnog_Mat,CPAPoff);
            
            displaytext=['Total AHI: ' num2str(AHIdata{n}(58),4) ' per hr'];
            disp(displaytext); 
            displaytext=['NREM supine AHI: ' num2str(AHIdata{n}(16),4) ' per hr'];
            disp(displaytext); 
   
        else
            continue
        end
    end 
    
catch Afault
    disp(Afault.getReport);
end