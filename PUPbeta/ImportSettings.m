function settings = ImportSettings(settings,AMasterSpreadsheet)


%% CONVERT AND ANALYSIS SETTINGS
if isfield(settings,'Filenames')
    Filenames=settings.Filenames; % settings.Filenames must exist
end

if ~(isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1) %default, using spreadsheet main page
    [~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
    try
        temp = MasterWorksheet(:,end); %get last col of data [in BF] to indicate EEG referencing
        settings.EverySecondEEGisRef = nancell2mat(temp)==1;   %NaN (empty) taken as zeros, negative
    catch me
        disp('EverySecondisRef column needs to be complete with 1s or 0s, no text')
        %  error();
    end
elseif (isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1)
    settings.EverySecondEEGisRef =0; % there is no default for EEGisRef
end


%% Bypass Xls spreadsheet settings?
TempSpreadsheetSettingsBypass=0;
if isfield(settings,'SpreadsheetSettingsBypass') && settings.SpreadsheetSettingsBypass==1 ||...
        (isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1)
    TempSpreadsheetSettingsBypass=1;
end

%% Look for local settings.csv
try
    if exist([settings.workdir '\PUPstart\settings.csv'],'file') %SS added
        disp('Loading additional settings from settings.csv file (overwrite=off), also bypasses Xls method');
        MSW = table2struct(readtable('settings.csv'));
        settings=mergestructs(settings,MSW,0);
        TempSpreadsheetSettingsBypass=1;
    end
end

%% Import from SpreadSheet
if TempSpreadsheetSettingsBypass==0
settingsSpreadsheet = ImportSettingsSpreadsheet(AMasterSpreadsheet);
settings = mergestructs(settings,settingsSpreadsheet,0); %no overwrite
end

%% POSITION SETTINGS
disp('Using the central position database, found at \PUPbeta\Position\PositionDatabase.xlsx')
PosDatabaseDir = [settings.codedir 'Position' filesep];   % to run in linux backward/forward slash issues
PosDatabaseSpreadsheet = [PosDatabaseDir, 'PositionDatabase.xlsx']; %
[~,~,settings.poscodesdatabase] = xlsread(PosDatabaseSpreadsheet,1,'B2:K55');
settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');

%% CONVERT WORKSHEET
if isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1 %use combination of starthere settings and defaults (below)
    ConvertMatFlag=1; %"List" that describes convert or skip
else %import edf/folder related info from spreadsheet
    [~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
    %consider turning MasterWorksheet into a Table here. Then import a pre-written / saved table (above) for a given cohort.
    
    ChannelNumbers_ = MasterWorksheet(:,17:end-1);
    ChannelNumbers = NaN*ones(size(ChannelNumbers_));
    for i = 1:size(ChannelNumbers_,1)
        for j = 1:size(ChannelNumbers_,2)
            if ischar(ChannelNumbers_{i,j})
                ChannelNumbers(i,j)=NaN;
                %ChannelNumbers_{i,j}=NaN;
            else
                ChannelNumbers(i,j)=ChannelNumbers_{i,j};
            end
        end
    end
    
    ConvertMatFlag = cell2mat(MasterWorksheet(:,9));
    
    %last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
    lastrow = find(1*(~isnan(ConvertMatFlag)),1,'last');
        MasterWorksheet(lastrow+1:end,:)=[];
        ChannelNumbers(lastrow+1:end,:)=[];
        ConvertMatFlag(lastrow+1:end,:)=[];
    settings.MasterWorksheet = MasterWorksheet;
    settings.Filenames = MasterWorksheet(:,2:8); %Source Filenames For Convert + +
    settings.ConvertMatFlag=ConvertMatFlag;
    settings.ChannelNumbers=ChannelNumbers;
    
    
    
    
end


%% ANALYSIS WORKSHEET
if isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1
    %ccould these just be made as defaults? then just do nothing here? 
    settings.analyzelist = 1; % set as always 1 since spreadsheet is not used
    settings.invertflowlist = 0; % set as default 0
    % manually assigning converted directory
    if isfield(settings,'ConvertedDirectory') && isfield(settings,'Filenames') % set in StartHere
        % generating converted filenames
        Filenamestemp=extractBefore(Filenames{:,1},'.edf');
        patients{:,1}=strcat(Filenamestemp,repmat('_XHz.mat',size(Filenamestemp,1),1));
        patients{:,2}=settings.ConvertedDirectory;
        patients{:,3}=settings.protocol; % set in StartHere; %is the SystemPos "column" of data
        settings.patients = patients;
    end
else %Standard Method:
    [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
    NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
    settings.analyzelist = logical(num(:,2));
    settings.invertflowlist = logical(num(:,1));
    settings.protocol = patients(:,3);          %also used in position
    settings.patients = patients;
    [num,~,~] = xlsread(AMasterSpreadsheet,1,'S4:S10003');
    settings.ID = num;
end

%% Default Settings

disp('DefaultSettings: Active');
settingsDefault = ImportSettingsDefault(settings);
settings = mergestructs(settings,settingsDefault,0);

%% Extra Housekeeping
if settings.skiparousalintensity==1
    settings.runarousalintensity=[0 0 0];
end

%detect if Source data directories are full or partial paths
    temp = nan(length(settings.Filenames),1);
    for i=1:size(settings.Filenames,1)
        temp(i)=any(settings.Filenames{i,4}==':');
    end
    SourceDataIsPartialPath=nanmean(temp)<0.5;
    if SourceDataIsPartialPath==1 && ~(isfield(settings,'FileNameNoOverwrite')&& settings.FileNameNoOverwrite==1) %if not a full path (missing ":"), then do this:
        for j=4:6
            for i=1:size(settings.Filenames,1)
                if isnan(settings.Filenames{i,j})
                    settings.Filenames{i,j} = settings.SourceDirectory; %write sourceDirectory to Filenames cols 4 to 6
                else
                    settings.Filenames{i,j} = [settings.SourceDirectory settings.Filenames{i,j}]; %concat sourceDirectory to Filenames cols 4 to 6
                end
            end
        end
    end
    
%% Set channel numbers to dummy variable here to avoid breaking GeneralImport
if isfield(settings,'UseHarmonizedChannelNumbers')&& settings.UseHarmonizedChannelNumbers==1
    settings.ChannelNumbers=0; %faked
end

%%
settings.ImportedSettingsComplete=1;



