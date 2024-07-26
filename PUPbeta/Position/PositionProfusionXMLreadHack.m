
%% global variables and settings
global AMasterSpreadsheet settings
%note phasing out F_samp, replace with settings.Fs in subfunctions.

%% master spreadsheet
convertfilename = [settings.workdir, filesep, 'AMasterSpreadsheet.xlsx'];
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = convertfilename;
end
%% read spreadsheet (Files worksheet)
if isfield(settings,'SpreadsheetWorksheetBypass') && settings.SpreadsheetWorksheetBypass==1 %use combination of starthere settings and defaults (below)
    % Do nothing
    MasterWorksheet = settings.MasterWorksheet; % loadfromelsewheresomehow
    
    %no need to reimport for Analysis/Summary etc.
else %import settings from spreadsheet
    
    [~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
    disp('Warning: XLS cols BD, BE now point to LOC, ROC (eyes); BF=Every2ndEEGref; update xls if needed');
    
    %consider turning MasterWorksheet into a Table here. Then import a pre-written / saved table (above) for a given cohort. 
end

ConvertMatFlag = cell2mat(MasterWorksheet(:,9));
%last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
lastrow = find(1*(~isnan(ConvertMatFlag)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];

Filenames = MasterWorksheet(:,2:8);

%%
cleardata=1; %%% set to 1 first time ?
if cleardata
    clear SystemPos Positions
end
for i=1:730 %2907:length(Filenames)
   filedir = [Filenames{i,6} Filenames{i,2}];
       
        S=xml2struct(filedir);
        try
        if size(S.CMPStudyConfig.StepChannels.StepChannel,2)>1
            temp = S.CMPStudyConfig.StepChannels.StepChannel{1}.Labels.Label;
        else
            temp = S.CMPStudyConfig.StepChannels.StepChannel.Labels.Label;
        end
        temp2 = [];
        temp3='';
        for j=1:length(temp)
            temp2{j}=lower(temp{j}.Text);
            temp2{j}(1)=upper(temp2{j}(1));
            if temp2{j}=="Front"
                temp2{j}='Prone';
            end
            if temp2{j}=="Back" || temp2{j}=="Back?"
                temp2{j}='Supine';
            end
            if temp2{j}=="Lefft" || temp2{j}=="Left?"
                temp2{j}='Left';
            end
            if j==1
                temp3 = ['Profusion' temp2{j}];
            else
                temp3 = [temp3 temp2{j}];
            end
        end
        SystemPos{i,1}=temp3;
        Positions{i,1}=temp2;
        disp(i)
        disp(temp3);
        catch
            disp(i)
            disp('Position Data Not Available');
            SystemPos{i,1}=NaN;
        Positions{i,1}=NaN;
            
        end
end
%%
I =  find(~cellfun(@isempty,SystemPos));
try
I2= find(cell2mat(cellfun(@(x)any(isnan(x)),SystemPos,'UniformOutput',false)));
I(I2)=[];
catch
end

unique(SystemPos(I))
if ~isempty(I2)
    if size(unique(SystemPos(I)),1)==1
       SystemPos(I2)= unique(SystemPos(I));
       I=sort([I;I2]);
    else
        disp('unknown position for studies:');
        disp(num2str(I2));
        I(I2)=I2;
        
    end
end
    

%SystemPos(I)

%% Make new database entries
% Positions(I,:)
%Supine	Left	Right	Prone [0 3 1 2]
[SystemPosUnique,Ia] = unique(SystemPos(I));

PositionsUnique = Positions(I(Ia));

Options{1} = {'Back','Supine'};
Options{2} = {'Left'};
Options{3} = {'Right'};
Options{4} = {'Prone','Front'};
Options{5} = {'Unknown'};
Options{6} = {'Upright','Up'};

for i=1:length(SystemPosUnique)
    %inside patient loop
    PositionsA = PositionsUnique{i}
    for k=1:length(Options)
        clear match
        for j=1:length(Options{k})
            match(j,:) = strcmpi(Options{k}(j),PositionsA);
        end
        match=sum(match,1)>0;
        I_ = find(match);
        if isempty(I_)
            Output(k)=NaN;
        else
            Output(k)=I_(1)-1; %map is one less than number found
        end
    end
    UniqueOutput(i,:) = Output;
end

%% to excel sheet system position
sheet='Master';
x1range='AF4';
xlswrite(AMasterSpreadsheet,SystemPos,sheet,x1range)

%% to position database
Tposdatabaseentries = array2table(UniqueOutput);
Tposdatabaseentries.Name = SystemPosUnique;
Tposdatabaseentries = Tposdatabaseentries(:,[7 1:6]);
Tposdatabaseentries.Properties.VariableNames={'SystemPos','Supine','Left','Right','Prone','Unknown','Upright'};

PosDatabaseFile=[settings.codedir 'Position\PositionDatabase.xlsx'];
opts = detectImportOptions(PosDatabaseFile);
opts.DataRange='A3';
PosT=readtable(PosDatabaseFile,opts);

TposdatabaseFinal=Tposdatabaseentries;
clear Flag
for ii=1:length(Tposdatabaseentries.SystemPos)
    if sum(strcmp(Tposdatabaseentries.SystemPos(ii),PosT.Var2))
        TposdatabaseFinal{ii,2:7}=NaN;
        Flag(ii,1)=logical(0);
    else
        Flag(ii,1)=logical(1);
    end
end

% need testing for the following code.
if sum(Flag)>=1
    Tposdatabaseentries2=TposdatabaseFinal(Flag,:);
    
if (PosT.Var2(end) == "")
    disp 'Last row has zero characters: begin to write new position code from here'
    x1range=['B' num2str(length(PosT.Var2)+2)] % 2 header lines
    writetable(Tposdatabaseentries2,PosDatabaseFile,'Range',x1range,'WriteVariableNames',false);
else
    x1range=['B' num2str(length(PosT.Var2)+3)]
    writetable(Tposdatabaseentries2,PosDatabaseFile,'Range',x1range,'WriteVariableNames',false);
end
else
    disp('position code already exists in database. Avoiding overwrite!')
    disp(TposdatabaseFinal{:,1})
end

