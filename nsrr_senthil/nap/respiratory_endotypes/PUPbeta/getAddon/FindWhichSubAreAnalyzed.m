clear all; close all; clc;

settings.workdir=('E:\MESA\');

% from the folder
mainFolder = uigetdir();
dirx = dir(mainFolder);
subFolder = {dirx.name};
dirx(ismember(subFolder, {'.','..'})) = [];
dirFlags = [dirx.isdir];
dirx(dirFlags)=[];
subNamesT = {dirx.name};
N=3;

% from excelspreadsheet
settings.AMasterdir = [settings.workdir '\PUPStart\'];
AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx']; %


%% for analyzed
if 0
    for ii=1:length(subNamesT)
        temp=extractBetween(subNamesT{ii},'_','.mat');
        %     subNames(ii,1)=str2double(temp(end-N:end));
        temp2=extractAfter(temp,'_');
        subNames(ii,1)=str2double(temp2);
    end
    subNames=sort(subNames);
    [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
    NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
    subNamesX=find(num(:,2)==1);
    
    numtemp=num(:,2);
    
    for ii=1:length(subNamesX)
        if ismember(subNamesX(ii),subNames)
            numtemp(subNamesX(ii))=0;
        end
    end
    
    xlswrite(AMasterSpreadsheet,numtemp,1,'AH4');
    
    %% for converted
else
    for ii=1:length(subNamesT)
        temp2=extractBefore(subNamesT{ii},'_');
        %     subNames(ii,1)=str2double(temp(end-N:end));
        subNames(ii,1)=str2double(temp2);
    end
    subNames=sort(subNames);
    [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AB4:AD10003');
    NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
    subNamestemp=extractBefore(patients(:,2),'_');
    subNamesX=str2double(subNamestemp);
    subNamesX(NaNlist) = 0;
    
    numtemp=ones(length(num),1);
    
    for ii=1:length(subNamesX)
        if ismember(subNamesX(ii),subNames)
            numtemp(ii)=0;
        end
    end
    
    xlswrite(AMasterSpreadsheet,numtemp,1,'AB4');
end



