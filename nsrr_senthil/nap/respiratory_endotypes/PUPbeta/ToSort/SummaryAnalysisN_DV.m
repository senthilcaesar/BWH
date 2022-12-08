%SummaryAnalysisN
% clear all
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');

%%
settings.directory='J:\PEOPLE\FACULTY\SANDS\OralApplianceMM2018\Analyzed';
settings.savename='MMOAT';
settings.getCIs=0;
settings.plotfigs=0;
settings.selecttimeofnight=0;
settings.selecttimeofnight_XTiles=2;
settings.selecttimeofnight_NthXTile=2;
settings.selectstate=4; %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
settings.selectposition=[]; %empty for unused, %1=supine in profusion

Mrange=1:32;
cutoffvalues=[0.7 0.4 12 40 115 85 95 10];
posclass    =[1 1 1 1 0 0 0 0]; 
units={'','','(s)','(%eupnea)','(%eupnea)','(%eupnea)','(%eupnea)','(%eupnea)'};
longname={'Loop gain, sensitivity','Loop gain, instability','Chemoreflex delay','Vresponse-to-arousal','Arousal threshold',...
    'Collapsibility, Vpassive','Collapsibility, Vactive','Compensation'}

clear Data DataCI DataN AHItotal success Fstates Veupnea medianV
for i=Mrange%1:M
    try
        close all
        [Data{i},DataCI{i},DataN{i},varlist,AHItotal(i),Fstates{i},Veupnea(i,1),medianV(i,1)] = SummaryAnalysisOne(i,settings);
        success(i)=1;
    catch me
        success(i)=0;
    end
    %pause
end
I=find(success==1);
Ifail=find(success==0);

%% Find missing arousal scoring
ArousalDataMissing=zeros(max(Mrange),1);
for i=Mrange
    if success(i)==1
    if DataN{i}(1,4)==0&DataN{i}(1,5)==0&(DataN{i}(1,1)>3)
        ArousalDataMissing(i)=1;
    end
    end
end
IArousalDataMissing=find(ArousalDataMissing==1);

%% Overwrite, n.b. changes index
if 0
    Data=Data(I);
    DataCI=DataCI(I);
    AHItotal=AHItotal(I);
    DB=DB(I,:);
end

%% 
clear DataArray 
for i=1:length(Data)
    try
        DataArray(i,:) = Data{i};
    catch me
        DataArray(i,:) = NaN*DataArray(1,:);
    end
    DataArray(i,DataN{i}<3)=NaN;
end
DataArray(IArousalDataMissing,:)=NaN;
%%
clear DataArrayN 
for i=1:length(DataN)
    try
        DataNArray(i,:) = DataN{i};
    catch me
        DataNArray(i,:) = zeros(1,length(DataNArray(1,:)));
    end
    %DataNArray(i,DataN{i}<3)=NaN;
end
DataNArray(IArousalDataMissing,:)=NaN;

%%
clear Upper Lower
for i=1:length(Data)
    try
        Upper(i,:) = DataCI{i}(2,:);
    catch me
        Upper(i,:) = Inf*DataArray(1,:);
    end
end
for i=1:length(Data)
    try
        Lower(i,:) = DataCI{i}(1,:);
    catch me
        Lower(i,:) = -Inf*DataArray(1,:);
    end
end
CI=(Upper-Lower)/2;
%%
clear FstatesArray

for i=1:length(Data)
    for j=1:length(varlist)
        try
            FstatesArray{j}(i,:) = Fstates{i}(:,j)';
        catch me
            FstatesArray{j}(i,:) = [NaN NaN NaN NaN];
        end
    end
end

%ean(CI)
% for i=1:size(DataArray,2)
%     AHIcorreluni(i)=corr(AHItotal(:),DataArray(:,i));
% end
% 1.96*std(DataArray)
LargerData = [DataArray medianV Veupnea];
%%