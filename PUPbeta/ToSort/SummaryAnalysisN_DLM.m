%SummaryAnalysisN
clear all
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');

%%
% try
%     load DB
%     catch me
% end
%%
%Notes

%%
settings.directory='C:\Users\szs88\Dropbox (Partners HealthCare)\MAD-OX\Traits';
settings.savename='MadoxN0';
settings.getCIs=1;
settings.plotfigs=1;
settings.selecttimeofnight=0;
settings.selecttimeofnight_XTiles=2;
settings.selecttimeofnight_NthXTile=2;
settings.selectstate=4; %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
settings.selectposition=[]; %empty for unused, %1=supine in profusion

Mrange=1; % was 9
cutoffvalues=[0.7 0.45 12 40 115 85 95 10];
posclass    =[1 1 1 1 0 0 0 0]; 
units={'','','(s)','(%eupnea)','(%eupnea)','(%eupnea)','(%eupnea)','(%eupnea)'};
longname={'Loop gain, sensitivity','Loop gain, instability','Chemoreflex delay','Vresponse-to-arousal','Arousal threshold',...
    'Collapsibility, Vpassive','Collapsibility, Vactive','Compensation'}

clear Data DataCI DataN AHItotal success Fstates
for i=Mrange%1:M
    try
        close all
        [Data{i},DataCI{i},DataN{i},varlist,AHItotal(i),Fstates{i}] = SummaryAnalysisOnePes(i,settings);
        success(i)=1;
    catch me
        success(i)=0;
    end
    %pause
end
success(1) = 1;

I=find(success==1);
Ifail=find(success==0);

%% Find missing arousal scoring
ArousalDataMissing=zeros(max(Mrange),1);
for i=Mrange
    if success(i)==1
    % DLM if DataN{i}(1,4)==0&DataN{i}(1,5)==0&(DataN{i}(1,1)>3)
    if DataN(1,4)==0&DataN(1,5)==0&(DataN(1,1)>3)
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
        % DLM DataArray(i,:) = Data{i};
        DataArray(i,:) = Data;
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
%%
%% Import quickly in one open/close from Excel

%O2PSGJune2016
col1='A';
col2='AZ'
row1=21; 
row2=2080;
sheet=1;
%read MAT filenames from xls
filexls='J:\PEOPLE\FACULTY\SANDS\MESA\Mesae5_sleeppolysomn_20140617_SS.xlsx';

loadnewvariables={'AHI','BMI','ESS','SBP','HTN','Age','Sex','minsat','ODI4','insomnia','ArI','Race','DBP','MABP'};
colsnewvariables=[ 6     13    11    9     10    22    21    27       28     29         30    31    32 33];

%rangerowsread = [max(rowsnewvariables)];
range=[col1 num2str(row1) ':' col2 num2str(row2)];
[num,txt,raw]=xlsread(filexls,sheet,range);

for i=1:length(loadnewvariables)
    eval(['clear ' loadnewvariables{i} ';'])
    col=colsnewvariables(i);    %Includes flow events
    for j=1:size(raw,1)
        try
            eval([loadnewvariables{i} '(j,1)=raw{j,col};']);
        catch me
            eval([loadnewvariables{i} '(j,1)=NaN;']);
        end
    end
    if 1 %trimming data to length of available PUPdata
        eval([loadnewvariables{i} '((M+1):end)=[];'])
    end
end

%% Analyze 

load('workspace_S1_PX_TX','DataArray','M','FstatesArray');
x=1;
Trait=DataArray(:,x);
Subject=(1:M)';
F1=FstatesArray{x}(:,1);
F3=FstatesArray{x}(:,3);
Frem=FstatesArray{x}(:,4);
Fwake=1-sum(FstatesArray{x},2);
tbl = table(Trait,Subject,F1,F3,Frem,Fwake,'VariableNames',{'Trait','Subject','F1','F3','Frem','Fwake'});
tbl1=tbl;

load('workspace_S2_PX_TX','DataArray','M','FstatesArray');
Trait=DataArray(:,x);
Subject=(1:M)';
F1=FstatesArray{x}(:,1);
F3=FstatesArray{x}(:,3);
Frem=FstatesArray{x}(:,4);
Fwake=1-sum(FstatesArray{x},2);
tbl = table(Trait,Subject,F1,F3,Frem,Fwake,'VariableNames',{'Trait','Subject','F1','F3','Frem','Fwake'});
tbl2=tbl;

load('workspace_S3_PX_TX','DataArray','M','FstatesArray');
Trait=DataArray(:,x);
Subject=(1:M)';
F1=FstatesArray{x}(:,1);
F3=FstatesArray{x}(:,3);
Frem=FstatesArray{x}(:,4);
Fwake=1-sum(FstatesArray{x},2);
tbl = table(Trait,Subject,F1,F3,Frem,Fwake,'VariableNames',{'Trait','Subject','F1','F3','Frem','Fwake'});
tbl3=tbl;

load('workspace_S5_PX_TX','DataArray','M','FstatesArray');
Trait=DataArray(:,x);
Subject=(1:M)';
F1=FstatesArray{x}(:,1);
F3=FstatesArray{x}(:,3);
Frem=FstatesArray{x}(:,4);
Fwake=1-sum(FstatesArray{x},2);
tbl = table(Trait,Subject,F1,F3,Frem,Fwake,'VariableNames',{'Trait','Subject','F1','F3','Frem','Fwake'});
tbl4=tbl;

tbl=[tbl1;tbl2;tbl3;tbl4];
tbl.Subject = nominal(tbl.Subject);

ExcludeVars=[];
lme = fitlme(tbl,'Trait~F1+F3+Frem+Fwake+(1|Subject)','Exclude',ExcludeVars)

[B,Bnames,stats]=randomEffects(lme);
subjectSE=median(double(stats(:,5)))

%% Predict "Outcome"

Outcome=SBP>nanmedian(SBP);
Outcome=ESS;
%vars={'Age','Sex','BMI','Race','Outcome'};
% for i=1:length(vars)
%     eval([vars{i} '_=[' vars{i} ';' vars{i} ';' vars{i} ';' vars{i} '];']);
%     eval(['tbl.' vars{i} '=' vars{i} '_;']); 
% end
x=1
load('workspace_S2_PX_TX','DataArray','M','FstatesArray');
Trait=DataArray(:,x);
Subject=(1:M)';
F1=FstatesArray{x}(:,1);
F3=FstatesArray{x}(:,3);
Frem=FstatesArray{x}(:,4);
Fwake=1-sum(FstatesArray{x},2);
tbl = table(Trait,Age,Sex,BMI,AHI,Outcome,'VariableNames',{'Trait','Age','Sex','BMI','AHI','Outcome'});
%lme = fitlme(tbl,'Outcome~Trait+Age+BMI+Sex+AHI','Exclude',ExcludeVars)
%fitglm(tbl,'Distribution','Binomial')
fitglm(tbl,'Distribution','Normal')


%% Plot versus subject
figure(1); clf(1);
dir = posclass;
    dir(dir==0)=-1;
    clear abnormal
for a=1:length(varlist)
subplot(4,2,a);
for i=1:length(Data)  
    x=(dir(a)*Data{i}(a))>(dir(a)*cutoffvalues(a));
    abnormal(i,a)=x;
    plot(i*[1 1],DataCI{i}(:,a),'-','color',[x 0 0]);
    hold('on');
    
    plot(i,Data{i}(a),'.','markersize',10,'color',[x 0 0]);
    xlim([0 length(Data)+1]); box('off');
    if a==1
        text(i-0.33,1.2,num2str(round(AHItotal(i))),'fontsize',8)
    end
end
%ylabel(varlist{a})
ylabel([longname{a} units{a}])
end
xlabel('subject');
%% multivariate regression
[b,dev,stats] = glmfit(DataArray(:,[1 5 6 8]),AHItotal)
Pvals=stats.p;
%% multivariate regression on CPAP level
TherapeuticPAP = DB(:,3);
Amatrix=DataArray(:,[2 5 6 8]); %Amatrix=DataArray(:,[2 5 6 8]);
if 1
Amatrix=[Amatrix AHItotal']
end
[b,dev,stats] = glmfit(Amatrix,TherapeuticPAP)
Pvals=stats.p

TherapeuticPAP_pred = [Amatrix]*b(2:end)+b(1);
figure(2)
plot(TherapeuticPAP_pred,TherapeuticPAP,'k.','markersize',15); box('off');


%% Abnormality grid 
abnormal(:,3)=[]
%% Abnormality grid 
reorder = [2 7 11 3 12 5 6 8 1 9 10 4]';
C1=1-abnormal(reorder,:);

x=1:size(C1,2);
y=1:size(C1,1);
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
%XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored

figure(9); clf(9);
pcolor(XGrid,YGrid,C2)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
colormap(gcf,'hot');
%cmap = get(gcf,'colormap');
%h=colorbar();
%set(gcf,'colormap',customcmap);
%set(gcf,'colormap','hot');
yticks([]);
xticks([]);
round(AHItotal(reorder)')

