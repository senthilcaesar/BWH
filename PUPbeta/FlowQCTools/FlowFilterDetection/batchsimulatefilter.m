
developmentmode=1;
cleardata=1;
lineofID=[1 15];
lineofIDH=[0.001 0.2];

if 0
    load('FilterDetectModels','mdlLPFadj','LPFtransmdls','predrange','mdlHPF','predrangeH');
end

% 'RICCADSA'
% 'Praxis'
% 'G:\Partners Healthcare Dropbox\SATP Group\SaraODB OAT\Converted2019\'
% 'G:\Partners Healthcare Dropbox\SATP Group\AtOStraits\Converted\'


if cleardata==0
    try
        load('FilterDetectTrainingData');
        
        Subject = T.Subject;
        Method = T.Method;
        LPF = T.LPF;
        
        SubjectH = TH.Subject;
        HPF = TH.HPF;
        MethodH = TH.MethodH;
    end
end

%can run script in a loop and it will add new data to the model
%% Select Directory
selpath = uigetdir();
dirx = dir(selpath);
dirx(1:2)=[];
%settings = ImportSettings(settings,AMasterSpreadsheet);

goodfile = ones(length(dirx),1);
for i=1:length(dirx)
    temp = dirx(i).name;
    if ~strcmp(temp(end-6:end),'XHz.mat')
        goodfile(i)=0;
    end
    
end

dirx(goodfile==0)=[];

%% settings
freqinput = [1:15 20 25];
methodinput = [1:4];
methodinputH = [1:2];
ploton=0;

freqinputH = [0 0.005 0.01:0.01:0.2];

%% Clear data
if ~exist('Subject')||cleardata
    Subject=[];
    logPratio=[];
    LPF = [];
    Method = [];
end

if ~exist('SubjectH')||cleardata
    SubjectH=[];
    logPratioH=[];
    HPF=[];
    MethodH = [];
end

%% Run loop and add data
for i=1:length(dirx)
    try
    i
    filedir = [selpath '\' dirx(i).name]
    matObj = matfile(filedir);
    varlist = who(matObj);
    load(filedir);
    
    Flow = DataEventHypnog_Mat(:,2);
    Time = DataEventHypnog_Mat(:,1);
    
    for k=1:length(methodinputH)
        for j=1:length(freqinputH)
            HFmethod=methodinputH(k);
            HFinputtestH=freqinputH(j);
            [temp,tempH]=FlowFilterSimulator(Flow,Time,0,[],HFmethod,HFinputtestH,ploton);
            SubjectH = [SubjectH;i];
            logPratioH = [logPratioH;tempH];
            HPF = [HPF;HFinputtestH];
            MethodH = [MethodH;HFmethod];
        end
        k
    end
    disp('completed Hpass')
    
    for k=1:length(methodinput)
        for j=1:length(freqinput)
            LFmethod=methodinput(k);
            LFinputtest=freqinput(j);
            [temp,tempH]=FlowFilterSimulator(Flow,Time,LFmethod,LFinputtest,0,0,ploton);
            Subject = [Subject;i];
            logPratio = [logPratio;temp];
            LPF = [LPF;LFinputtest];
            Method = [Method;LFmethod];
        end
        k
    end
    disp('completed Lpass')
    end
end

%%
% temp=[repmat(1,17,1);repmat(2,17,1);repmat(3,17,1);repmat(4,17,1)];
% temp2 = repmat(temp,155,1);
% Method = temp2

%%
if 0
    load workspace
end
%% Management
LPF(Method==0)=25;

T = table(LPF,Method,Subject);

Tbackup = T;

Irem = T.LPF>20 | T.LPF>12 & T.LPF<15;
T(Irem,:)=[];
% Irem = T.LPF>12 & T.LPF<15;
% T(Irem,:)=[];

if size(logPratio,1)~=size(T,1)
    if ~exist('logPratiobackup')
        logPratiobackup = logPratio;
    end
    logPratio(Irem,:)=[];
end

predrange=[1 3 6 9 12 15];
Iexclude=any(T.Method==[1,3],2); %train data on clear options: 2 and 4
mdlLPF = fitglm(logPratio(:,[predrange]),log10(T.LPF),'Exclude',Iexclude);
R0 = corrcoef(logPratio(:,[predrange])); VIFs=diag(inv(R0))'; % correlation matrix
T.LPFpred = 10.^predict(mdlLPF,logPratio(:,[predrange]));

figure(33)
for i=1:4
    subplot(2,2,i);
    I = any(T.Method==[i],2);
    scatter(T.LPFpred(I),T.LPF(I)+0.2*randn(length(T.LPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2);
    set(gca,'xscale','log','yscale','log');
    hold on
    plot(lineofID,lineofID,'k');
end

title('Single Model All Data')

%% ManagementH

TH = table(HPF,SubjectH,MethodH);

I = TH.HPF==0;
TH.HPF(I)=0.001;

THbackup = TH;
Irem = TH.HPF>0.1 & TH.HPF<0.15 |  TH.HPF>0.15 & TH.HPF<0.2;
TH(Irem,:)=[];
if size(logPratioH,1)~=size(TH,1)
    if ~exist('logPratiobackupH')
        logPratiobackupH = logPratioH;
    end
    logPratioH(Irem,:)=[];
end



predrangeH=[1:3:10];
IexcludeH=any(TH.MethodH==[1],2); %train data on option 2
mdlHPF = fitglm(logPratioH(:,[predrangeH]),log10((TH.HPF)),'Exclude',IexcludeH);
R0 = corrcoef(logPratioH(:,[predrangeH])); VIFs=diag(inv(R0))'; % correlation matrix
TH.HPFpred = 10.^predict(mdlHPF,logPratioH(:,[predrangeH]));

%     clear Ydata Xdata YdataU95 YdataL95
%     PredCiles = prctile(TH.HPFpred,[0:5:100]);
%     for i=1:length(PredCiles)-1
%         I = TH.HPFpred>=PredCiles(i) & TH.HPFpred<PredCiles(i+1);
%         Ydata(i,1)=prctile(TH.HPF(I),50);
%         Xdata(i,1)=prctile(TH.HPFpred(I),50);
%         YdataU95(i,1)=prctile(TH.HPF(I),97.5);
%         YdataL95(i,1)=prctile(TH.HPF(I),2.5);
%     end

figure(133); clf(133);
for i=1:2
    subplot(1,2,i);
    I = any(TH.MethodH==[i],2);
    scatter(TH.HPFpred(I),TH.HPF(I)+0.00002*randn(length(TH.HPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2);
    set(gca,'xscale','log','yscale','log');
    hold on
    plot(lineofIDH,lineofIDH,'k');
end

title('Single Model All Data')

%%
% save('FilterDetectTrainingData','T','TH','logPratio','logPratioH');

%% Separate predictions from separate models
%predrange=[1 3 6 9 12 15];
%Method 1
figure(34); clf(34);

i=1
I = any(T.Method==[i],2);
R0 = corrcoef(logPratio(:,[predrange])); VIFs=diag(inv(R0))'; % correlation matrix
mdlLPF_LP1 = fitglm(logPratio(I,[predrange]),log10(T.LPF(I)))
[Pred1,temp] = predict(mdlLPF_LP1,logPratio(:,[predrange]));
Pred1=10.^Pred1;
Pred1u95=10.^temp(:,2);

subplot(2,2,i);
scatter(Pred1(I),T.LPF(I)+0.2*randn(length(T.LPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log');
hold on
plot(lineofID,lineofID,'k');
i=2
I = any(T.Method==[i],2);
R0 = corrcoef(logPratio(:,[predrange])); VIFs=diag(inv(R0))'; % correlation matrix
mdlLPF_LP2 = fitglm(logPratio(I,[predrange]),log10(T.LPF(I)));
[Pred2,temp] = predict(mdlLPF_LP2,logPratio(:,[predrange]));
Pred2=10.^Pred2;
Pred2u95=10.^temp(:,2);
subplot(2,2,i);
scatter(Pred2(I),T.LPF(I)+0.2*randn(length(T.LPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log');
hold on
plot(lineofID,lineofID,'k');
i=3
I = any(T.Method==[i],2);
R0 = corrcoef(logPratio(:,[predrange])); VIFs=diag(inv(R0))'; % correlation matrix
mdlLPF_LP2 = fitglm(logPratio(I,[predrange]),log10(T.LPF(I)));
[Pred3,temp] = predict(mdlLPF_LP2,logPratio(:,[predrange]));
Pred3=10.^Pred3;
Pred3u95=10.^temp(:,2);
subplot(2,2,i);
scatter(Pred3(I),T.LPF(I)+0.2*randn(length(T.LPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log');
hold on
plot(lineofID,lineofID,'k');
i=4
I = any(T.Method==[i],2);
R0 = corrcoef(logPratio(:,[predrange])); VIFs=diag(inv(R0))'; % correlation matrix
mdlLPF_LP2 = fitglm(logPratio(I,[predrange]),log10(T.LPF(I)));
[Pred4,temp] = predict(mdlLPF_LP2,logPratio(:,[predrange]));
Pred4=10.^Pred4;
Pred4u95=10.^temp(:,2);
subplot(2,2,i);
scatter(Pred4(I),T.LPF(I)+0.2*randn(length(T.LPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log');
hold on
plot(lineofID,lineofID,'k');
%
T.Pred1=Pred1;
T.Pred2=Pred2;
T.Pred3=Pred3;
T.Pred4=Pred4;

T.Pred1u95=Pred1u95;
T.Pred2u95=Pred2u95;
T.Pred3u95=Pred3u95;
T.Pred4u95=Pred4u95;


%% Translate between methods -- all converted to Full Model
% Convert what 5 Hz 2nd order should be called to allow comparison with 5 Hz first order (~2.7 Hz)
% Y values are model/method 1

xticks = [1 2 3 6 9 15 25];

figure(41); clf(41);
I = any(T.Method==[1],2);
mdl1toX = fitglm(log10(T.Pred1(I)),log10(T.LPFpred(I)));
%mdl2to1 = fitglm((T.Pred2u95(I)),(T.Pred1u95(I)));
subplot(2,2,1);
scatter(T.Pred1(I),T.LPFpred(I),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log','xtick',xticks,'ytick',xticks);
hold on
plot(lineofID,lineofID,'k');

I = any(T.Method==[2],2);
mdl2toX = fitglm(log10(T.Pred2(I)),log10(T.LPFpred(I)));
subplot(2,2,2);
scatter(T.Pred2(I),T.LPFpred(I),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log','xtick',xticks,'ytick',xticks);
hold on
plot(lineofID,lineofID,'k');

I = any(T.Method==[3],2);
mdl3toX = fitglm(log10(T.Pred3(I)),log10(T.LPFpred(I)));
subplot(2,2,3);
scatter(T.Pred3(I),T.LPFpred(I),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log','xtick',xticks,'ytick',xticks);
hold on
plot(lineofID,lineofID,'k');

I = any(T.Method==[4],2);
mdl4toX = fitglm(log10(T.Pred4(I)),log10(T.LPFpred(I)));
subplot(2,2,4);
scatter(T.Pred4(I),T.LPFpred(I),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log','xtick',xticks,'ytick',xticks);
hold on
plot(lineofID,lineofID,'k');

%translations
T.LPFadj=T.LPF;
I = T.Method==1;
T.LPFadj(I)=10.^predict(mdl1toX,log10(T.LPF(I)));
I = T.Method==2;
T.LPFadj(I)=10.^predict(mdl2toX,log10(T.LPF(I)));
I = T.Method==3;
T.LPFadj(I)=10.^predict(mdl3toX,log10(T.LPF(I)));
I = T.Method==4;
T.LPFadj(I)=10.^predict(mdl4toX,log10(T.LPF(I)));

LPFtransmdls{1}=mdl1toX;
LPFtransmdls{2}=mdl2toX;
LPFtransmdls{3}=mdl3toX;
LPFtransmdls{4}=mdl4toX;

%% rerun single model, with translations

%predrange=[1 3 6 9 12 15];
%Iexclude=any(T.Method==[1,3],2);

disp(mdlLPF);
mdlLPFadj = fitglm(logPratio(:,[predrange]),log10(T.LPFadj),'Exclude',Iexclude);
LPFpredadj = 10.^predict(mdlLPFadj,logPratio(:,[predrange]));

figure(36)
colorY = 0.5*rand(1,3);
for i=1:4
    subplot(2,2,i);
    
    I = any(T.Method==[i],2);
    scatter(LPFpredadj(I),T.LPFadj(I)+0.02*randn(length(T.LPFadj(I)),1),8,colorY,'filled','markerfacealpha',0.2)
    set(gca,'xscale','log','yscale','log','xtick',xticks,'ytick',xticks);
    hold on
    plot(lineofID,lineofID,'k');
end

%% HF: Separate predictions from separate models

figure(134);

i=1
I = any(TH.MethodH==[i],2);
R0 = corrcoef(logPratioH(:,[predrangeH])); VIFs=diag(inv(R0))'; % correlation matrix
mdlLPF_HP1 = fitglm(logPratioH(I,[predrangeH]),log10(TH.HPF(I)))
[Pred1,temp] = predict(mdlLPF_HP1,logPratioH(:,[predrangeH]));
Pred1=10.^Pred1;

subplot(1,2,i);
scatter(Pred1(I),TH.HPF(I)+0.002*randn(length(TH.HPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log');
hold on
plot(lineofIDH,lineofIDH,'k');

i=2
I = any(TH.MethodH==[i],2);
R0 = corrcoef(logPratioH(:,[predrangeH])); VIFs=diag(inv(R0))'; % correlation matrix
mdlLPF_HP2 = fitglm(logPratioH(I,[predrangeH]),log10(TH.HPF(I)));
[Pred2,temp] = predict(mdlLPF_HP2,logPratioH(:,[predrangeH]));
Pred2=10.^Pred2;

subplot(1,2,i);
scatter(Pred2(I),TH.HPF(I)+0.002*randn(length(TH.HPF(I)),1),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log');
hold on
plot(lineofIDH,lineofIDH,'k');

% keep this:
TH.Pred1=Pred1;
TH.Pred2=Pred2;

%% Translate between methods -- all converted to Full Model
% Convert what 5 Hz 2nd order should be called to allow comparison with 5 Hz first order (~2.7 Hz)
% Y values are model/method 1

xticksH = [0.001 0.005 0.01 0.02 0.05 0.1 0.2];

figure(141); clf(141);
I = any(TH.MethodH==[1],2);
mdl1toY = fitglm(log10(TH.Pred1(I)),log10(TH.HPFpred(I)));
%mdl2to1 = fitglm((T.Pred2u95(I)),(T.Pred1u95(I)));
subplot(1,2,1);
scatter(TH.Pred1(I),TH.HPFpred(I),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log','xtick',xticksH,'ytick',xticksH);
hold on
plot(lineofIDH,lineofIDH,'k');

I = any(TH.MethodH==[2],2);
mdl2toY = fitglm(log10(TH.Pred2(I)),log10(TH.HPFpred(I)));
subplot(1,2,2);
scatter(TH.Pred2(I),TH.HPFpred(I),8,[0.2 0.5 0.6],'filled','markerfacealpha',0.2)
set(gca,'xscale','log','yscale','log','xtick',xticksH,'ytick',xticksH);
hold on
plot(lineofIDH,lineofIDH,'k');

%translations
TH.HPFadj=TH.HPF;

I = TH.MethodH==1;
TH.HPFadj(I)=10.^predict(mdl1toY,log10(TH.HPF(I)));

I = TH.MethodH==2;
TH.HPFadj(I)=10.^predict(mdl2toY,log10(TH.HPF(I)));

HPFtransmdls{1}=mdl1toY;
HPFtransmdls{2}=mdl2toY;

%% HF no need to rerun fully, but show translations;

%predrange=[1 3 6 9 12 15];
%Iexclude=any(T.Method==[1,3],2);

disp(mdlHPF);
%mdlLPFadj = fitglm(logPratio(:,[predrange]),log10(T.LPFadj),'Exclude',Iexclude);
HPFpredadj = 10.^predict(mdlHPF,logPratioH(:,[predrangeH]));

figure(136); clf(136);
colorY = 0.5*rand(1,3);
for i=1:2
    subplot(1,2,i);
    
    I = any(TH.MethodH==[i],2);
    scatter(HPFpredadj(I),TH.HPFadj(I)+0.00002*randn(length(TH.HPFadj(I)),1),8,colorY,'filled','markerfacealpha',0.2)
    set(gca,'xscale','log','yscale','log','xtick',xticksH,'ytick',xticksH);
    hold on
    plot(lineofIDH,lineofIDH,'k');
end

%% How to use model %%

%save('FilterDetectModels','mdlLPFadj','LPFtransmdls','HPFtransmdls','predrange','mdlHPF','predrangeH');

Flow = DataEventHypnog_Mat(:,2);
Time = DataEventHypnog_Mat(:,1);

%copy this code into somewhere, Flow and Time must exist, output is FlowFilterDetect "structure"

ploton=1;
verbose=1;
if verbose
   disp('    ------Flow quality analysis------    ');
end
FlowFilterDetect = FlowFilterDetector(Flow,Time,ploton,verbose);

% add clip detection:
ClipThresholdFmax=0.90;
ClipFthresholdFnormal=0.002; %higher value removes more (i.e. false) clipping (0.002)
[~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[],ClipThresholdFmax,ClipFthresholdFnormal,1);

FlowFilterDetect.FclippedI=FclippedI;
FlowFilterDetect.FclippedE=FclippedE;
FlowFilterDetect.FclippedTotal=FclippedE+FclippedI;

if verbose
if FlowFilterDetect.FclippedTotal(1)>0.005
    disp('Warning: Flow appears clipped');
else
    disp('Checked: Flow appears free of clipping');
end
disp(['   Clipping fraction: ' num2str(100*FlowFilterDetect.FclippedTotal,2) ' %']);
end












