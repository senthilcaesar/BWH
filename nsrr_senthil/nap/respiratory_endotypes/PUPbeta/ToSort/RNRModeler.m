%% Model

% Settings
hypop = 1; % use hypopnea data
flowlim = 0; removeApnea = 0; % use flow limitation data / remove apnea
doPlot = 0;

%load subject info
if ~exist('subjectInfo')
    subjectInfo = readtable('AnalyzeDataSpreadsheet.xlsx','Sheet', 1,'Range','B2:CT34');
    subjectInfo.labels = (subjectInfo.TotalSupAHIBaseline - ...
        subjectInfo.TotalSupAHITreatment)./subjectInfo.TotalSupAHIBaseline >= 0.70;
end

%Figure folder
figfolder = 'J:\PEOPLE\POST DOC\VENA\Flow shape\Results and Processed Data\Figures\PSG MM Scope';

%Load data
if hypop == 1 && ~exist('BreathDataTableFinal')
    load('HypopTables.mat') %hypopnea breath data
    
    % Temporary - remove rise time 
    BreathFLDataTableFinal.RiseTime50_O = [];
    BreathFLDataTableFinal.RiseTime50_T = [];
    
elseif flowlim == 1 && ~exist('BreathDataTableFinal')
    load('FLTables.mat') %flow limited breath data
end

% If flow limited, remove apneas

if removeApnea
    apneaIdx = BreathDataTableFinal.Etype == 2;
    BreathFLDataTableFinal(apneaIdx,:) = [];
    InspArrayFinal(apneaIdx,:) = [];
    BreathDataTableFinal(apneaIdx,:) = [];
end

%generate label array and idx of subs with AHI >= 20
labelArray = nan(size(BreathDataTableFinal,1),1);
AHICutIdx = nan(size(BreathDataTableFinal,1),1);
for subnum = 1:size(subjectInfo,1)
    subIdx = BreathDataTableFinal.Subject == str2num(subjectInfo.MATFilename{subnum}(1:end-8));

    labelArray(subIdx) = subjectInfo.labels(subnum);
    AHICutIdx(subIdx) = subjectInfo.TOTAL_SUP_AIH(subnum) >=20;
end

% temporary
ExpArrayFinal = normalizeExp(ExpArrayFinal);

% plot flow shapes
if doPlot
    ylimit = [-2 1.05];
    exploreFlowShapes(BreathDataTableFinal, BreathArrayFinal, labelArray, AHICutIdx, ylimit)
    close all 
end

% Remove breaths based on criteria
% Remove breaths with high AsymIndex and Discontinuite - these have severe 
% palatal prolapse or severe epiglottic collapse which impacts the overall
% flow shape

%Subject level means
%Initialize vars
subcount = 1;
numsubs = sum(subjectInfo.TOTAL_SUP_AIH >=20);
BreathFLDataSub = zeros(numsubs, size(BreathFLDataTableFinal,2));
LabelsSub = zeros(numsubs, 1);
dAHIper = zeros(numsubs, 1);
%Run loop
for subnum = 1:size(subjectInfo,1)
    subIdx = BreathDataTableFinal.Subject == ...
        str2num(subjectInfo.MATFilename{subnum}(1:end-8)) &...
        AHICutIdx == 1;
    
    if sum(subIdx) == 0
        continue
    end
    
    %MIFL data
    BreathFLDataSub(subcount,:) = nanmean(BreathFLDataTableFinal{subIdx,:},1);
    
    %Traits
    
    %Labels
    LabelsSub(subcount,1) = subjectInfo.labels(subnum);
    dAHIper(subcount,1) = subjectInfo.PercentReduction(subnum);
    
    SubsIn(subcount,1) = str2num(subjectInfo.MATFilename{subnum}(1:end-8)); 
    %Count
    subcount = subcount + 1;
end

%Subjects to include
AHIgt20 = subjectInfo.TotalSupAHIBaseline >= 20;

%Traits

SummaryAnalysisN_DV;

pcrit = subjectInfo.PcritBaseline(AHIgt20 == 1);
collapse = LargerData(AHIgt20 == 1,9);
Vpassive = LargerData(AHIgt20 == 1,6);
Vactive = LargerData(AHIgt20 == 1,7);
Compensation = LargerData(AHIgt20 == 1,8);
LG1 = LargerData(AHIgt20 == 1,1);
LGn = LargerData(AHIgt20 == 1,2);
ArTh = LargerData(AHIgt20 == 1,5);
Veupnea = LargerData(AHIgt20 == 1,10);
SubsIn(Veupnea < 5.5)

% Explore features
% exploreFeatures(subjectInfo, BreathDataTableFinal,...
%                     BreathFLDataTableFinal, labelArray, AHICutIdx)

figure
scatter(collapse(LabelsSub == 0), dAHIper(LabelsSub == 0), 'rx')
hold on
scatter(collapse(LabelsSub == 1), dAHIper(LabelsSub == 1), 'bo')
ylabel('Percent change in AHI')
xlabel('Collapsibility (%Veupnea)')

for ii = 1:length(SubsIn)
    text(collapse(ii), dAHIper(ii), num2str(SubsIn(ii)))
end

close all
FeatName = 'AsymIndex_O'
FeatNum = find(strcmp(FeatName, BreathFLDataTableFinal.Properties.VariableNames))

FeatStatsTable.(FeatName) = computeFeatStats(subjectInfo,...
    FeatNum, BreathDataTableFinal, BreathFLDataTableFinal, AHICutIdx);

stats2plot = FeatStatsTable.(FeatName).Properties.VariableNames;
for ii = 2:length(stats2plot)
    plotFeatStats(FeatStatsTable.(FeatName), dAHIper, LabelsSub,...
        stats2plot{ii}, FeatName, SubsIn)
end

plotRelation

% plotFeatureDist(BreathFLDataTableFinal, FeatNum, labelArray, AHICutIdx,...
%     BreathDataTableFinal);


%% Logistic regresion
warning('off', 'all')
% All data
maxdev = chi2inv(.858,1);

opt = statset('display','iter',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');
          
inmodel = sequentialfs(@critfun,BreathFLDataSub,LabelsSub,...
                   'cv','none',...
                   'nullmodel',true,...
                   'nfeatures', 1,...
                   'options',opt,...
                   'direction','forward');
               
model0 = fitglm(BreathFLDataSub(:,inmodel),LabelsSub,'Distribution','binomial');

labelPred = round(model0.Fitted.Probability);
acc = 100*(1-sum(abs(labelPred - LabelsSub))/length(LabelsSub))

for varnum = 1:size(BreathFLDataTableFinal,2)             
    model0 = fitglm(BreathFLDataSub(:,varnum),LabelsSub,'Distribution','binomial');
    results(varnum,:) = model0.Coefficients(2,:); 
end




% Cross-validation with feature selection
clear labelPred predProb
for testSub = 1:numsubs
    % training indices
    trainIdx = [1:numsubs]' ~= testSub;

    % Train with feature selection
    maxdev = chi2inv(.851,1);

    opt = statset('display','off',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');

    inmodel = sequentialfs(@critfun,BreathFLDataSub(trainIdx,:),LabelsSub(trainIdx),...
                           'cv','none',...
                           'nullmodel',true,...
                           'options',opt,...
                           'direction','forward');
    find(inmodel)
    modelTr = fitglm(BreathFLDataSub(trainIdx,inmodel),LabelsSub(trainIdx),'Distribution','binomial');
    testData = BreathFLDataSub(~trainIdx,inmodel);
    predProb(testSub,1) = predict(modelTr, testData);
    labelPred(testSub,1) = round(predProb(testSub,1)); 
end
acc = 100*(1-sum(abs(labelPred - LabelsSub))/length(LabelsSub))


% Correlations with AHI
 mdl = stepwiselm(BreathFLDataSub(:, 95:188),dAHIper,'PEnter',0.1, 'PRemove', 0.25, 'NSteps', 4, 'Verbose', 2)
 [R, P] = corrcoef(BreathFLDataSub(:,48), dAHIper)