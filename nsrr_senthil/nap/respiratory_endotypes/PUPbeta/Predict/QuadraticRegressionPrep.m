% QuadraticRegressionPrep for "AtoOxyPredict"

%% Preparation
addpath(genpath('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PUPbeta_git\PUPbeta'));
addpath(genpath('C:\Users\sarao\Partners Healthcare Dropbox\SATP Group (satpgroup)\SaraODB_HGNS2\ATS abstract 2019'));
addpath(genpath('G:\Partners Healthcare Dropbox\SATP Group\PUPbeta_git\PUPbeta'));

converttop = @(x) 200*(x/100)./((x/100)+1); % 1 line function that backtransform %change in AHItransfrom to actual %reduction

%%
        
        YcontinuousStr = 'DeltaAHIpT';
        
        T.Responder = T.DeltaAHIp>50;
        N_RandNR = [sum(T.Responder==1) sum(T.Responder==0)]
        %YcontinuousStr = 'Responder';
        
linregnotlogreg=1; %0 = logistic, 1 = linear; note logistic regression yields perfect separation so model yields undefined results.
if linregnotlogreg %not sure if plot works for lin reg fyi
    distr='Normal';
    OutcomeStr = YcontinuousStr;
else
    distr='Binomial';
    OutcomeStr = 'Responder';
end

squaredtermexcl=[];

experiment=1
%Choose which variables you want to include in the model HERE
switch experiment
    case 1
        xvalueslist = {'LGn','VpassiveT','AHIbaseline','Vcomp','ArThresT'};  %% {'AHIbaseline'} %
        Altlist = xvalueslist(3);
        linearonlyterms=[3]; %remove all interactions/squares with this variable e.g. AHI
        %linearonlyterms=[1:5]; %to skip all interactions etc.
        forcespecificterms=[1:5]; %forcing all main variables prevents inflated significance       
    case 2
        xvalueslist = {'VpassiveT','LGn','Vcomp','ArThresT','AHIbaseline','Study'};
        Altlist = xvalueslist(end-1:end);
        linearonlyterms=[5 6]; %remove all interactions/squares with this variable e.g. AHI
        forcespecificterms=[1:6]; %forcing all main variables prevents inflated significance      
    case 3
        xvalueslist = {'VpassiveT','LGn','Vcomp','ArThresT'};
        Altlist = xvalueslist(1);
        linearonlyterms=[]; %remove all interactions/squares with this variable e.g. AHI
        forcespecificterms=[1:4]; %forcing all main variables prevents inflated significance      case 3
    case 4
        xvalueslist = {'dVpassiveT','dLGn','dVcomp','dArThresT','AHIbaseline'};
        Altlist = xvalueslist(1);
        linearonlyterms=[]; %remove all interactions/squares with this variable e.g. AHI
        forcespecificterms=[1:4]; %forcing all main variables prevents inflated significance           
end
warning('off');

Nloops = 50; % to make it nearly a while loop
skipleaveoneout = 0;

modifymaxP = 0;
maxP = 0.157; %0.3
PpenaltyQuadTerms=1; %set to 1 to turn this off 

minterms=1; %
minfwdsstartup=1; %6

minlinearterms=5; %
maxlinearterms=5; %6
minquadraticterms = 2; %cant use this without forcing enough linear terms, needs edits
maxquadraticterms = 2; %6 no good for CV of LW using 6 best, 12 ok

enableFullFwdsFullBwds=0;
    BwdsOnly=1;
    
    StartWithAllLinearTerms=0;
    if StartWithAllLinearTerms
        minfwdsstartup=0;
        enableFullFwdsFullBwds=0;
    end
    
Yvar = T.(YcontinuousStr);

Amatrix = T{:,xvalueslist};

N = height(T);

Exclude = zeros(N,1);
Exclude = isnan(T.Responder);

weights = ones(N,1);
useweights=1; %usually minimal effect of this
if useweights
    weights = zeros(N,1);
    weights(T.Responder==1) = sum(T.Responder==0)./sum(Exclude==0);
    weights(T.Responder==0) = sum(T.Responder==1)./sum(Exclude==0);
end
weights = weights./nanmean(weights(Exclude==0));

protectbaseterms=1;
includequadraticterms=1;
ignoresquaredterms=0;


preliminaryscreenoutP = maxP;
preliminaryscreenoutP = 0;

rangevalidate = [1:N]';
criteriaR = T.Responder*1;
%% More Options
standardizeearlyforplots=1; %1
standardizelateforBetastd=0; %0

standardizefinalmodel=0;

%% 

run QuadraticRegression
run QuadraticRegression2DAllViews