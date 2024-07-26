%% FlowDrive analysis was performed in the following stepwise order

% code is currently located at
%mypath = 'C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO';
mypath = 'C:\Users\uqdmann\Dropbox\QAO\FlowAndEdiPesQAO_CodeForPublication';
addpath(mypath);
cd(mypath);
clear

% 1
PUPbeta();

% 1.1
% This was used to convert the files from .mat spike export to the improved
% DataEventHypnog_Mat format
ConvertToMat();

% 1.2
% This was used to analyse each study, for VE, Drive and FlowShape
% assessment on each breath.
Analysis();
Analysis_2(); % this is an alternate method, with slightly differnet pnasal handling

% The folowing subroutines may be required (depending on options set in Analysis)
% 2.1
% Use this function if the output from Analysis was set to save individual files
CombineIndivPtAnalysisMats(); 

% 2.2
% Use this function to produce a simple table that details which channels 
% are available in each study. This is easily copied into the Analysis
% spreadsheet, where the Edi column can be used in the next step.
WhichPatientWhatData();

% 3.1
% This was used to create the FeatureSpace
% note, the production of the FeatureSpace for FlowDrive data has the
% option of using automated reference breaths, or importing Gold-Standard
% manually scored reference breaths.
MakeFeatureSpace();

% 3.2
% Then we need to tidy the FeatureSpace in readiness for the next analysis.
% this involves handling apneaB and low flow breaths, removing NaN's etc.
TidyFeatureSpace();

% 4
% Then we conduct the summary analysis
% this involves another multitude of processing, but in summary we:
%  1. Reduce features by reverse stepwise logistic regression, binary classification
%  2. The N best features are then transformed (sqrt and sq), yeilding Nx3 features
%  3. The highest performing 3N features are selected by linear regression
%  4. The M best features are used to produce a single model

% 4.1 
% FeatureSelection, find the top N (N=50) features,
% using logistic regression of natural features only
% This is primarily intended for the OA data
L1O_LogisticRegression();

% 4.2
% FeatureReduction, identify the top N (50) features,
% using linear regression of natural features, plus transforms
% This is primarily intended for the FD data
L1O_LinearRegression();
L1O_LinearRegression_PlotsforPaper(); % this is the one for the main analysis we have used so far
L1O_LinearRegression_TrainPnasal_TestPnasal_old();
L1O_LinearRegression_TrainFlow_TestPnasal_old();

% these are the most recent versions, this replaces the main analysis above, and the '_old' ones above
L1O_LinearRegression_TrainFlow_TestFlowAndPnasal(); % this one does the L1O run
L1O_LinearRegression_TrainFlow_TestFlowAndPnasal_Analysis(); % this one does the post-run analysis

L1O_LinearRegression_TrainPnasal_TestPnasal(); % this is an updated version of the pnasal part of the main analysis
L1O_LinearRegression_TrainPnasal_TestPnasal_Analysis(); % and the associated post-run analysis

% this script runs the external validation step
ExternalValidation_OA();
