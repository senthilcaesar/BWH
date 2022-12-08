function StartHere(resp_output,id,edfname1,xmlname1)

%% debug 
debug=0;

if debug==1 %debug mode compatibile with windows testing
    % eg format--comes upstream from NAP
    resp_output='U:/working/puptest/nap/';
    id= 'XE52YJK42AAXF32-201909-psg-XE52YJK42AAXF32/';
    edfname1='U:/working/puptest/nap/XE52YJK42AAXF32-201909-psg-XE52YJK42AAXF32/data/resp.edf';
    xmlname1='U:/working/puptest/nap/XE52YJK42AAXF32-201909-psg-XE52YJK42AAXF32/annots/harm.annot';
    codedir= 'U:/luna/testing/nsrr/nap/respiratory_endotypes/';
end



%% NAP 
outfold=[resp_output filesep id filesep 'respiratory']
if ~(exist(outfold, 'dir') == 7)
    mkdir(outfold);
end

mydir  = extractBefore(mfilename('fullpath'),'StartHere')
codedir=[mydir filesep]


temp=extractBefore(id,'/');
if ~isempty(temp)
    id=temp;
else
    id =id;
end

[filepath1,edftemp,~] = fileparts(edfname1);
edffold=[filepath1 filesep]
edfname=[edftemp '.edf']

[filepath2,xmltemp,xmlext] = fileparts(xmlname1)
xmlfold=[filepath2 filesep]
xmlname=[xmltemp xmlext]

if debug==1
    xmlfold=[filepath2 '/']
    edffold=[filepath1 '/']
end


sys='NSRR';
syspos='Unknown';

clear global
clc


% This script sets up the necessary directory structure for the analysis
% Initially, the only directory is PUPbeta,which contains the code
% After running this script (by pressing F5), the following directories will be made:
% Analyzed                  - Results of Analysis are output here
% Analyzed\EventAnalyzed    - Results of Event analysis are output here
% Converted                 - Files that have been converted to the new standard format will appear here
% Source                    - This is a space where you can put your original PSG files if you wish
% Summary                   - The output of summary analysis goes here
%
%
% You can use this script after setting up the directories,
% or you can just run PUPbeta.m from the code directory.

global settings PUPbetaGUI AMasterSpreadsheet handletext HarmonizedChannelSpreadsheet

%% Set Working Directory, based on this script running, and add to path
mydir  = mfilename('fullpath');
idcs   = strfind(mydir,filesep);
settings.workdir = mydir(1:idcs(end-1));
addpath(genpath(settings.workdir));

%% Set Code Directory and Version, and add to path
settings.CurrentCodeVersion =  'PUPbeta';
settings.codedir = [codedir filesep settings.CurrentCodeVersion filesep];
if exist(settings.codedir,'dir')~=7 % if we still fail, let the user know
    disp('code directory not found');
    disp(settings.codedir);
    %     keyboard
end
addpath(genpath(settings.codedir));


%% settings for sequential processing

HarmonizedChannelSpreadsheet=[settings.workdir, 'HarmonizedChannelLabels.xlsx'];
settings.HarmonizedChannelSpreadsheetDir= extractBefore(mfilename('fullpath'),'StartHere');

settings.SpreadsheetBypass=1;
settings.RunPUPA2Z=1;
settings.NSRRSave=1;

settings.SpreadsheetSettingsBypass=0; % bypass settings--determine if we need to use or not in NSRR.



%% settings for FileNames/Dir/System/SystemPos

Filenames={edfname,xmlname,xmlname,edffold,xmlfold,xmlfold,sys}
settings.Filenames=Filenames;
settings.FileNameNoOverwrite=1;

settings.Nsrroutfold=outfold;
settings.ConvertedDirectory = [outfold filesep 'Converted' filesep];
settings.AnalyzedDirectory = [outfold filesep 'Analyzed' filesep];
settings.SummaryDirectory = [outfold filesep 'Summary' filesep];
settings.protocol =syspos;
settings.napfolder=[resp_output filesep id filesep]
settings.nsrrid=id;
settings.NoKeepConvertLog=1;
settings.NoKeepAnalysisLog=1;

%% make directories in workdir if they do not exist
if 1
    if ~(exist([outfold filesep 'Analyzed' filesep], 'dir') == 7)
        mkdir([outfold filesep 'Analyzed' filesep]);
    end
    
    if ~(exist([outfold filesep 'Converted' filesep], 'dir') == 7)
        mkdir([outfold filesep 'Converted' filesep]);
    end
    
    if ~(exist([outfold filesep 'Summary' filesep], 'dir') == 7)
        mkdir([outfold filesep 'Summary' filesep]);
    end
end


%% Specialized (optional) settings
%CONVERTED
settings.UseHarmonizedChannelNumbers=1;
settings.BlockEDFload=0;
settings.processEEG = 1; %always
settings.allpowers = 1; %always
settings.EKGplot = 0;
settings.PnasalUprightAuto=1; %Set to 1 if you do not know if inspiration is up or down, the code will guess. Is ~99% accurate. 


%ANALYSIS
settings.useCentralPositionDatabase = 1;
settings.intentionalclipping = 0;
settings.ApplyClippingCorrection = 0; %Analyze--LGfromFlow, use for APPLES (set to 1) not RICCADSA
settings.computeFlowDrive10 = 0; %use for data with signals low pass cut at 10Hz e.g. MESA.


settings.rerunspecificwindows = [];
settings.SaveAnalysisOff=0;
settings.doAHIonly=0; %default 0 
settings.parallelAnalysis=0; %default 0
settings.savelongbreathtables=1; %default 0
settings.plotHBfigs=0;


% SUMMARY
settings.Boot=0;
settings.useCentralPositionDatabase=1;

%% Automated detection settings

settings.WSanalysis = 1; %Wake-sleep analysis
settings.WSArVersion=2;
settings.useWSanalysisToReplaceAr = 1; %0=use original scoring, 1=use best EEG-EvtsAuto, 2=use "predicted best" EEG-EvtsAutoB.
%At startup simply rename new Ar channel as "EventsAr"; rename original/manual
%"EventsAr" to "EventsArManual" during Analyze so that all functions
%use the desired arousal scoring.
settings.AllowOneEEG = 1; % DLM testing for limited EEG channel data

settings.AutoScoreRespEvents=1; % its turned on by default in ImportSettingsDefault; set to 1 here just incase
settings.UseAutoScoredRespEventsForLG = 1;  %turn this on to use autoscored events in endotyping

settings.nanifnoarousals=0; % 1= NaN summary values if arousals are absent


%% launch the front end

if ~isfield(settings,'RunPUPA2Z') && ~settings.RunPUPA2Z
    [PUPbetaGUI] = PUPbeta(); % read from code dir, it's on the path now
    
else % for sequential running- bypass GUI
    set(handletext,'Enable','off')
    settings.Mrange=1; % for nsrr
    settings.SummaryMrange=1;
    
    ConvertN();
    disp('finished conversion')
    
    AnalysisN();
    disp('finished analysis')
    
    % Run summary for NREM--default
    settings.selectstate=4;
    SummaryAnalysisN()
    
    % Run summary for REM
    settings.selectstate=5;
    SummaryAnalysisN()
    
    % Run summary for AllSleep
    settings.selectstate=8;
    SummaryAnalysisN()
    disp('finished summary')
end

exit


