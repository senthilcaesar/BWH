function settingsDefault = ImportSettingsDefault(settings)

%% ChannelCheck
settingsDefault.LabelsOnly=0; % default=0, 0 will write Label|Transducer|Fs to excel, 1 will just write the Label to excel

%% Directories
settingsDefault.SourceDirectory = [settings.workdir 'Source\'];
settingsDefault.ConvertedDirectory = [settings.workdir 'Converted\'];
settingsDefault.AnalyzedDirectory = [settings.workdir 'Analyzed\'];
settingsDefault.SummaryDirectory = [settings.workdir 'Summary\'];

%% GENERAL
settingsDefault.HaltOnErrors = 0;       % default = 0, 0 or no setting does nothing, 1 invokes keyboard at catch statements
settingsDefault.verbose=logical(1);     % useful during trouble shooting

%% CONVERT
settingsDefault.UseHarmonizedChannelNumbers=0;
settingsDefault.Fs=125;
settingsDefault.reportEDFprefilters = 0; 	% default = 0, 0 or no setting does nothing, 1 reports edf channel prefilters to cmd window
settingsDefault.saveFlowQc=1;
settingsDefault.RunInReverse = 0; % always
settingsDefault.EverySecondEKGisRef = logical(0);
settingsDefault.FsLeak = 0;
settingsDefault.skiparousalintensity = 0;
settingsDefault.AnalyzeSnore = 0;
settingsDefault.scoredarousalsinwake = 0;
settingsDefault.ProfusionPositionUpdate = 0; % OFF by default. Use requires summary export from Profusion named as Filename_posP.txt
settingsDefault.processEEG = 1; % always
settingsDefault.allpowers = 1; % always
settingsDefault.EKGplot = 0; % not important
settingsDefault.WSanalysis = 1; %converted
settingsDefault.WSArVersion = 2; %Always use 2; 1 was arousal detection based on PrW changes; 2 is newer detection that gives less weight to falling theta/delta, more to rising alpha/beta 
settingsDefault.findStableBreathing = 1; % this finds stable breathing periods and adds to DataHypnogMat
settingsDefault.runarousalintensity = [1 0 0]; %leave secondary arousal intensity analysis off
settingsDefault.phasicREMfromEOG=1;
settingsDefault.PositionCorrectionTimeIsSinceStartRecording=1;
settingsDefault.plotHBfigs=0;
settingsDefault.parallelConvert = 0;

%% ANALYSIS
settingsDefault.savename = 'PUP';
% settingsDefault.AnalyzedDirectory = char(raw{18});
settingsDefault.Pnasaldownisinsp=1; %Default but could be wrong!
settingsDefault.LGfromFlowVersion = 'LGfromFlowBeta';
settingsDefault.sqrt_scaling=1; 
settingsDefault.scalingexponent=0.67;
settingsDefault.saveplots=logical(0);
settingsDefault.plotfigure=logical(0);
if settingsDefault.saveplots % force plotfigure on if saveplots is on
    settingsDefault.plotfigure = true;
end
settingsDefault.plotfiguresqrtscaling=0;
settingsDefault.plotfigureLGparameters=0;
settingsDefault.savefigas='saveasPNG';  % options are 'saveasTIFF' | 'saveasFIG' | 'saveasPNG'
settingsDefault.usescoredcentralapneas=logical(1);
settingsDefault.eventsarebreathsfullywithinmargins=logical(1);
settingsDefault.havescoredcentralhypops=logical(0);
settingsDefault.manualscoringtouchups=logical(0);
settingsDefault.maxdelaybreaths=5;
settingsDefault.windowlength=7;
settingsDefault.WindowStep=120;
settingsDefault.ignoreCPAPdata=logical(1);
settingsDefault.exportresultstoxls=logical(0);
settingsDefault.handlemixedeventsseparately=logical(1);
settingsDefault.longestwakeduration=30;
settingsDefault.ARmodel=logical(1);
settingsDefault.AnalyzeNREMonly=logical(0);
settingsDefault.findcentralhypopneasandapneas=logical(0);
settingsDefault.PreLowPass=0;   % default 0, if value, then applies LP filter (in Hz) prior to analysis processing
settingsDefault.PreHighPass=0;  % generally, this should never be switched on, this is for testing only
settingsDefault.minabsPmaskforCPAPoff=1.2;
settingsDefault.flowshapesonly=0; %SS: not boolean
settingsDefault.GetLocalSignals='None';
settingsDefault.SavePerSubject=1;
settingsDefault.FlowSignal = 'Normal';
settingsDefault.downsampledFs = 25;
settingsDefault.spikedetector = logical(0);
settingsDefault.plotbreathdetectionfigures = logical(0); 
settingsDefault.modBB_i_start = logical(1);
settingsDefault.nearerzerocrossings = 0;
settingsDefault.AnalyzeSnore = 0;
settingsDefault.AutoScoreRespEvents=1;              % turned on by default 8/20/21
settingsDefault.minEventDuration = 10;
settingsDefault.ApplyClippingCorrection = 1; 		% turned on by default 8/20/21
settingsDefault.MouthBreathingDetect = 1; 			% EAS added on 2022-01-04 to detect mouth breathing
settingsDefault.MouthBreathingPDThreshold = 55.719; % EAS added on 2022-01-04 to detect mouth breathing
%settings.PnasalUprightAuto==0 						%RMA: should we use this yet?
settingsDefault.UseAutoScoredRespEventsForLG=0;     %turn this on to score respiratory events
%settings.useWSanalysisToReplaceAr==1 				%uses autoscored arousals, standard (EEG channel selection is based on agreement [wake v sleep] with visual staging)
%settings.useWSanalysisToReplaceAr==2 				%uses autoscored arousals, uses autoselected "best' EEG (use for when there is no staging, or visual staging is worse than no info) 
settingsDefault.DriftEstimation=1;
settingsDefault.parallelAnalysis = 0;


%% Other settings that are intentionally not defined by default, but listed here
%settingsDefault.lowVEareObstructiveEvents = 0.5; %will force small breaths to be scored as "obstructed" (unknown drive) for phenotyping
    %use if using manual scoring that might be unreliable


%% SUMMARY ANALYSIS

%settingsDefault.Mrange=1; %single patient; 
%removed by SS 7/12/21 because this clashes %with spreadhseet read

%settingsDefault.SummaryDirectory=char(raw{11}); %SummaryAnalysisDirectory
%settings.AnalyzedDirectory=char(raw{1});
%settings.directory=settings.AnalyzedDirectory;

% summaryanalysisMulti % outputs a ridiculous table of combos requested
% supress outputs until end, need to keep N for each combo
% Pos       State
% all       all
% all       NREM
% supine    all
% supine    NREM
% nonsupine all
% currently no plans to do N1, N2, N3 as unique analyses
% suggest this could be done on as-needs basis if ever req'd


settingsDefault.getCIs=0;
settingsDefault.plotfigs=1;
settingsDefault.selecttimeofnight=0;
settingsDefault.selecttimeofnight_XTiles=double(2);
settingsDefault.selecttimeofnight_NthXTile=double(2);

settingsDefault.selectstate=double(4);      %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
settingsDefault.selectposition = 'All';
% consider renaming selectstateSummaryAnalysis, so it is clearly distinct from selectstateEventAnalysis
% consider renaming selectposSummaryAnalysis, so it is clearly distinct from selectposEventAnalysis

settingsDefault.nanifnoarousals=double(1);
settingsDefault.Vpassive1d=double(1);
settingsDefault.method = 1; %1=standard analysis (estimated drive), 2=analysis with Edi drive
settingsDefault.comb=[0 1 30]; %disabled by default: [enable, odd1even2, combwidth(min)]

settingsDefault.minNeventsLG = 1; 
settingsDefault.minNeventsUA = 1; 
settingsDefault.maxwakethresLG = 30; 
settingsDefault.maxwakethresUA = 300;       % more tolerant to more wake for UA measures

settingsDefault.Nciles = 10;                %default value, alternative--for smoother data--is 131

settingsDefault.runcompare = 0;
settingsDefault.Comparename = 'CompareXY';
settingsDefault.constantEupnea=0; %1
settingsDefault.normalizeusingconstantEupnea=0; %1
settingsDefault.runcompare2 = 0;

settingsDefault.minNeventsLG = 1;
settingsDefault.minNeventsUA = 1;           %used to be zero, not clear when it was changed to 1 (before July 2019).
settingsDefault.maxwakeLG = 30;
settingsDefault.maxwakeUA = 300;

%settings.scoredarousalsinwake = 1;
settingsDefault.useBDT = 1;
settingsDefault.breathlevelinclusion = 0;
settingsDefault.DriveincmH2OnotEupnea = 0;

settingsDefault.drivecol=19;
settingsDefault.minwindows=3;
settingsDefault.EdimtaCalibration=0;        %0 default; set to 1 for unbiased comparison between states
settingsDefault.EdiAutoLinearization=1;


%% EVENT ANALYSIS EXTRA SETTINGS
settingsDefault.selecteventtype=1;          %1=all 2=hyp only 3=apnea only
settingsDefault.PlotEventData2Option=3;
settingsDefault.savelongbreathtables=1;
settingsDefault.BreathEnsembleMethod=1;
settingsDefault.selectstateEventAnalysis=8; % query: confirm this as different to selectstate above, yes, different.
% should we also have selectposEventAnalysis ?
settingsDefault.UseAutoRespInEventAnalysis=0; %set this to decide which event set to use [0,1,2,3,4,5], used in EventAnalysis and in getData currently.
settingsDefault.CalculateDesatAndObstructionSeverity = 1; % EAS added on 2022-01-21 to calculate Desaturation Severity, Obstruction Severity and other traditional parameters. Parameters used by UEF group.
settingsDefault.CalculateEventAreaDepthDuration = 1; % added by EAS on 2022-02-18



