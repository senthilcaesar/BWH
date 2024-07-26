function settingsDefault = ImportSettingsDefault(settings)

%% ChannelCheck
settingsDefault.LabelsOnly=0;

%% Directories
settingsDefault.SourceDirectory = [settings.workdir 'Source\'];
settingsDefault.ConvertedDirectory = [settings.workdir 'Converted\'];
settingsDefault.AnalyzedDirectory = [settings.workdir 'Analyzed\'];
settingsDefault.SummaryDirectory = [settings.workdir 'Summary\'];

%% CONVERT
settingsDefault.saveFlowQc=1;
settingsDefault.RunInReverse = 0; % always
settingsDefault.EverySecondEKGisRef = logical(0);
settingsDefault.FsLeak = 0;
settingsDefault.skiparousalintensity = 0;
settingsDefault.AnalyzeSnore = 0;
settingsDefault.scoredarousalsinwake = 0;
settingsDefault.ProfusionPositionUpdate = 0; %Must stay in the off position because otherwise the filename.txt naming convention used by BE/DM will conflict with other systems.
settingsDefault.processEEG = 1; % always
settingsDefault.allpowers = 1; % always
settingsDefault.EKGplot = 0; % not important
settingsDefault.WSanalysis = 1; %converted
settingsDefault.WSArVersion=2;
settingsDefault.findStableBreathing = 1; % this finds stable breathing periods and adds to DataHypnogMat
settingsDefault.runarousalintensity = [1 0 0]; %leave secondary arousal intensity analysis off
settingsDefault.phasicREMfromEOG=1;

%% ANALYSIS
settingsDefault.savename = 'PUP';
% settingsDefault.AnalyzedDirectory = char(raw{18});
% settingsDefault.Pnasaldownisinsp=1; %Default but could be wrong!!
settingsDefault.Fs=125;
settingsDefault.LGfromFlowVersion = 'LGfromFlowBeta';
settingsDefault.sqrt_scaling=1; 
settingsDefault.saveplots=logical(0);
settingsDefault.plotfigure=logical(0);
if settingsDefault.saveplots % force plotfigure on if saveplots is on
    settingsDefault.plotfigure = true;
end
settingsDefault.usescoredcentralapneas=logical(1);
settingsDefault.eventsarebreathsfullywithinmargins=logical(1);
settingsDefault.havescoredcentralhypops=logical(0);
settingsDefault.manualscoringtouchups=logical(0);
settingsDefault.maxdelaybreaths=5;
settingsDefault.windowlength=7;

settingsDefault.ignoreCPAPdata=logical(1);
settingsDefault.exportresultstoxls=logical(0);
settingsDefault.handlemixedeventsseparately=logical(1);
settingsDefault.longestwakeduration=30;
settingsDefault.ARmodel=logical(1);
settingsDefault.AnalyzeNREMonly=logical(0);
settingsDefault.findcentralhypopneasandapneas=logical(0);
settingsDefault.PreLowPass=0;
settingsDefault.scalingexponent=0.67;
settingsDefault.minabsPmaskforCPAPoff=1.2;
settingsDefault.flowshapesonly=0; %SS: not boolean
settingsDefault.GetLocalSignals='None';
settingsDefault.SavePerSubject=1;
settingsDefault.WindowStep=120;
settingsDefault.FlowSignal = 'Normal';
settingsDefault.downsampledFs = 25;
settingsDefault.spikedetector = logical(0);
settingsDefault.plotbreathdetectionfigures = logical(0);
settingsDefault.modBB_i_start = logical(1);
settingsDefault.nearerzerocrossings = 0;
settingsDefault.PreHighPass=0;
settingsDefault.verbose=logical(1);
settingsDefault.savefigas='saveasPNG';    % 'saveasTIFF' | 'saveasFIG' | 'saveasPNG'
settingsDefault.AnalyzeSnore = 0;
settingsDefault.savelongbreathtables = 1;

settingsDefault.plotfiguresqrtscaling=0;
settingsDefault.plotfigureLGparameters=0;

settingsDefault.AutoScoreRespEvents=1; %turned on by default 8/20/21
settingsDefault.UseAutoScoredRespEventsForLG=0; %only influences codeflow if AutoScoredRespEvents==1

settingsDefault.parallelAnalysis = 0;

settingsDefault.minEventDuration = 10;
settingsDefault.ApplyClippingCorrection = 1; %turned on by default 8/20/21

%% Other settings that are intentionally not defined by default, but listed here
%settingsDefault.lowVEareObstructiveEvents = 0.5; %will force small breaths to be scored as "obstructed" (unknown drive) for phenotyping
    %use if using manual scoring that might be unreliable

%% SUMMARY ANALYSIS

%settingsDefault.Mrange=1; %single patient; 
%removed by SS 7/12/21 because this clashes %with spreadhseet read

%settingsDefault.SummaryDirectory=char(raw{11}); %SummaryAnalysisDirectory
%settings.AnalyzedDirectory=char(raw{1});
%settings.directory=settings.AnalyzedDirectory;

settingsDefault.getCIs=logical(0);
settingsDefault.plotfigs=logical(1);
settingsDefault.selecttimeofnight=logical(0);
settingsDefault.selecttimeofnight_XTiles=double(2);
settingsDefault.selecttimeofnight_NthXTile=double(2);
settingsDefault.selectstate=double(4); %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
settingsDefault.selectposition = 'All';

settingsDefault.nanifnoarousals=double(1);
settingsDefault.Vpassive1d=double(1);
settingsDefault.method = 1;
settingsDefault.comb=[0 1 30]; %disabled by default: [enable, odd1even2, combwidth(min)]

settingsDefault.minNeventsLG = 1; 
settingsDefault.minNeventsUA = 1; 
settingsDefault.maxwakethresLG = 30; 
settingsDefault.maxwakethresUA = 300; % more tolerant to more wake for UA measures

settingsDefault.Nciles = 10;%default value, alternative--for smoother data--is 131

settingsDefault.runcompare = 0;
settingsDefault.Comparename = 'CompareXY';
settingsDefault.constantEupnea=0; %1
settingsDefault.normalizeusingconstantEupnea=0; %1
settingsDefault.runcompare2 = 0;

settingsDefault.minNeventsLG = 1;
settingsDefault.minNeventsUA = 1; %used to be zero, not clear when it was changed to 1 (before July 2019).
settingsDefault.maxwakeLG = 30;
settingsDefault.maxwakeUA = 300;

%settings.scoredarousalsinwake = 1;
settingsDefault.useBDT = 1;
settingsDefault.breathlevelinclusion = 0;
settingsDefault.DriveincmH2OnotEupnea = 0;

settingsDefault.drivecol=19;
settingsDefault.minwindows=3;
settingsDefault.EdimtaCalibration=0; %0 default; set to 1 for unbiased comparison between states
settingsDefault.EdiAutoLinearization=1;


%% EVENT ANALYSIS EXTRA SETTINGS
settingsDefault.selecteventtype=1; %1=all 2=hyp only 3=apnea only
settingsDefault.PlotEventData2Option=5;
settingsDefault.savelongbreathtables=1;
settingsDefault.BreathEnsembleMethod=1;
settingsDefault.selectstateEventAnalysis=8;


