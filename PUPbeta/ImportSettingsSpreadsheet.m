function settings = ImportSettingsSpreadsheet(AMasterSpreadsheet)

try % DLM added this try catch because sheetnames was throwing error
    sheets = sheetnames(AMasterSpreadsheet); % turns out, sheetnames was only introduced in R2019b, hence fails for me running 2018b
catch
    [status,sheets] = xlsfinfo(AMasterSpreadsheet); % xlsinfo only needs R2006a, or later. (Great success.)
end

    if length(sheets)<2
        settings=struct([]);
        return
    end
        
    %%%%%%%%%%% put in settings = ImportSettingsConvert(settings,AMasterSpreadsheet): %%%%%%%%%%%%%%%
    disp('Importing settings from spreadsheet');
    [~,~,options] = xlsread(AMasterSpreadsheet,2,'C3:C11');
    
    settings.Fs = options{1};
    settings.ConvertedDirectory = char(options(2));
    if settings.ConvertedDirectory(end)~=filesep
        settings.ConvertedDirectory=[settings.ConvertedDirectory filesep];
    end
    
    settings.scoredarousalsinwake = options{3};
    % options{4} was EverySecondEEGisRef, but this is per file now, so not used
    % dlm reused options{4} at one point for Profusion position update, but
    % I don't think this is used anymore either...
    % settings.ProfusionPositionUpdate = options{4};
    settings.RunInReverse = logical(options{5});
    settings.EverySecondEKGisRef = logical(options{6});
    settings.FsLeak = double(options{7});
    settings.skiparousalintensity = options{8};
    settings.AnalyzeSnore = options{9};
    if isnan(settings.AnalyzeSnore)
        settings.AnalyzeSnore = 0;
    end
    
    %analysis settings from spreadsheet
    [~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C58');
    settings.savename = char(raw{1});
    settings.AnalyzedDirectory = char(raw{18});
    if settings.AnalyzedDirectory(end)~=filesep
        settings.AnalyzedDirectory=[settings.AnalyzedDirectory filesep];
    end
    %settings.AnalyzedDirectory = char(raw{18}); %reload/overwrite
    
    settings.Pnasaldownisinsp=logical(raw{2});
    settings.LGfromFlowVersion = char(raw{3});
    settings.sqrt_scaling=logical(raw{4});
    settings.saveplots=logical(raw{5});
    settings.plotfigure=logical(raw{6});
    if settings.saveplots % force plotfigure on if saveplots is on
        settings.plotfigure = true;
    end
    
    settings.usescoredcentralapneas=logical(raw{7});
    settings.eventsarebreathsfullywithinmargins=logical(raw{8});
    settings.havescoredcentralhypops=logical(raw{9});
    settings.manualscoringtouchups=logical(raw{10});
    settings.maxdelaybreaths=raw{11};
    settings.windowlength=raw{12};
%     settings.Fs=raw{13};
    settings.ignoreCPAPdata=logical(raw{14});
    settings.exportresultstoxls=logical(raw{15});
    settings.handlemixedeventsseparately=logical(raw{16});
    settings.longestwakeduration=raw{17};
    settings.OutputDataDirectory = char(raw{18});
    settings.ARmodel=logical(raw{19});
    settings.AnalyzeNREMonly=logical(raw{20});
    settings.findcentralhypopneasandapneas=logical(raw{21});
    settings.PreLowPass=raw{22};
    settings.scalingexponent=raw{23};
    settings.minabsPmaskforCPAPoff=raw{24};
    settings.flowshapesonly=raw{25}; %SS: not boolean
    settings.GetLocalSignals=raw{26};
    settings.SavePerSubject=raw{27};
    settings.WindowStep=raw{28};
    settings.FlowSignal = char(raw{29});
    settings.downsampledFs = raw{30};
    settings.spikedetector = logical(raw{31});
    settings.plotbreathdetectionfigures = logical(raw{32});
    settings.modBB_i_start = logical(raw{33});
    settings.nearerzerocrossings = raw{34};
    settings.PreHighPass=raw{35};
    settings.verbose=logical(raw{36});
    settings.savefigas=raw{37};    % 'saveasTIFF' | 'saveasFIG' | 'saveasPNG'
    settings.AnalyzeSnore = raw{38};
    settings.savelongbreathtables=raw{39};
    
    %SUMMARY ANALYSIS SETTINGS
    [~,~,raw] = xlsread(AMasterSpreadsheet,2,'C61:C75'); %settings all
    
    % Read spreadsheet (position codes worksheet)
    if isfield(settings,'useCentralPositionDatabase') && settings.useCentralPositionDatabase==0
        [~,~,settings.poscodesdatabase] = xlsread(AMasterSpreadsheet,3,'B2:K55');
    end
    
    settings.AnalyzedDirectory=char(raw{1});
    %was called directory
    settings.savename=char(raw{2});
    settings.getCIs=logical(raw{3});
    settings.plotfigs=logical(raw{4});
    settings.selecttimeofnight=logical(raw{5});
    settings.selecttimeofnight_XTiles=double(raw{6});
    settings.selecttimeofnight_NthXTile=double(raw{7});
    settings.selectstate=double(raw{8}); %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
    settings.selectposition = char(raw{9});
    
    
    disp(['Position selected: ', settings.selectposition])
    eval(['settings.SummaryMrange=[' num2str(raw{10}) '];']);
    
    settings.SummaryDirectory=char(raw{11}); %SummaryAnalysisDirectory
    settings.nanifnoarousals=double(raw{12});
    settings.Vpassive1d=raw{15};
    temp=double(raw{13});
        if ~isnan(temp), settings.method=temp; end    
    try settings.comb=eval(raw{14}); end
    

    
