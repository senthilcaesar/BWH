function LoadAndPlot(M,zoomsize,PlotLayout)
    global settings
    ConvertedFileDir = string(strcat(settings.ConvertedDirectory,settings.patients(:,1),'.mat'));
    AnalyzedFileDir = string(strcat(settings.AnalyzedDirectory,settings.savename,'_',string([1:M]'),'.mat'));
    load(ConvertedFileDir(M));
    try
        load(AnalyzedFileDir(M),'SigT2');
    end
    if ~exist('PlotLayout')
        PlotLayout=13;
    end
    PlotXHzData;
    if exist('zoomsize') && ~isempty(zoomsize)
        StartTime = SigT.Time(1);
        Duration = SigT.Time(end) - StartTime;
        xlim(StartTime + Duration/2 + zoomsize*[-0.5 0.5]*60); %may not work when plotting in datetime mode
    end
end
