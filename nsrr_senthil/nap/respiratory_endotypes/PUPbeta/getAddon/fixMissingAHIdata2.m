function fixMissingAHIdata2(Mrange)
global settings AMasterSpreadsheet


    settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
    
M = length(settings.patients(:,1));

SourceDir=strcat(string(settings.Filenames(:,4)),string(settings.Filenames(:,1)));
SourceDirE=strcat(string(settings.Filenames(:,5)),string(settings.Filenames(:,2)));
SourceDirH=strcat(string(settings.Filenames(:,6)),string(settings.Filenames(:,3)));
ConvertedFileDir = string(strcat(settings.ConvertedDirectory,settings.patients(:,1),'.mat'));
AnalyzedFileDir = string(strcat(settings.AnalyzedDirectory,settings.savename,'_',string([1:M]'),'.mat'));

for n=Mrange
    
    try
disp(n)
    try
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
    catch
        % above code failing if settings.protocol is char array.
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(n,:));
    end
    
    settings.supinepositioncode = settings.positioncodes(1);
    
    
    C = load(ConvertedFileDir(n),'SigT');
    A = load(AnalyzedFileDir(n),'Evts');
    %SigT2 = load(AnalyzedFileDir(n),'SigT2');
%add AHIdata2 to Evts.EvtsAutoRespOnly

%%
%A.Evts.EvtsAutoRespOnly.RespT;
    if ~isfield(A.Evts.EvtsAutoRespOnly,'AHIdata2')
            ROImask = 1 + 0*C.SigT.Time; %update
            [A.Evts.EvtsAutoRespOnly]=getAHIAll(C.SigT,ROImask,A.Evts.EvtsAutoRespOnly,C.SigT.Properties.VariableNames,settings,'EventsRespAuto','EventsAr');
            try A.Evts.EvtsAutoRespOnly=getEvtPRMainFn(A.Evts.EvtsAutoRespOnly,C.SigT); end
            
            
       %save(AnalyzedFileDir(n),'-append');
       save(AnalyzedFileDir(n),'-struct','A','-append');
    end
    
    %[EventsRespAuto]=FindApneaHypopnea(SigT2.SigT2.VEpeupnea,C.SigT,settings.Fs,C.SigT.Properties.VariableNames,6);
    disp(['finished ' num2str(n)]);
    catch
    disp(['failed ' num2str(n)]);
    end
    
end
            