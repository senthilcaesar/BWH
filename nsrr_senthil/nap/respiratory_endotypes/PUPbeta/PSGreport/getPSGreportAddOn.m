
%% run StartHere

%%
settings = ImportSettings(settings,AMasterSpreadsheet);

%what patient do we want, select this:
n=1;
fname = extractBefore(settings.patients{n,1},'_XHz');

if ~isfield(settings,'positioncodes')
    % get Position code
    settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
    settings.supinepositioncode = settings.positioncodes(1);
end

load([settings.ConvertedDirectory settings.patients{n,1} '.mat'])
A_Summary_T=getPSGreport(DataEventHypnog_Mat,Evts,ChannelsList,settings,fname);

