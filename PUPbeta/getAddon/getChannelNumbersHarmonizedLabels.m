function ChannelNumbers = getChannelNumbersHarmonizedLabels(EDFfilenamedir,HarmonizedChannelLabels)
try
    [Label,~,Fs] = EDFChannelLabels(EDFfilenamedir);
    Label=strtrim(Label);
    HarmonizedChannelT=readtable(HarmonizedChannelLabels);
    SigList=HarmonizedChannelT.Properties.VariableNames;
    ChannelNumbers=NaN(1,size(SigList,2));
    for ii=1:length(SigList)
        
        SigToLoad=HarmonizedChannelT(:,ii);
        loc=cellfun('isempty',SigToLoad{:,SigList{ii}});
        SigToLoad(loc,:)=[];
        clear SigColTemp
        for jj=1:height(SigToLoad)
            ind = find(strcmpi(SigToLoad{jj,1},Label));
            if ~isempty(ind)
                SigColTemp=ind;
                ChannelNumbers(ii)=SigColTemp;
            end
        end
        
        
    end
catch
end
