function ChannelNumbers = getChannelNumbersHarmonizedLabels(EDFfilenamedir,HarmonizedChannelLabels,Label)
try
    if ~exist('Label')
    [Label,~,Fs] = EDFChannelLabels(EDFfilenamedir);
    Label=strtrim(Label);
    end
    HarmonizedChannelT=readtable(HarmonizedChannelLabels);
    SigList=HarmonizedChannelT.Properties.VariableNames;
    ChannelNumbers=NaN(1,size(SigList,2));
    for ii=1:length(SigList)
        
        SigToLoad=HarmonizedChannelT(:,ii);
        loc=cellfun('isempty',SigToLoad{:,SigList{ii}});
        SigToLoad(loc,:)=[];
        clear SigColTemp
        clear ind
        for jj=1:height(SigToLoad)
           
            ind = find(strcmpi(SigToLoad{jj,1},Label));
            if length(ind)>1
                disp(['Warning: Found two signals called: ' char(SigToLoad{jj,1}) '; Chose first one'])
            end
            if ~isempty(ind)
                SigColTemp=ind(1);
                ChannelNumbers(ii)=SigColTemp;
               break
            end
        end
        
        
    end
catch
end
