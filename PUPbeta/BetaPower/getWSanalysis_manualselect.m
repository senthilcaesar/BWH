function getWSanalysis_manualselect
%first run starthere.m
global settings AMasterSpreadsheet
    selpath = uigetdir();
    dirx = dir(selpath);
    dirx(1:2)=[];
    %settings = ImportSettings(settings,AMasterSpreadsheet);
    settings.scoredarousalsinwake = input('Enter scoredarousalsinwake [0 or 1], typical clinical scoring is 0: ');
    for i=1:length(dirx)
        try
        filedir = [selpath '\' dirx(i).name]
        matObj = matfile(filedir);
        varlist = who(matObj);
        load(filedir);
        settings.Fs = 1./(DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1));
        figure(1);
        [DataEventHypnog_Mat,ChannelsList,ChannelsFs]=runWSanalysis(DataEventHypnog_Mat,ChannelsList,ChannelsFs);
        close all;
        ChannelsFs(length(ChannelsList)+1:end)=[]; %correction of historical error
        save([selpath '\' dirx(i).name],varlist{:},'-v7.3')
            i
        'success'
        catch me
            disp(me.message);
            i
            'failure'
        end
    end
end