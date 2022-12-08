function getFOT_manualselect()
% global settings AMasterSpreadsheet
    selpath = uigetdir();
    dirx = dir(selpath);
    dirx(1:2)=[];
    %settings = ImportSettings(settings,AMasterSpreadsheet);
    
    for i=1:length(dirx)
        try
            i
        filedir = [selpath '\' dirx(i).name]
        matObj = matfile(filedir);
        varlist = who(matObj);
        load(filedir);
        
%         settings.Fs = 1./(DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1));
        figure(1);
        [DataEventHypnog_Mat,ChannelsList,ChannelsFs,Info.FOTinfo]=FOTfromXHz(DataEventHypnog_Mat,ChannelsList,ChannelsFs,1);
        %[DataEventHypnog_Mat,ChannelsList,ChannelsFs]=runWSanalysis(DataEventHypnog_Mat,ChannelsList,ChannelsFs);
        close all;
        ChannelsFs(length(ChannelsList)+1:end)=[]; %correction of historical error
        if 1 %%enable when ready
            save([selpath '\' dirx(i).name],varlist{:},'-v7.3')
        else
            disp('Saving disabled');
        end
        
        disp('success')
        catch me
            disp(me.message);
            i
            disp('failure')
        end
    end
end