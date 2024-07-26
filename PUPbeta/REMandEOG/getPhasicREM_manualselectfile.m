function getPhasicREM_manualselectfile
    selpath = uigetfile();
        filedir = [selpath]
        matObj = matfile(filedir);
        varlist = who(matObj);
        load(filedir);
        [DataEventHypnog_Mat,ChannelsList,ChannelsFs]=getPhasicREM(DataEventHypnog_Mat,ChannelsList,ChannelsFs);
        ChannelsFs(length(ChannelsList)+1:end)=[];
        save([selpath],varlist{:},'-v7.3')
end