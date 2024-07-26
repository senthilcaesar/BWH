function getPhasicREM_manualselect_folder
    selpath = uigetdir();
    dirx = dir(selpath);
    dirx(1:2)=[];
    for i=1:length(dirx)
        filedir = [selpath '\' dirx(i).name]
        matObj = matfile(filedir);
        varlist = who(matObj);
        load(filedir);
        [DataEventHypnog_Mat,ChannelsList,ChannelsFs]=getPhasicREM(DataEventHypnog_Mat,ChannelsList,ChannelsFs);
        ChannelsFs(length(ChannelsList)+1:end)=[];
        save([selpath '\' dirx(i).name],varlist{:},'-v7.3')
    end
end