function hackAddStableBreathingToConverted_XHz_manualselect()
% The purpose of this function is to add stable breathing periods
% to already converted data (i.e. Converted_XHz)
% It adds a column at the end of the DataEventHypnogMat that describes the
% number of minutes of continuous stable breathing (0-10mins)
%
% the file handling process in this fn is based on "getWSanalysis_manualselect"
% Manually select the folder, and then try to process everything in the folder
%
close
clear
clc

errorlog = []; % keeping a record of files that fail, and cause of failure

selpath = uigetdir(); dirx = dir(selpath); dirx(1:2)=[];
for i=1:length(dirx)
    try
        disp('.'); % add line space in command window for ease of viewing
        filedir = [selpath '\' dirx(i).name];
        matObj = matfile(filedir); varlist = who(matObj);
        % before loading the big parts of the file, test here if SB already exists
        ChList = load(filedir, 'ChannelsList');
        if isempty(find(strcmp(ChList.ChannelsList,'StableBreathing')==1))
            str = ['File: ', dirx(i).name, ', loading ...'];  disp(str);
            load(filedir);
            
            try
                [StableBB_mins, StableBB_Table] = getStableBreathingPeriods(DataEventHypnog_Mat,ChannelsList);
                DataEventHypnog_Mat(:,end+1)=StableBB_mins;
                ChannelsList{end+1} = 'StableBreathing';
                ChannelsFs(end+1) = 1./(DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1));
                Info.StableBB = StableBB_Table; % add to Info struct if exists, otherwise, just make it now.
                save([selpath '\' dirx(i).name],varlist{:},'-v7.3');
                str = ['File: ', dirx(i).name, ', great success']; disp(str);
            
            catch FaultFindingorSavingStableBB
                str = ['File: ', dirx(i).name, ', error finding or saving stable breathing data'];  disp(str);
                % disp(FaultFindingorSavingStableBB.message);
                errorlog{end+1,1} = dirx(i).name;
                errorlog{end,2} = FaultFindingorSavingStableBB.message;
            end
            
        else
            str = ['File: ', dirx(i).name, ', already contains stable breathing data']; disp(str);
        end
        
    catch FaultLoadingFile        
        str = ['File: ', dirx(i).name, ', could not be loaded']; disp(str);
        % disp(FaultLoadingFile.message);
        errorlog{end+1,1} = dirx(i).name;
        errorlog{end,2} = FaultLoadingFile.message;
    end
end

save([selpath, '\addSB_errorlog.mat'],'errorlog');

end
