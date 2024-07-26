function hackAddStableBreathingToConverted_XHz_FixLoggedErrors()
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

selpath = uigetdir(); 
load([selpath, '\addSB_errorlog.mat']);
dirx = errorlog(:,1);
errorlog2=[];
for i=1:length(dirx)
    try
        disp('.'); % add line space in command window for ease of viewing
        filedir = [selpath '\' char(dirx(i))];
        matObj = matfile(filedir); varlist = who(matObj);
        % before loading the big parts of the file, test here if SB already exists
        ChList = load(filedir, 'ChannelsList');
        if isempty(find(strcmp(ChList.ChannelsList,'StableBreathing')==1))
            str = ['File: ', char(dirx(i)), ', loading ...'];  disp(str);
            load(filedir);
            
            try
                [StableBB_mins, StableBB_Table] = getStableBreathingPeriods(DataEventHypnog_Mat,ChannelsList);
                DataEventHypnog_Mat(:,end+1)=StableBB_mins;
                ChannelsList{end+1} = 'StableBreathing';
                ChannelsFs(end+1) = 1./(DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1));
                Info.StableBB = StableBB_Table; % add to Info struct if exists, otherwise, just make it now.
                save([selpath '\' char(dirx(i))],varlist{:},'-v7.3');
                str = ['File: ', char(dirx(i)), ', great success']; disp(str);
            
            catch FaultFindingorSavingStableBB
                str = ['File: ', char(dirx(i)), ', error finding or saving stable breathing data'];  disp(str);
                % disp(FaultFindingorSavingStableBB.message);
                errorlog2{end+1,1} = char(dirx(i));
                errorlog2{end,2} = FaultFindingorSavingStableBB.message;
            end
            
        else
            str = ['File: ', char(dirx(i)), ', already contains stable breathing data']; disp(str);
        end
        
    catch FaultLoadingFile        
        str = ['File: ', char(dirx(i)), ', could not be loaded']; disp(str);
        % disp(FaultLoadingFile.message);
        errorlog2{end+1,1} = char(dirx(i));
        errorlog2{end,2} = FaultLoadingFile.message;
    end
end

save([selpath, '\addSB_errorlog2.mat'],'errorlog2');

end
