%% Export Spike Data

    retest=1;
    if retest
        cedpath = [cd '\External\CEDS64ML2017'];
        addpath(cedpath); % so CEDS64LoadLib.m is found
        CEDS64LoadLib(cedpath); % load ceds64int.dll
        %fhand1 = CEDS64Open(file);
    end
    
    Fs=125;
    ChannelsToSave = {'EMGdiX','EMGdiCh','EMGdiBest','EMGdi1','EMGdi2','EMGdi3','EMGdi4','EMGdi5'}';
    ChannelsToSaveNames = {'EMGdi','EMGdiCh','EMGdiBest','EMGdi1','EMGdi2','EMGdi3','EMGdi4','EMGdi5'}';
    
    
  %  EMGdiX(isnan(EMGdiX))=0;
    
    if 1 %base new timebase file on old file
  %  file = 'J:\PEOPLE\FACULTY\SANDS\O2PSG\AdPhenotype\_Scored\DPW2023Adpheno\2023Adpheno_SpikeData - scored.smrx'
        fhand1 = CEDS64Open([path spikefile]);
        tbase = CEDS64TimeBase(fhand1);%CEDS64TimeBase( fhand1 );%%%obtain time divisions from the original file
    else
        tbase = 0.0005;
    end

        file2 = [path 'EdiData.smrx'];
        fhand2 = CEDS64Create(file2,32,2);
        
        CEDS64TimeBase(fhand2,tbase); %% set the time division for new file

    %hardcoded based on other example recordings, works fine
    
    for i=1:length(ChannelsToSave)
        if ~exist(ChannelsToSave{i})
            continue
        end
        
       eval([ChannelsToSave{i} '(isnan(' ChannelsToSave{i} '))=0;']);
         
            ch=i;
        iOk(i,1)=CEDS64SetWaveChan(fhand2,ch,1/(Fs*tbase),9); %%% set channel 1-- 9 is for channel type "waveform"
        % CEDS64WriteWave(fhand2,ch,Channel_7_125Hz,0);
        eval(['CEDS64WriteWave(fhand2,ch,' ChannelsToSave{i} ',0);']);
        iOk(i,2) = CEDS64ChanTitle( fhand2, ch,ChannelsToSaveNames{i});
    end
    disp(iOk==0)
    CEDS64Close(fhand2);
    

disp('Completed Spike export');
% Close Spike
CEDS64CloseAll(); % close all the files
unloadlibrary ceds64int; % unload ceds64int.dll