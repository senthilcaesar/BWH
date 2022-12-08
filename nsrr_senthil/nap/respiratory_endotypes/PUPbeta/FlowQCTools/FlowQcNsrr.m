function FlowQcNsrr(infold,outfold)

% infold= 'U:/EndotypeAnalysis/Source/mesa-sleep-0001-remap.edf';
% outfold= 'U:/EndotypeAnalysis/Source/';
% matlab -nodisplay -r "FlowQcNsrr($infold $outfold)"

%% settings
mydir  = extractBefore(mfilename('fullpath'),'FlowQcNsrr');
addpath(genpath(mydir));

disp (mydir)

settings.Fs=64;
settings.savefigure=0;

%% read edf
ChannelFound = 0;
EDFfilenamedir = infold;
% StrPos=strsplit(EDFfilenamedir,'\');
StrPos=strsplit(EDFfilenamedir,'/');
StrPos=StrPos{end};
ID=string(extractBefore(StrPos,'.edf'));


try
    [Label,~,Fs] = EDFChannelLabels(EDFfilenamedir);
    Label=strtrim(Label);
    SigToLoad={'nas_pres','csCAN','Cannula Flow','Cannulaflow','CannulaFlow','Flow','PNasal','Nasal Pressure','NASAL PRESS','NASAL PRES','NASAL PRESSURE'};
   % use canonical signal output. original sampling freq
    for jj=1:length(SigToLoad)
        ind = find(strcmpi(SigToLoad(jj),Label));
        if ~isempty(ind)
            SigColTemp(jj)=ind;
            
        end
    end
    ChNum=(nonzeros(SigColTemp));
catch
    FlowQcT=table(ID,NaN,NaN,NaN,NaN,NaN,NaN);
end

try
    [Flow,fs_flow,~,~,LabelEDF{1},~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChNum-1,0,Inf);
    ChannelFound=1;
catch me
    ChannelFound=0;
end

if ChannelFound==1
    %% Time information
    N_Flow=length(Flow);
    N_timeXHz = round((N_Flow/fs_flow*settings.Fs));
    fid = fopen(EDFfilenamedir,'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    try
        StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
        if StartTime<43200; StartTime=StartTime+86400; end
    catch me
        StartTime=0;
        disp('Failed to import EDF Start Time, set to 0 sec');
    end
    Time=(StartTime:(1/settings.Fs):StartTime+(N_timeXHz-1)*(1/settings.Fs))'; % This is the time vector associated with the _XHz Flow data.
    EndTime=Time(end);
    
    %If needed for handling flow etc:
    TimeFlow=(StartTime:(1/fs_flow):StartTime+(N_Flow-1)*(1/fs_flow))'; % This is the time vector associated with the _XHz Flow data.
    
    
    %% flow QC
    try
        
        
        clear FlowQcT
        
        % 1. find if flow is inverted & overall noise level
        [prob_upright,fraction_noise]=FlowInvertedDetectTool(Flow,TimeFlow);
        prob_upright=round(prob_upright*100)/100;
        fraction_noise=round(fraction_noise*100)/100;
        
        if fraction_noise>=0.1
            disp('noisy or absent flow signal')
        end
        
        
        % 2. find if a high pass/low pass filter is applied to flow signal
        verbose=0;
        ploton=0;
        FlowRS = interp1(TimeFlow,Flow,Time,'linear');
        FlowFilterDetect = FlowFilterDetector(FlowRS,Time,ploton,verbose);
        clear FlowRS;
        
        filter_lp_predict=FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
        filter_lp_predict=round(filter_lp_predict*100)/100;
        filter_hp_predict=FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
        filter_hp_predict=round(filter_hp_predict*10000)/10000;
        
        
        if FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1)<8
            disp('Warning: Flow appears smoothed or downsampled');
        end
        
        if FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1)>0.011
            disp('Warning: Flow appears distorted by baseline adjustment (high pass)');
        end
        
        % 3. find if flow is clipped or not
        [~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[]);
        fraction_clipping=FclippedE+FclippedI;
        fraction_clipping=round(fraction_clipping*10000)/10000;
        
        if fraction_clipping>0.005
            disp('Warning: Flow appears clipped');
        end
        disp(['   Clipping fraction: ' num2str(100*fraction_clipping,2) ' %']);
        
        FlowQcT = table(ID,fs_flow,prob_upright,fraction_noise,filter_lp_predict,filter_hp_predict,fraction_clipping);
        
    catch
        
        FlowQcT=table(ID,NaN,NaN,NaN,NaN,NaN,NaN);
        FlowQcT.Properties.VariableNames={'ID','fs_flow','prob_upright','fraction_noise','filter_lp_predict','filter_hp_predict',...
            'fraction_clipping'};
    end
    
else
    FlowQcT=table(ID,NaN,NaN,NaN,NaN,NaN,NaN);
    FlowQcT.Properties.VariableNames={'ID','fs_flow','prob_upright','fraction_noise','filter_lp_predict','filter_hp_predict',...
        'fraction_clipping'};
end


% mainfolder2=extractBefore(EDFfilenamedir,id);
mainfolder2=outfold;
fnamensrr=[mainfolder2 'flow_qc.txt'];
writetable(FlowQcT, fnamensrr,'WriteVariableNames',true,'Delimiter',' ')




