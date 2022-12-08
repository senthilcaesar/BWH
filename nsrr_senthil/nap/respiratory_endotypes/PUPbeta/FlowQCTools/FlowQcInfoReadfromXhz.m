function FlowQcInfoReadfromXhz(MrangeOverride)
% Run start here first
global AMasterSpreadsheet settings
settings.plotfigure=0;

[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
Filenames = MasterWorksheet(:,11:12);

I=find(cell2mat(cellfun(@(x)any(~isnan(x)),Filenames(:,1),'UniformOutput',false)));
Filenames((I(end)+1):end,:)=[]; %remove data below last legit entry

Mrange=1:size(Filenames,1);
if exist('MrangeOverride')
    Mrange = MrangeOverride;
end

FlowQcT=[];
for ii=Mrange
    ii
    clear Info LowPassPredict HighPassPredict temp temp1 SigT
    try
        Convname=[[settings.workdir,'Converted\'] Filenames{ii,1} '.mat']
        temp1=load(Convname);
        if isfield(temp1, 'DataEventHypnog_Mat') % convert DataEventHypnog_Mat to SigT
            SigT=array2table(temp1.DataEventHypnog_Mat);
            SigT.Properties.VariableNames = temp1.ChannelsList;
            temp1 = rmfield(temp1, {'DataEventHypnog_Mat','ChannelsList'});
        else
            SigT=temp1.SigT;
        end
        
        try
            Info=temp1.Info;
        end
        
        if ~(exist('Info')==1)
            Info=struct();
        end
        if isfield(Info,'FlowQ') && isfield(Info.FlowQ,'IEratio')
            LowPassPredict=Info.FlowQ.FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
            HighPassPredict=Info.FlowQ.FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
            Info.FlowQ=rmfield(Info.FlowQ,'FlowFilterDetect'); % remove for concatenation as a table
            temp=[table(string(Filenames{ii,1})) struct2table(Info.FlowQ) table(LowPassPredict,HighPassPredict,ii)];
            FlowQcT=[FlowQcT;temp];
        else
            clear Flow Time Fs_Flow StartTime N_Flow TimeFlow
            ChannelsFs=temp1.ChannelsFs;
             Flow=SigT.Flow;
             Time=SigT.Time;
             Fs_Flow=ChannelsFs(find(strcmp(SigT.Properties.VariableNames,'Flow')));
             StartTime=SigT.Time(1);
             N_Flow=size(SigT.Flow,1);
             TimeFlow=(StartTime:(1/Fs_Flow):StartTime+(N_Flow-1)*(1/Fs_Flow))'; % This is the time vector associated with the _XHz Flow data.

            % 1. find if flow is inverted & overall noise level
            clear PrUpright FnoiseAll 
            [PrUpright,FnoiseAll,Info]=FlowInvertedDetectTool(Flow,TimeFlow,Info);
            Info.FlowQ.PrUpright = PrUpright;
            Info.FlowQ.FnoiseAll = FnoiseAll;
            
            if FnoiseAll>=0.1
                disp('noisy or absent flow signal')
            end
             % 2. find if a high pass/low pass filter is applied to flow signal
             clear FlowRS FlowFilterDetect
             verbose=0;
             ploton=0;
             FlowRS = interp1(TimeFlow,Flow,Time,'linear');
             FlowFilterDetect = FlowFilterDetector(FlowRS,Time,ploton,verbose);
             
             clear FlowRS;
             Info.FlowQ.FlowFilterDetect = FlowFilterDetect;
             
             if FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1)<8
                 disp('Warning: Flow appears smoothed or downsampled');
             end
             
             if FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1)>0.011
                 disp('Warning: Flow appears distorted by baseline adjustment (high pass)');
             end
             % 3. find if flow is clipped or not
             clear FclippedI FclippedE
             [~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[]);
             
             Ttemp.FclippedI=FclippedI;
             Ttemp.FclippedE=FclippedE;
             Ttemp.FclippedTotal=FclippedE+FclippedI;
             
             Info.FlowQ.FclippedTotal = Ttemp.FclippedTotal;
             
             if Ttemp.FclippedTotal>0.005
                 disp('Warning: Flow appears clipped');
                 disp(['   Clipping fraction: ' num2str(100*Ttemp.FclippedTotal,2) ' %']);
                 
             end
             LowPassPredict=Info.FlowQ.FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
             HighPassPredict=Info.FlowQ.FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
             Info.FlowQ=rmfield(Info.FlowQ,'FlowFilterDetect'); % remove for concatenation as a table
             temp=[table(string(Filenames{ii,1})) struct2table(Info.FlowQ) table(LowPassPredict,HighPassPredict,ii)];
             FlowQcT=[FlowQcT;temp];
        
        end
        
    catch
       
        N=(NaN(1,16));
        N=array2table(N);
        temp=[table(string(Filenames{ii,1})) N];
        
        temp.Properties.VariableNames={'Var1','IEratio','IEratio2','IEratioEstimated','IEratioEstimated2','SNRwindow','SNRwindowNoiseMean','FSNRover15','FSNRover20','FSNRover25','FSNRover30','PrUpright','FnoiseAll','FclippedTotal','LowPassPredict','HighPassPredict','ii'};
        FlowQcT=[FlowQcT;temp];
    end
end

save([settings.workdir '\Summary\FlowQcT.mat'],'FlowQcT');
