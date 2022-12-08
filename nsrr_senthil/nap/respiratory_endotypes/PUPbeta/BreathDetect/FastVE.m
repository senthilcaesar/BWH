function FastVE(SigT,settings,n)
ploton=0;
settings.sqrt_scaling=2;
settings.modBB_i_start=0; %work on speeding this up, 

%% downsample and invert flow 
%make a function
if (settings.downsampledFs~=0)&&(settings.downsampledFs~=settings.Fs)
    if strcmp(settings.FlowSignal,'SecondPnasal')
        Flow = SigT.Pnasal;
        disp('Using Pnasal signal in this analysis');
    else
        Flow = SigT.Flow;
    end
    displaytext=['Resampling: Flow data from ' num2str(settings.Fs) ' to ' num2str(settings.downsampledFs) 'Hz...(and back)'];
    disp(displaytext);
    TotalTime = (length(Flow)-1)*(1/settings.Fs);
    Time_ = [0:(1/settings.Fs):TotalTime]';
    TimeDS=[0:(1/settings.downsampledFs):TotalTime]';
    FlowDS = interp1(Time_,Flow,TimeDS,'linear'); % downsample
    Flow = interp1(TimeDS,FlowDS,Time_,'linear');  % upsample back to original length.
    Time_=[];
    TimeDS=[];
    FlowDS=[];
    
    % ensure correct length
    mismatch = length(Flow)-length(SigT.Flow);
    if mismatch>=0  % truncate Flow
        Flow = Flow(1:end-mismatch);
    else  % add spacer to Flow
        Flow = [Flow; zeros(abs(mismatch),1)];
    end
    if strcmp(settings.FlowSignal,'SecondPnasal')
        SigT.Pnasal = Flow; %write it back into DataEventHypnog_Mat
    else
        SigT.Flow = Flow;
    end
    
    %Invert flow as needed; was previously inside LGfromFlowPart1
    SigT.Flow=((-1)^settings.Pnasaldownisinsp)*((-1)^settings.invertflowlist(n))*SigT.Flow; %%% added invert flow
    
end

%% Noise Wav
AllowNoisierFlowFactor=[]; % 2021-11-02 10:45:35 Tuesday EAS changed from AllowNoisierFlow to AllowNoisierFlowFactor (bug fix)
if isfield(settings,'AllowNoisierFlow') && settings.AllowNoisierFlow==1
    AllowNoisierFlowFactor=settings.AllowNoisierFlowFactor;
end
noisewav = FlowSignalToNoise(SigT.Time,SigT.Flow,settings.plotfigure,AllowNoisierFlowFactor); 

%% Analysis, Window

numwind=ceil((height(SigT)-settings.windowlength*60*settings.Fs)/(settings.WindowStep*settings.Fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
winnumrange=1:1:numwind;

%%
for winNum=winnumrange
    disp(winNum);
%     try
%         progressbar([],winNum/(length(winnumrange)),[]); % update progress bar
%     catch me
%         me.message
%     end
    
    lefti=(winNum-1)*settings.WindowStep*settings.Fs+1;
    righti=(winNum-1)*settings.Fs*settings.WindowStep+settings.windowlength*60*settings.Fs;
    
    % Dan Vena added an condition to get end index of last partial window
    % 7/22/2021
    if righti > size(SigT,1)
        righti = size(SigT,1);
    end
    
    %hypnogwind=hypnog(lefti:righti);
    %CPAPwind=CPAP(lefti:righti);
    noisewind=noisewav(lefti:righti);
    %findest longest durations of wakefulness
    %clear I I2 I3 I4
    %I=hypnogwind==4; I2=[I(1)==1;diff(I)]; I3=find(I2==1); I4=find(I2==-1);
    %if length(I4)<length(I3), I4(length(I3),1)=length(hypnogwind)+1; end
    %lengthwake=(I4-I3)/settings.Fs;
    %if isempty(lengthwake),lengthwake=0; end
    %LongestWake=max(lengthwake);
    %Position data in window
    %     PositionData(winNum) = mode(Evtsdata{n}.PositionVariable(lefti:righti));
    
    % If it is N1, N2 or N3 sleep the whole window, calculate loop gain:
    %allsleep(winNum)=(sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0);
    %CPAPmedian=median(CPAPwind);
    %CPAPstd=std(CPAPwind);
    %CPAPabs95=max([prctile(CPAPwind,95),-prctile(CPAPwind,5)]);
    %if 1
    %    CPAPoff_(winNum)=CPAPabs95<settings.minabsPmaskforCPAPoff;%CPAPabs95<settings.minabsPmaskforCPAPoff;
    %else
    %    CPAPoff_(winNum)=prctile(CPAPoff(lefti:righti),1); %original method
    %end
    %if round(CPAPmedian*100)/100==-1, CPAPoff=0;end
    %CPAPData(winNum,:)=[CPAPoff_(winNum) CPAPmedian CPAPstd CPAPabs95];
    %Fnoise3=sum(noisewind==3)/length(noisewind);
    Fnoise2(winNum)=sum(noisewind>=2)/length(noisewind);
    %Fnoise1=sum(noisewind>=1)/length(noisewind);
    %StoNData.Fnoise(winNum,:)=[Fnoise1 Fnoise2(winNum) Fnoise3];
    %FNREM1=sum(hypnogwind==2)/length(hypnogwind);
    %FNREM2=sum(hypnogwind==1)/length(hypnogwind);
    %FNREM3=sum(hypnogwind==0)/length(hypnogwind);
    %FREM=sum(hypnogwind==3)/length(hypnogwind);
    %FNREM=FNREM1+FNREM2+FNREM3;
    %FWAKE=sum(hypnogwind==4)/length(hypnogwind);
    %SleepData(winNum,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
    %ContainsUnknownSleep(winNum)=sum(hypnogwind==-1|hypnogwind==8|isnan(hypnogwind))>0;
    
    %REMdetected(winNum)=(sum(hypnogwind==3)>0);
    
    %choose criteria; mostly just run everything here and tidy up in SummaryAnalysis
    criteria(winNum) = Fnoise2(winNum)<0.1;
    if Fnoise2(winNum)>=0.1
        displaytext=['Skipping ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)...
            ', noisy or absent flow signal'];
        disp(displaytext);

    end
    
    AnalysisIndex(winNum,:)=[(winNum-1)*settings.WindowStep*settings.Fs+1 (winNum-1)*settings.Fs*settings.WindowStep+settings.windowlength*60*settings.Fs];
        
    if criteria(winNum)
        displaytext=['Analyzing ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)];
        disp(displaytext);
        %try
        tic
            [BreathDataTable{winNum},Vflow_out{winNum}]= ...
            VEFlowAll(SigT(lefti:righti,:),settings,ploton); %,settings,n,winNum,ChannelsList
%         catch me
%             disp(['error evaluating VEFlowAll or saving its data: ' me.message])
%             BreathDataTable{winNum} = NaN;
%             Vflow_out{winNum} = [];
%         end
        toc
    else
        BreathDataTable{winNum} = NaN;
    end % if criteria for LGfromFlow analysis
	
end



%%
function [bT,Vflow_out] = VEFlowAll(SigT,settings,ploton) 

    Flow=SigT.Flow; %%% added invert flow
    Time=SigT.Time; %%% added invert flo
VI=NaN;

    FlowOrig = Flow;
    %settings.plotfiguresqrtscaling=1;
    %% Estimate Baseline Drift (SS 2/6/2022)
    
    try
        
        medianlength=18;
        medianstepsize=1;
        centileoptions=[40:1:60];
        Noverlaps=round(medianlength*settings.Fs);
        Nstepsize=round(medianstepsize*settings.Fs);
        FlowBuffer = buffer(FlowOrig,Noverlaps,Noverlaps-Nstepsize,'nodelay');
        temp = prctile(FlowBuffer,centileoptions); %moving time centiles
        tempminmax=prctile(FlowBuffer,[0 100]);
        %temp = movmean(temp',3*[1 1])'; %faster
        temp2 = std(temp');
        [~,Itemp] = min(temp2); 
        FlowDrift = movmean(temp(Itemp,:),3*[1 1]);
        TimeDrift = Time(1) + medianlength/2 + [0:medianstepsize:medianstepsize*(length(FlowDrift)-1)]';
        
        FlowDriftB = interp1(TimeDrift,FlowDrift,Time);
            FlowDriftB(Time<TimeDrift(1))=FlowDrift(1);
            FlowDriftB(Time>TimeDrift(end))=FlowDrift(end);
            FlowDriftAmplitudeA = 100*temp2(Itemp)/std(Flow - FlowDriftB); %Per Percentage of Flow SD
            FlowDriftAmplitudeB = 100*temp2(Itemp)/std(tempminmax(1,:) - tempminmax(end,:)); %Per Percentage of Flow Envelope SD
            disp(['Zero Flow Drift Estimation: ' num2str(FlowDriftAmplitudeA,2) '%, ' num2str(FlowDriftAmplitudeB,2) ,'%']);
        if isfield(settings,'DriftEstimation') && settings.DriftEstimation>0 && settings.plotfiguresqrtscaling==1
            %using settings.plotfiguresqrtscaling==1 as a guide
            if ploton
            figure(88); set(gcf,'color',[1 1 1]);
            plot(centileoptions,temp2/std(Flow - FlowDriftB)*100,'-');
            ylabel('Drift SD (%Flow SD)');
            xlabel('Zero Flow %From Exp to Insp');
            
            figure(891); clf(891);
            title('Zero Flow Drift Estimation');
            ax(1)=subplot(2,1,1);
            plot(Time,FlowOrig);
            hold on;
            plot(Time,prctile(Flow,centileoptions(Itemp))+Flow*0,'k-');
            %plot(TimeDrift,temp(1,:));
            %plot(TimeDrift,temp(end,:));
            box off;
            plot(Time,FlowDriftB,'r-');
            ax(1)=subplot(2,1,2);
            plot(Time,Flow - FlowDriftB)
            hold on 
            plot(Time,Flow*0,'k');
            linkaxes(ax(1),'x');
            box off;
            end
        end    
        if settings.DriftEstimation==2
            Flow = Flow-FlowDriftB;
        end
    catch
        FlowDriftAmplitudeA = NaN;
        FlowDriftAmplitudeB = NaN;
    end
    
    
    %% Flow optional additional filters
    
    if settings.PreLowPass>0
        disp('warning: Pre low pass on Flow');
        F_lowcut = settings.PreLowPass;
        filter_order1 = 4;
        [B_butterPre,A_butterPre] = butter(filter_order1,[F_lowcut]/(1/dt/2),'low');
        Flow=filtfilt(B_butterPre,A_butterPre,Flow);
        if 0&&settings.plotfigure
            figure(8);  clf(8);
            plot(Time,FlowOrig,'color',0.8*[1 1 1]); hold('on');
            plot(Time,Flow,'k'); hold('off');
        end
    end
    
    if settings.PreHighPass>0
        disp('warning: Pre high pass on Flow');
        F_highcut = settings.PreHighPass;
        filter_order1 = 4;
        [B_butterPre,A_butterPre] = butter(filter_order1,[F_highcut]/(1/dt/2),'high');
        Flow=filtfilt(B_butterPre,A_butterPre,Flow);
        if 0&&settings.plotfigure
            figure(8);  clf(8);
            plot(Time,FlowOrig,'color',0.8*[1 1 1]); hold('on');
            plot(Time,Flow,'k'); hold('off');
        end
    end
    
    %% Actual breath detection 
    
    bT=table();
    [bT,bSigT] = VEfromFlowT(Time,Flow,settings,ploton);
    Vflow_out = bSigT.Vflow_out;
    
%     if 0
%         figure(789)
%         ax789(1)=subplot(2,1,1);
%         plot(Time,Flow);
%         
%         ax789(3)=subplot(2,1,2);
%         plot(Time,Vflow_out);
%         
%         Vflow_out2 = Flow-leak;
%         exponent = settings.scalingexponent;
%         Vflow_out2(Vflow_out2>0)=(Vflow_out2(Vflow_out2>0).^(exponent))/(IEratio^0.5);
%         Vflow_out2(Vflow_out2<0)=(-(-Vflow_out2(Vflow_out2<0)).^(exponent))*(IEratio^0.5);
%         leaksignalfilt = interp1(Time(BB_i_mid),leak_B,Time,'nearest','extrap');
%         Vflow_out2 = Vflow_out2-leaksignalfilt(:);
%         plot(Time,[Vflow_out Vflow_out2]);
%     end
    
    
    %DM to make breath-by-breath (mean) values of leak and IEratio
    %(repetitive), put into table, and concatenate table with
    %BreathDataTable at the end of function
    
% catch me
%     disp(me.message);
%     
% end







