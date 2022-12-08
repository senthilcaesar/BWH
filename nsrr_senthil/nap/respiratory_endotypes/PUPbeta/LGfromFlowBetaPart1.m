% This version uses LG calculation software neatened for paper publication....
function [BreathDataTable,BreathFLDataTable,BreathSnoreTable,LocalSignals,Vflow_out] = ...
    LGfromFlowBetaPart1(SigT)
%% global variables
global n winNum settings ChannelsList FlowDriveModels

%% local settings
clocktimestart = clock;
settings.plotfiguresqrtscaling=0;
settings.plotfigureLGparameters=0;
WindowDuration = settings.windowlength*60;

%if sum(n==settings.invertflowlist)         % edited out by dlm
%if ismember(n, settings.invertflowlist)     % replacement code, take one
if settings.invertflowlist(n)               % replacement code, take two
    invertflow=1;
else
    invertflow=0;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if settings.plotfigure==1
    FDfig = 99;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Control parameters for operation of algorithm
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ArousalShift=settings.ARmodel;

if ArousalShift==0 %value of 1 should lead to normal use of VRA; was inverted from ?? until 10/7/19 -SS
    % was correct in 20190205 code version, was incorrect in 20190404 code version and 20190629 until 10/7/19 as above. (DLM)
    VraOn=0;
else
    VraOn=1;
end
% mark_onset_of_arousal_only=0;
delete_arousal_in_events=1;
seteventsbasedonvlessthaneupnea2=7;
eventslessthanX=0.7; %default<0.7
findcentralhypopneasandapneas=settings.findcentralhypopneasandapneas;
if settings.havescoredcentralhypops
    removecentralhypops=0;
else
    removecentralhypops=1;
end

numArousalBreaths=2;
dt=1/settings.Fs;
filter_LFcutoff_butter1 = 0.03;
filter_HFcutoff_butter1 = 10;
filter_order1 = 2;

if 0
[B_butter1,A_butter1] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Extraction and Preprocessing of flow data for calculation of LG:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Flow=((-1)^settings.Pnasaldownisinsp)*((-1)^invertflow)*SigT.Flow; %%% added invert flow

% Oral flow
if strcmp('FlowOral',SigT.Properties.VariableNames)
    FlowOral = ((-1)^invertflow)*SigT.FlowOral;
end

if 0
    Flow_direct = SigT.Flow;
    figure(965); clf(figure(965));
    plot(Flow_direct); hold on
    plot(Flow);
    legend('Direct','Processed');
end

PosChan = find(strcmp(ChannelsList,'Position')==1);
if ~isempty(PosChan)
    PositionSig = 0.1*round(10*SigT.Position);
else
    PositionSig = NaN*SigT{:,1};
end

Time=SigT.Time;
if isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==1
    Arousal=SigT.EventsArWS; %EEG selection uses scoring
elseif isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==2
    Arousal=SigT.EventsArWSB; %EEG selection does not use scoring
else
    Arousal=SigT.EventsAr; %Use original arousal scoring
end

hypnog=SigT.Epochs;

spo2Chan=find(strcmp(ChannelsList,'SpO2')==1);
if ~isempty(spo2Chan)
    spo2=SigT.SpO2;
    %     filter_HFcutoff_butter0 = 2; filter_order0 = 2;
    %     [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    %     spo2 = filtfilt(B_butter0,A_butter0,spo2);
    %     spo2 = round(spo2);
    %     spo2(spo2>100)=100;
    %     spo2(spo2<45)=NaN;
else
    spo2 = NaN*SigT{:,1};
end
%
% % Load in the manually scored events.
% %obstructive if obstrucive apnea / hypopnea / or mixed apnea
% EventsChan = find(strcmp(ChannelsList,'EventsResp')==1);
% Events=sum(DataEventHypnog_Mat(:,EventsChan)'==[2 4 5]')';
% EventsScoredCentralApnea=1*[DataEventHypnog_Mat(:,EventsChan)==3];
% EventsScoredCentralHyp=1*[DataEventHypnog_Mat(:,EventsChan)==6];
%
% %first half of mixed event is central apnea
% if settings.handlemixedeventsseparately
%     try
%         EventsMixed = 1*[DataEventHypnog_Mat(:,EventsChan)==5];
%         if sum(EventsMixed)>0
%             Imixed = find([NaN;diff(EventsMixed)]==1); %index of arousal onset
%             Imixed2 = find(diff(EventsMixed)==-1); %index of pre arousal onset
%             if Imixed2(1)<Imixed(1) %tidy up in case they don't match
%                 Imixed = [1;Imixed];
%             end
%             if Imixed2(end)<Imixed(end)
%                 Imixed2 = [Imixed2;length(EventsMixed)];
%             end
%             lengthmixedE=Imixed2-Imixed+1;
%             for i=1:length(Imixed)
%                 ImixedMID = round(lengthmixedE(i)/2)+Imixed(i);
%                 EventsScoredCentralApnea(Imixed(i):ImixedMID)=1;
%                 Events(Imixed(i):ImixedMID)=0;
%             end
%         end
%     catch me
%     end
% end


if length(Flow)==1
    MiscEpochData=NaN;
    BreathDataTable = NaN;
    BreathFLDataTable = NaN;
    BreathSnoreTable = NaN;
    LocalSignals = NaN;
    disp('exiting 1');
    return
end

%% RIP analysis
% included as (1) placeholder for future development, (2) interesting, (3) cool plot
% assign variables here, they are used below, even if not here for TAA
try
Ab = SigT.Abdomen;
Th = SigT.Thorax;
catch
Ab = [];
Th = [];
end
%if there are no Th or Ab data the channel will be set to zeros
if isempty(Th)
    Th = 0*Time;
end
if isempty(Ab)
    Ab = 0*Time;
end

%% RIP TAA analysis by DM
% note this using parallel processing toolbox, and takes a significant
% amount of time to run for an entire study (i.e. a couple of minutes)
TAA_=NaN*Time';
if 0
    try
        TAA_ = CalcTAA([SigT.Time,Ab,Th]);
        if 0&&settings.plotfigure % TAA figure
            T_ = nanmedian(Th);
            A_ = nanmedian(Ab);
            figure(9); clf(figure(9));
            ax(1)=subplot(3,1,1);
            plot(Time, Flow); % Flow for plot only
            ylabel('Flow');
            ax(2)=subplot(3,1,2);
            plot(Time,Th-T_,'r-',...
                Time,Ab-A_,'b-');
            ylabel('RIPS');
            ax(3) = subplot(3,1,3);
            TAA_(isnan(TAA_))=0;
            plot(Time, TAA_);
            ylim([-1.1 1.1]);
            ylabel('TAA');
            linkaxes(ax, 'x');
        end
    catch me
        disp('failed RIP TAA analysis');
    end
end

%% Flow analysis

%leak=NaN;

clocktimestart2 = clock; % delta in seconds

VI=NaN;
try
    FlowOrig = Flow;
    
    %Estimate Baseline Drift (SS 2/6/2022)
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
            
            box off;
        end    
        if settings.DriftEstimation==2
            Flow = Flow-FlowDriftB;
        end
    catch
        FlowDriftAmplitudeA = NaN;
        FlowDriftAmplitudeB = NaN;
    end
    
    if (isfield(settings, 'ManuallyEditBreathTimes') && settings.ManuallyEditBreathTimes)
        Snore = SigT.Snore;
        SnoreDB = SigT.SnoreDB;
    end
    
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
    
    % ToDo: Vflow_out is shifted, and is causing odd breath starts.
    
    [time,Vflow,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot_B,leak,...
        IEratio,VT,Vpeak,Vpeakmean,Apnea_B,Vflow_out,VTi,VTe,Ti,Te,leak_B,IEratioEstimated] =...
        VEfromFlow(Time,Flow);
    
    MiscEpochData=[leak];
    
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
    
catch me
    disp(me.message);
    
end


if length(VI)==1
    MiscEpochData=NaN;
    BreathDataTable = NaN;
    BreathFLDataTable = NaN;
    BreathSnoreTable = NaN;
    LocalSignals = NaN;
    Vflow_out = NaN.*Time;
    disp('exiting 2: problem with breath detection analysis');
    return
end

if settings.verbose==2
    delta_t = etime(clock, clocktimestart); % delta in seconds
    str = ['Analysis duration (PostFlow) : ', num2str(round(delta_t,1)), ' s']; disp(str);
end

deltatimeFlow = etime(clock, clocktimestart2); % delta in seconds
disp(['BreathDetectionTimeElapsed: ' num2str(deltatimeFlow,3) ' s']);

%% Intentional Clipping - for test purposes only
ClipThresholdFmax=0.90;
ClipFthresholdFnormal=0.002; %higher value removes more (i.e. false) clipping (0.002)

if isfield(settings,'intentionalclipping')==1 && settings.intentionalclipping==1
    
    baseline = median(FlowOrig(BB_i_start));
    fractionforclipping = 2;
    
    for i = 1:length(BB_i_start)
        PeakExp(i) = baseline - min(FlowOrig(BB_i_start(i):BB_i_end(i)));
        PeakInsp(i) = max(FlowOrig(BB_i_start(i):BB_i_end(i))) - baseline;
    end
    
    temp = (PeakInsp + PeakExp) /2;
    I = temp<prctile(temp,95) & temp>prctile(temp,5);
    temp = temp(I);
    peakflowtypical = mean(temp.^settings.scalingexponent)^(1/settings.scalingexponent);
    upperline = fractionforclipping*peakflowtypical + baseline;
    lowerline = baseline - fractionforclipping*peakflowtypical;
    
    Flowclipped = FlowOrig;
    Flowclipped(Flowclipped>upperline)=upperline;
    Flowclipped(Flowclipped<lowerline)=lowerline;
    
    
    for i = 1:length(BB_i_start)
        tempI = FlowOrig(BB_i_start(i):BB_i_mid(i)-1);
        tempE = FlowOrig(BB_i_mid(i):BB_i_end(i)-1);
        FclippedUpper(i) = sum(tempI >= upperline)/length(tempI);
        FclippedLower(i) = sum(tempE <= lowerline)/length(tempE);
    end
    
    %Inoclipping = FclippedUpper==0 & FclippedLower==0;
    
    [time,Vflowclipped,BB_i_startclipped,BB_i_midclipped,BB_i_endclipped,...
        BB_tclipped,VIclipped,VEclipped,Ttot_Bclipped,leakclipped,...
        IEratioclipped,VTclipped,Vpeakclipped,Vpeakmeanclipped,...
        Apnea_Bclipped,Vflow_outclipped,VTiclipped,VTeclipped] =...
        VEfromFlow(Time,Flowclipped);
    
    if 0 %testing effect of gentle filtering on clipping detection.
        filter_order1 = 2;
        [B_butterPre,A_butterPre] = butter(filter_order1,[15]/(1/dt/2),'low');
        Flowclipped=filtfilt(B_butterPre,A_butterPre,Flowclipped);
    end
    
    if 1
        [FclippedIB,FclippedEB,FclippedI,FclippedE,~]= ClipDetection(Flowclipped,...
            BB_i_startclipped,BB_i_midclipped,BB_i_endclipped,ClipThresholdFmax,ClipFthresholdFnormal,1);
        
        %testing specificity for clipping by passing in unclipped signal
        [~,~,FclippedI_,FclippedE_,~]= ClipDetection(FlowOrig,...
            BB_i_startclipped,BB_i_midclipped, BB_i_endclipped,ClipThresholdFmax,ClipFthresholdFnormal,0);
        
        detectedClippingInOriginalData = 100*[FclippedI_ FclippedE_];
    end
    
    if 0
        figure(110); clf(110);
        subplot(2,1,1)
        plot(Time,FlowOrig);
        hold on
        plot(Time,lowerline + 0*Time);
        plot(Time,upperline + 0*Time);
        plot(Time,Flowclipped);
        
        FlowClipped2 = Flowclipped;
        FlowClipped2(Flowclipping_clipped99==0)=NaN;
        plot(Time,FlowClipped2,'r');
        
        subplot(2,1,2)
        plot(Time,Vflow);
        hold on
        plot(Time,Vflowclipped);
    end
    
    % Align clipped breaths with non-clipped breaths
    
    temp = NaN*zeros(1,length(BB_i_start));
    VTi_clipped = temp;
    VTe_clipped = temp;
    sumerrI_t = temp;
    sumerrE_t = temp;
    Ib = temp;
    for i = 1:length(BB_i_start)
        [~, Ib(i)] = min(abs(BB_i_midclipped(:)-BB_i_mid(i)));
        sumerrI_t(i) = dt*(abs(BB_i_start(i) - BB_i_startclipped(Ib(i))) + abs(BB_i_mid(i) - BB_i_midclipped(Ib(i))));
        sumerrE_t(i) = dt*(abs(BB_i_end(i) - BB_i_endclipped(Ib(i))) + abs(BB_i_mid(i) - BB_i_midclipped(Ib(i))));
        VTi_clipped(i) = VTiclipped(Ib(i));
        VTe_clipped(i) = VTeclipped(Ib(i));
    end
    
    IinsErr = sumerrI_t>0.1;
    IexpErr = sumerrE_t>0.1;
    %temp = [VTi(:) VTi_clipped(:) Fclipped_IB(:) Fclipped_EB(:) sumerrI_t sumerrE_t];
    VTi_clipped(IinsErr==1)=NaN;
    VTe_clipped(IexpErr==1)=NaN;
    InoclippingI = FclippedUpper==0;
    InoclippingE = FclippedLower==0;
    VTiratioNoClip = VTi_clipped(InoclippingI)./VTi(InoclippingI);
    VTeratioNoClip = VTe_clipped(InoclippingE)./VTe(InoclippingE);
    
    VTiratioNoClip=VTiratioNoClip(:);
    VTeratioNoClip=VTeratioNoClip(:);
    
    VTe_clipped = VTe_clipped/nanmedian(VTeratioNoClip);
    VTi_clipped = VTi_clipped/nanmedian(VTiratioNoClip);
    
    FclippedUpper = FclippedUpper(:);
    FclippedLower = FclippedLower(:);
    
    FVTiclippedout = 1 - VTi_clipped(:)./VTi(:);
    FVTeclippedout = 1 - VTe_clipped(:)./VTe(:);
    
    VTiCorrection = VTi(:)./VTi_clipped(:);
    VTeCorrection = VTe(:)./VTe_clipped(:);
    
    
    BreathDataTableClipped = table(FclippedUpper,FclippedLower,VTiCorrection,VTeCorrection);
    
    if 1
        figure(998); clf(998); set(gcf,'color',[1 1 1]);
        plot(FclippedUpper,VTiCorrection,'.')
        hold on
        plot(FclippedLower,VTeCorrection,'.')
        
        Fline = 0:0.01:0.9;
        
        Yline = 1./( cos(0.5*pi()*Fline).^0.33);
        
        plot(Fline,Yline,'r')
        
        hold off
        box off
        
    end
end

%% Clipping detection (detects percentage clipped per window and per breath)

%  [percentageclippedInsp95, percentageclippedExp95, NumberClipped95, Flowclipping95,RobustFlowclipping95, ClipFractionBreath95] = ClipDetection(Flow, BB_i_start, BB_i_end,0.95);
%  [percentageclippedInsp99, percentageclippedExp99, NumberClipped99, Flowclipping99,RobustFlowclipping99,ClipFractionBreath99] = ClipDetection(Flow, BB_i_start, BB_i_end,0.99);
%  [percentageclippedInsp, percentageclippedExp, NumberClipped, Flowclipping,RobustFlowclipping,ClipFractionBreath] = ClipDetection(Flow, BB_i_start, BB_i_end,0.99);

if isfield(settings,'intentionalclipping')==1 && settings.intentionalclipping==1
    [FclipIB,FclipEB,FclipI,FclipE,~]= ClipDetection(Flowclipped,...
        BB_i_start,BB_i_mid,BB_i_end,ClipThresholdFmax,ClipFthresholdFnormal,1);
    if 0
        figure(996);
        plot(FclipIB,FclippedUpper,'.');
        hold on;
        plot(FclipEB,FclippedLower,'.');
    end
else
    plotclippingfigure=0;
    [FclipIB,FclipEB,FclipI,FclipE,~] = ClipDetection(FlowOrig,...
        BB_i_start,BB_i_mid,BB_i_end,ClipThresholdFmax,ClipFthresholdFnormal,settings.plotfigure & plotclippingfigure);
end

if sum([FclipI FclipE]*100)>0.5 %displaying clipping if more than 1% total Insp + Exp
    detectedClippingInOriginalData = [FclipI FclipE]*100;
    disp([num2str(n) ':' num2str(winNum) ') ' 'DetectedClipping: ' num2str(FclipI) '%, ' num2str(FclipE) '%']);
end

VIOrig = VI; %is the uncorrected VI and if there is cliping correction then may differ from VI per here:
if isfield(settings,'ApplyClippingCorrection') && settings.ApplyClippingCorrection==1
    if 0
        load('ClippingCorrectionF')
    else
        VTiCorrectionF = @(x)(((0.02063*tan(x*1.5))+0.5^0.1).^10+0.5);
        VTeCorrectionF = @(x)(((0.0527*tan(x*1.2))+0.5^0.1).^10+0.5);
        
    end
    if isfield(settings,'intentionalclipping')==1 && settings.intentionalclipping==1
        
        VTi_corrected = VTi_clipped.*VTiCorrectionF(FclipIB(:)');
        VTe_corrected = VTe_clipped.*VTeCorrectionF(FclipEB(:)');
        temp = [VTi_clipped(:) VTi_corrected(:) VTi(:)  FclipIB]; %for evaluation
        
    else
        %correct the data for clipping
        VTiOrig = VTi;
        VTeOrig = VTe;
        VTOrig = VT;
        VTi = VTi.*VTiCorrectionF(FclipIB(:)');
        VTe = VTe.*VTeCorrectionF(FclipEB(:)');
        VT = (VTi.*Te+VTe.*Ti)./(Ttot_B);
        VT(VT<0)=0; %just in case
        VE = VT./Ttot_B;
        VI = VT./Ttot_B;
        temp = [VI(:) VIOrig(:) FclipIB FclipEB]; %for evaluation
    end
end


if settings.verbose==2
   delta_t = etime(clock, clocktimestart); % delta in seconds
    str = ['Analysis duration (PostClipDetect): ', num2str(round(delta_t,1)), ' s']; disp(str);
end

%% Adjust BB times to be nearer zero crossings
BB_original = [BB_i_start BB_i_mid BB_i_end]; % make a backup of original times
BB_Times = BB_original; %default
try
if settings.nearerzerocrossings
    Vflow_out_temp = Vflow_out;
    Vflow_out_temp(isnan(Vflow_out_temp))=0; %Vflow_out as used below can not contain NaNs
    try
        BBs = length(BB_i_start); BB_i_max = NaN(BBs,1); BB_i_min = NaN(BBs,1);
        for bb = 1:BBs
            % find indx of max in insp
            [~, BB_i] = max(Vflow_out_temp(BB_i_start(bb):BB_i_mid(bb)));
            BB_i_max(bb) = BB_i + BB_i_start(bb);
            % find indx min in exp
            [~, BB_i] = min(Vflow_out_temp(BB_i_mid(bb):BB_i_end(bb)));
            BB_i_min(bb) = BB_i + BB_i_mid(bb);
        end
        % find zero crosses
        zero_crosses=(diff(sign(Vflow_out_temp)) & Vflow_out_temp(2:end)) .* sign(Vflow_out_temp(2:end));
        % find indx of upwards and downwards zero crossings
        up_x = 1+find(zero_crosses>0); down_x = 1+find(zero_crosses<0);
        % for each up_x, check if the value before it is closer to zero
        for u=1:length(up_x)
            t = abs(Vflow_out_temp(up_x(u)-1));
            if t<abs(Vflow_out_temp(up_x(u))); up_x(u) = up_x(u)-1; end
        end
        % for each down_x, check if the value before it is closer to zero
        for u=1:length(down_x)
            t = abs(Vflow_out_temp(down_x(u)-1));
            if t<abs(Vflow_out_temp(down_x(u))); down_x(u) = down_x(u)-1; end
        end
        % reset BB_i_mid
        for bb = 1:BBs
            if 1
                if Apnea_B(bb)
                    continue
                end
            end
            % New mid point is the dn_x that sits at the intersect of:
            % dn_x after insp max (BB(1:end)), and dn_x before exp min (BB(1:end)),
            % if multiple found, it is the closest in time to the original BB_i_mid
            BB_i_mid_candidate = intersect(find(down_x>BB_i_max(bb)), find(down_x<BB_i_min(bb)));
            if ~isempty(BB_i_mid_candidate)
                if length(BB_i_mid_candidate)==1
                    BB_i_mid(bb) = down_x(BB_i_mid_candidate);
                else
                    [~,indx]= min(abs(BB_i_mid(bb) - down_x(BB_i_mid_candidate)));
                    BB_i_mid(bb) = down_x(BB_i_mid_candidate(indx));
                end
                % else, we didn't find a single candidate, just use the original time
            end
        end
        % reset BB_i_start for breaths 2 to end, and BB_i_end for breaths 1 to end-1.
        for bb = 2:BBs
            if 1
                if Apnea_B(bb)
                    continue
                end
            end
            % New start/end point is the up_x that sits at the intersect of:
            % up_x before insp max (BB(2:end)), and up_x after exp min (BB(1:end-1)).
            % if multiple found, it is the closest in time to the original values
            BB_i_candidate = intersect(find(up_x<BB_i_max(bb)), find(up_x>BB_i_min(bb-1)));
            if ~isempty(BB_i_candidate)
                if length(BB_i_candidate)==1
                    BB_i_start(bb) = up_x(BB_i_candidate);
                    BB_i_end(bb-1) = up_x(BB_i_candidate);
                else
                    [~,indx]= min(abs(BB_i_start(bb) - up_x(BB_i_candidate)));
                    BB_i_start(bb) = up_x(BB_i_candidate(indx));
                    BB_i_end(bb-1) = up_x(BB_i_candidate(indx));
                end
                % else, we didn't find a single candidate, just use the original time
            end
        end
    catch me
        % if this process fails, then reset back to the original timing
        BB_i_start = BB_original(:,1); BB_i_mid = BB_original(:,2); BB_i_end = BB_original(:,3);
        disp(me.message); % me.getReport
    end
    BB_Times = [BB_i_start BB_i_mid BB_i_end];
end

BB_delta = BB_Times - BB_original; %unused but keeps track of the change in breath timing with this function
catch me
    disp(['failed to run nearerzerocrossings']);
end
Iflow.starti = BB_i_start; Iflow.midi = BB_i_mid; Iflow.endi = BB_i_end;    

%% Check if inverted trace
VTeprev = VTe(1:end-1);
VTicurr = VTi(2:end);
VTecurr = VTe(2:end);
temp = Apnea_B(2:end);
VTeprev(temp==1)=[];
VTicurr(temp==1)=[];
VTecurr(temp==1)=[];
diffcurr = median(abs((VTicurr-VTecurr)))/mean(VT);
diffprev = median(abs((VTicurr-VTeprev)))/mean(VT);
InvertedFlow = zeros(size(VTe));
if diffcurr>diffprev
    InvertedFlow = InvertedFlow + 1; 
    if settings.verbose
        disp(['Warning: Flow signal appears upside down, F=' num2str(-100*(diffprev-diffcurr),2)]);
    end
end


if settings.verbose==2
   delta_t = etime(clock, clocktimestart); % delta in seconds
    str = ['Analysis duration (CheckInvt): ', num2str(round(delta_t,1)), ' s']; disp(str);
end

%% Flow shape section
BreathFLDataTable=NaN;
if settings.flowshapesonly~=-1&&settings.flowshapesonly~=2
    
    %% Ttrans timing
    % Use the existing breath timing, and determine modified timing using the
    % Ttran threshold method
    
    range = 0.95; % 97 looks better (completely unvalidated, observation only)
    [BB_Ttrans,TiTrans,TeTrans] = TTransFromFlow(time, BB_Times, Vflow_out, range, Apnea_B, dt);
    if sum(isnan(BB_Ttrans)) > 0
        if settings.verbose; disp('NaN found in TTran timing'); end
    end
    
    %% Plot Flows and Breath Timings
    if 0&&settings.plotfigure
        figure(3); clf(figure(3));
        plot(time, Flow, 'g'); hold on;
        plot(time, Vflow, 'r');
        plot(time, Vflow_out, 'k');
        refline(0,0);
        flow_for_plot = Vflow_out;
        % transition time
        plot(Time(BB_Ttrans(:,1)), flow_for_plot(BB_Ttrans(:,1)), 'm^');
        plot(Time(BB_Ttrans(:,2)), flow_for_plot(BB_Ttrans(:,2)), 'm^');
        plot(Time(BB_Ttrans(:,3)), flow_for_plot(BB_Ttrans(:,3)), 'mv');
        plot(Time(BB_Ttrans(:,4)), flow_for_plot(BB_Ttrans(:,4)), 'mv');
        % original times
        plot(Time(BB_original(:,1)), flow_for_plot(BB_original(:,1)), 'ro');
        plot(Time(BB_original(:,2)), flow_for_plot(BB_original(:,2)), 'rx');
        plot(Time(BB_original(:,3)), flow_for_plot(BB_original(:,3)), 'r*');
        % good times
        plot(Time(BB_Times(:,1)), flow_for_plot(BB_Times(:,1)), 'bo');
        plot(Time(BB_Times(:,2)), flow_for_plot(BB_Times(:,2)), 'bx');
        plot(Time(BB_Times(:,3)), flow_for_plot(BB_Times(:,3)), 'b*');
    end
    
    %% Flow shape analysis
    if settings.verbose
        t_FLanalysis_start = clock;
        %disp(['Running FlowShape analysis']);%: settings.flowshapesonly=' num2str(settings.flowshapesonly)]);
    end
    
    try
        % BreathFLData is a table, with breaths per row, and each col is a feature,
        % the names of the individual features are provided in FtrName.       
        
        if 1 %one large table with multiple options
            [BreathFLDataTable1] = ComputeBreathFeatures(time, Vflow_out, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[0 1 0]); %[0: downsampleHz; 1=original timing, no resample]
            [BreathFLDataTable2] = ComputeBreathFeatures(time, Vflow_out, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[0 0 0]); %[0: downsampleHz, 0=DM timing, no resample]
            %[BreathFLDataTable3] = ComputeBreathFeatures(time, Vflow_out, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[0 1 250]); %[0: downsampleHz, 1=original timing, 250 = resample]
            %[BreathFLDataTable4] = ComputeBreathFeatures(time, Vflow_out, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[0 1 100]); %[0: downsampleHz, 1=original timing, 100 = resample]
            
            for i=1:size(BreathFLDataTable2,2)
                BreathFLDataTable1.Properties.VariableNames{i}=[BreathFLDataTable1.Properties.VariableNames{i} '_O'];
                BreathFLDataTable2.Properties.VariableNames{i}=[BreathFLDataTable2.Properties.VariableNames{i} '_T'];
                %BreathFLDataTable3.Properties.VariableNames{i}=[BreathFLDataTable3.Properties.VariableNames{i} '_O_250'];
                %BreathFLDataTable4.Properties.VariableNames{i}=[BreathFLDataTable4.Properties.VariableNames{i} '_O_100'];
            end
            %BreathFLDataTable = [BreathFLDataTable1 BreathFLDataTable2 BreathFLDataTable3];% BreathFLDataTable4];
            BreathFLDataTable = [BreathFLDataTable1 BreathFLDataTable2];
        else
            BreathFLDataTable = ComputeBreathFeatures(time, Vflow2, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[25 1]); %[25: downsampleHz, 1=original timing]
        end
        if settings.verbose
            delta_t = etime(clock, t_FLanalysis_start); % delta in seconds
            D_FL = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
            str = ['FlowShape analysis complete: ', char(D_FL), ' (hh:mm:ss)']; disp(str);
        end
    catch me
        disp(me.message); me.getReport
    end
end

%% Snore analysis
BreathSnoreTable=NaN;
if 0% settings.AnalyzeSnore==1 % I moved this out of this script
    if settings.verbose
        t_FLanalysis_start = clock;
        %disp(['Running FlowShape analysis']);%: settings.flowshapesonly=' num2str(settings.flowshapesonly)]);
    end
    
    try 
        % window-level analysis not compatible with SnoreInterpStruct
        load([settings.workdir,'Converted\',settings.filename(1:end-4),'_snore.mat'],...
            'SnoreChannelsList','SnoreInterpStruct','SnoreStruct')
        [BreathSnoreTable] = ComputeSnoreFeaturesOld(Time,SigT,SnoreInterpStruct,SnoreStruct,BB_Times,0);
        [BreathSnoreTable_,~,~] = ComputeSnoreFeatures(SigT,...
                SnoreInterpStruct,SnoreStruct,Time(BB_Times),0);
            
            % SigT replaces DataEventHypnog_Mat to run in AnalysisOne.m-RMA-5/3/2021
            % please make modifications in any main/sub functions as needed
        
        
        if settings.verbose
            delta_t = etime(clock, t_FLanalysis_start); % delta in seconds
            D_FL = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
            str = ['Snore analysis complete: ', char(D_FL), ' (hh:mm:ss)']; disp(str);
        end
    catch me
        disp(me.message); me.getReport
    end
end

%% Compute FL breaths and plot

if settings.flowshapesonly ~= -1
    try
    
    try
        [FlowDrive,~] = computeFlowDrive(FlowDriveModels.FinalModelTable,BreathFLDataTable);
    catch me
        disp('failed to compute flowdrive');
        FlowDrive=NaN*VI(:);
    end
    FlowDriveOut = table(FlowDrive); %aka "FlowDrive25"
    
    if isfield(settings,'computeFlowDrive10') && settings.computeFlowDrive10==1
        % if we do this successfully, only keep this 10 hz output as the
        % returned values. the 25 hz values made above are used for plots
        try
            [FlowDrive10,~] = computeFlowDrive(FlowDriveModels.FinalModelTable10,BreathFLDataTable);
            FlowDrive=FlowDrive10; %overwrite for plot
            FlowDriveOut.FlowDrive10=FlowDrive10;
        catch me
            disp('failed to compute flowdrive10');
            FlowDrive=NaN*VI(:);
        end
    end
    
    %DM to add certainty info to table 
    try
        [FlowLimCertOut,pihat] = computeFLcertainty(FlowDriveModels.FLcert,BreathFLDataTable);
        %FlowDriveOut.etc = etc
        %FlowDriveOut = [FlowDriveOut FlowLimCertOut];
    end
    
    if 1&&settings.plotfigure
        try
            
            
            Nsubplots=2;
            figure(FDfig); clf(FDfig); FDfigID = gcf;
            FDfigID.Color = [1 1 1];
            FDfigID.Units = 'Inches';
            FDfigID.Position = [2 1 12 9];
            
            ax99(1)=subplot(Nsubplots,1,1);

            plotbreathgrid=1;
            flowplotpos=nanmean(Vflow_out);
            gridheight=2*nanstd(Vflow_out);
            if plotbreathgrid
                xgridlist = Time(BB_i_start)';
                plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^2);
                hold('on');
                xgridlist = Time(BB_i_mid)';
                plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^4);
                %zero flow baseline:
                plot(Time,0*Time,'color',0.8*[1 0.5 0.5]);
            end
            
            plot(Time,Vflow_out,'k'); hold on
            
            if 0 % ToDo: Sort out why vflow_out shifts start times poorly
                %dlm debugging
                plot(Time,Vflow_out*10,'k'); hold on
                plot(Time,Vflow*10,'r');
                plot(Time, Th,'g');
                %plot(Time, Ab,'b');
                RIPflow = diff(Th); RIPflow=smooth(RIPflow,30);RIPflow=[RIPflow;NaN];
                plot(Time, -1*(RIPflow*20), 'm');
            end
            
            set(gca,'xtick',[],'xcolor',[1 1 1],'box','off','tickdir','out');
            ylabel('Flow');
            
            % DLM - if we did both 10 and 25 hz flowdrive, then show the 25hz
            % data in plot as well. plot 25hz first so it's underneath
            if isfield(settings,'computeFlowDrive10') && settings.computeFlowDrive10==1
                FlowDriveOutPlot = FlowDriveOut.FlowDrive*100; %plot original underneath
                FlowDriveOutPlot(FlowDriveOutPlot>105)=105;
                FlowDriveOutPlot(FlowDriveOutPlot<0)=0;
                FlowDriveOutPlot2 = FlowDriveOutPlot;
                FlowDriveOutPlot2(Apnea_B==1)=0;
                ax99(2)=subplot(Nsubplots,1,2);
                stairs(Time(BB_i_start),FlowDriveOutPlot2,'--', 'Color', [0.7 0 0]); hold('on'); %plot earlier 25Hz version
                stairs(Time(BB_i_start),FlowDriveOutPlot, 'Color', [0.7 0 0]);
            end
            
            FlowDriveOutPlot = FlowDrive*100;
            FlowDriveOutPlot(FlowDriveOutPlot>105)=105;
            FlowDriveOutPlot(FlowDriveOutPlot<0)=0;
            FlowDriveOutPlot2 = FlowDriveOutPlot;
            FlowDriveOutPlot2(Apnea_B==1)=0;
            ax99(2)=subplot(Nsubplots,1,2);
            stairs(Time(BB_i_start),FlowDriveOutPlot2,'--', 'Color', [0 0 1]); hold('on');
            stairs(Time(BB_i_start),FlowDriveOutPlot, 'Color', [0 0 1]);
            plot(Time,100+0*Time,'k--');
            plot(Time,50+0*Time,'k--');
            plot(Time,0+0*Time,'k--');
            linkaxes(ax99,'x');
            ylim([-5 106]);
            set(gca,'box','off','tickdir','out');
            ylabel('Flow:Drive, %');
            
            if Nsubplots>2
                set(gca,'xtick',[],'xcolor',[1 1 1]);
                
                EFLindex=BreathFLDataTable.AsymIndex_O;
                EFLindex(EFLindex<0.5)=0.5;
                ax99(3)=subplot(3,1,3);
                stairs(Time(BB_i_start),EFLindex); hold('on');
                plot(Time,Time*0+0.8);
                linkaxes(ax99,'x');
                ylim([0.5 1]);
                set(gca,'box','off','tickdir','out');
                ylabel('Exp. Prolapse');
            end
        catch me
            'failed to plot FlowDrive Figure'
        end
    end   
    catch me
       disp('failing FlowDrive, skipped')
    end
end

%% save flowdrive plot
if (settings.saveplots) && (settings.flowshapesonly ~= -1)
    % save in folder for each patient
    saveloc=[settings.OutputDataDirectory settings.savename filesep settings.filename];
    if ~(exist(saveloc, 'dir') == 7)
        mkdir(saveloc);
    end
    % save FDfigID, in the file format specified.
    switch settings.savefigas
        case 'saveasTIFF'
            print(FDfigID, [saveloc, '\FD_window=',num2str(winNum)], '-dtiff', '-r300');
        case 'saveasPNG'
            saveas(FDfigID, [saveloc, '\FD_window=',num2str(winNum)], 'png');
        case 'saveasFIG'
            savefig(FDfigID, [saveloc, '\FD_window=',num2str(winNum)],'compact');
        % 2021-11-02 12:15:33 Tuesday EAS added this option to save both .fig and .png plots
        case 'saveasFIGandPNG'
            FDfigID.WindowState = 'maximized';
            saveloc1 = [saveloc '\fdfigs'];
            saveloc2 = [saveloc '\fdpngs'];
            if ~(exist(saveloc1, 'dir') == 7)
                mkdir([saveloc1]);
            end                    
            if ~(exist(saveloc2, 'dir') == 7)
                mkdir([saveloc2]);
            end                  
            savefig(FDfigID, [saveloc1, '\FD_window=',num2str(winNum)],'compact');
            % saveas(FDfigID, [saveloc, '\FD_window=',num2str(winNum)], 'png');
            % exportgraphics is better than saveas because there's no
            % blank space around the figure
            exportgraphics(FDfigID, [saveloc2, '\FD_window=',num2str(winNum),'.png']);
    end
end


%%
if 0
    % harsh cleaner
    Table1 = CleanAwayNaNs(BreathFLDataTable1);
    Table3 = CleanAwayNaNs(BreathFLDataTable3);
    mat1 = Table1{:,:};     % standard duration
    mat2 = Table3{:,:};     % double duration
    FeatureNames = BreathFLDataTable1.Properties.VariableNames;
    
    % correlations
    corrVals = NaN(90,1);
    for ftr=1:90
        corrVals(ftr) = corr(mat1(:,ftr), mat2(:,ftr));
    end
    oddities = find(corrVals<=0.95);
    oddlist = [num2cell((oddities)) , num2cell(corrVals(oddities)), FeatureNames(oddities)'];
    oddTable = cell2table(oddlist, 'VariableNames', {'Ftr', 'Correlation', 'FtrName'});  oddTable
    
    % absolute deltas, and percent change
    mat_diff = (mat1 - mat2) ./ abs(mat1);
    percentchange = abs(mean(mat_diff,1));
    odds2 = find(percentchange>0.1);
    
    oddlist2 = [num2cell((odds2))', num2cell(corrVals(odds2)), num2cell(percentchange(odds2))', FeatureNames(odds2)'];
    oddTable2 = cell2table(oddlist2, 'VariableNames', {'Ftr', 'Correlation', 'Change', 'FtrName'}); oddTable2
    
    % scatter plot
    if 0
        ftr = 5;
        figure(1); clf(figure(1));
        scatter(mat1(:,ftr), mat2(:,ftr), 3,'filled','markerfacealpha',0.8);
        xlabel('normal duration'); ylabel('double duration');
    end
    
    %     figure(999); clf(figure(999));
    %     subplot(1,2,1);
    %     NaNSpace = isnan(BreathFLDataTable{:,:}); NaNSpaceMap = NaNSpace*90;
    %     image(NaNSpaceMap); xlabel('Table columns'); ylabel('breath number');
    %     title('All breaths');
    %     subplot(1,2,2);
    %     BreathFLDataTable_NoApnea = BreathFLDataTable{:,:};
    %     BreathFLDataTable_NoApnea(Apnea_B==1,:) = [];
    %     NaNSpace = isnan(BreathFLDataTable_NoApnea); NaNSpaceMap = NaNSpace*90;
    %     image(NaNSpaceMap); xlabel('Table columns'); ylabel('breath number');
    %     title('Excluding apnoea');
end

%% Oral and Nasal separate Flow Analysis
if any(SigT.Properties.VariableNames == "FlowOral")
    try
        [VToral,VTioral,VTeoral,VIoral,VEoral] = VEfromFlowWithoutBreathDetection(SigT.FlowOral, BB_i_start,BB_i_mid,BB_i_end,Ti,Te,Apnea_B,dt);
        
        FlowNasal = SigT.Flow - SigT.FlowOral; % same treatment for oral flow as total flow 
        [VTnasal,VTinasal,VTenasal,VInasal,VEnasal] = VEfromFlowWithoutBreathDetection(FlowNasal, BB_i_start,BB_i_mid,BB_i_end,Ti,Te,Apnea_B,dt);

        FVTinasal = VTinasal./(VTinasal+VTioral);
        FVTinasal(Apnea_B(:)==1)=NaN;
        FVTenasal = VTenasal./(VTenasal+VTeoral);
        FVTenasal(Apnea_B(:)==1)=NaN;
        
        FlowCompartments = table(VToral,VTioral,VTeoral,VIoral,VEoral,...
            VTnasal,VTinasal,VTenasal,VInasal,VEnasal,FVTinasal,FVTenasal);
    end
end

%% Stim ON
if any(SigT.Properties.VariableNames == "StimON")
    StimON = SigT.StimON*1;
    StimONinterp = interp1(1:length(StimON),StimON,BB_i_mid,'previous');
    StimONtbl = array2table(StimONinterp);
    StimONtbl.Properties.VariableNames = {'StimON'};
end

if any(SigT.Properties.VariableNames == "HGNSon")
    HGNSon = SigT.HGNSon*1;
    HGNSoninterp = interp1(1:length(HGNSon),HGNSon,BB_i_mid,'previous');
    HGNSonTbl = array2table(HGNSoninterp);
    HGNSonTbl.Properties.VariableNames = {'HGNSon'};
end

if any(SigT.Properties.VariableNames == "StimAmpl")
    StimAmpl = SigT.StimAmpl*1;
    StimAmplinterp = interp1(1:length(StimAmpl),StimAmpl,BB_i_mid,'previous');
    StimAmplTbl = array2table(StimAmplinterp);
    StimAmplTbl.Properties.VariableNames = {'StimAmpl'};
end

if any(SigT.Properties.VariableNames == "StimSelect")
    StimSelect = SigT.StimSelect*1;
    StimSelectinterp = interp1(1:length(StimSelect),StimSelect,BB_i_mid,'previous');
    StimSelectTbl = array2table(StimSelectinterp);
    StimSelectTbl.Properties.VariableNames = {'StimSelect'};
end

%% Pes analysis (if it exists)
if sum(strcmp(ChannelsList,'Pes')==1)
    settings.PmusEcw=30; %changed from 10 to 30, 1/25/21
    Pes = SigT.Pes; % not indexed, 'find' stays
    %Filtered Pes (gentle low pass)
    [B_butter0,A_butter0] = butter(2,5/(1/dt/2),'low');
    Pes = filtfilt(B_butter0,A_butter0,Pes);
    
    if sum(isnan(Pes))==0
        [DeltaPes,DeltaPmus,Ipes,VIpes,Pesbaselineest]=PesAnalysis(Pes,Vflow_out,Time,Iflow,settings.PmusEcw,1); %pause;
        DeltaPes=AlignWithBreaths(DeltaPes,Time,Iflow,Ipes);
        BaselinePes=AlignWithBreaths(Pes(Ipes.starti),Time,Iflow,Ipes); %for plot only
        DeltaPmus=AlignWithBreaths(DeltaPmus,Time,Iflow,Ipes);
        VIpes=AlignWithBreaths(VIpes,Time,Iflow,Ipes);
        try
            % what figure number should this be?
            if settings.plotfigure&&0
                figure(211); clf(211);
                set(gcf,'color',[1 1 1]);
                hold('on')
                ylabel('Pes');
                plot(Time,0*Pes,'k');
                plot(Time,Pes-Pesbaselineest,'b');
                stairs(Time(Ipes.starti),Pes(Ipes.starti)-Pesbaselineest,'r');
                stairs(Time(Iflow.starti),BaselinePes-DeltaPes-Pesbaselineest,'r'); %not exactly the right baseline
                hold('off');
            end
        catch me
        end
    else
        DeltaPes = NaN*VI';
        DeltaPmus = NaN*VI';
        VIpes = NaN*VI';
    end
    
else
    DeltaPes = NaN*VI';
    DeltaPmus = NaN*VI';
    VIpes = NaN*VI';
end

%% Edi analysis (if it exists)
DeltaEdi = NaN*VI';
PeakEdi = NaN*VI';
VIedi = NaN*VI';
if sum(strcmp(ChannelsList,'Edi')==1)
    try
        Edi = SigT.Edi;
        
        %check for runs of exact zeros in Edi, after, window completely
        I = SigT.Edi==0;
        if sum(I)>0 %if any zeros found
            disp('warning: found exact zeros in Edi')
            I1=find(diff(I)>0);
            I2=find(diff(I)<0);
            [I1,I2]=TidyStartEndEventList(I1,I2,length(Edi));
            ErrorLength = I2-I1;
            for i=1:length(I1)
                if ErrorLength(i)>3 %more than 3 sec of exact zeros are removed.
                    disp('warning: found more than 3 s of exact zeros in Edi')
                    Edi(I1(i):I2(i))=NaN;
                end
            end
        end
        
        if sum(isnan(Edi))==0
            [DeltaEdi,~,Iedi,VIedi,Edibaselineest,PeakEdi]=PesAnalysis(-Edi,Vflow_out,Time,Iflow,0,0);
            [DeltaEdi]=AlignWithBreaths(DeltaEdi,Time,Iflow,Iedi);
            [PeakEdi]=AlignWithBreaths(PeakEdi,Time,Iflow,Iedi);
            BaselineEdi=AlignWithBreaths(-Edi(Iedi.starti),Time,Iflow,Iedi); %for plot only
            VIedi=AlignWithBreaths(VIedi,Time,Iflow,Iedi);
            try
                % what figure number should this be?
                if settings.plotfigure
                    figure(22); clf(22);
                    hold('on')
                    ylabel('-Edi');
                    plot(Time,-Edi-Edibaselineest,'b');
                    stairs(Time(Iedi.starti),-Edi(Iedi.starti)-Edibaselineest,'r');
                    stairs(Time(Iflow.starti),BaselineEdi-DeltaEdi-Edibaselineest,'r'); %not exactly the right baseline
                    hold('off');
                end
            catch me
            end
        end
    catch me  
        disp('warning, failed Edi analysis; results will be NaN');
    end
end

%% FlowEdi analysis (if it exists)
if sum(strcmp(ChannelsList,'FlowEdi')==1)
    FlowEdi = SigT.FlowEdi;
    if sum(isnan(Edi))==0
        [~,~,BB_i_startE,BB_i_midE,BB_i_endE,~,FlowEdi_VI,~,~,~,~,~,~,~,~,~,~,~] = ...
            VEfromFlow_sqrt_V16(Time,FlowEdi);
        Iflowedi.starti = BB_i_startE;
        Iflowedi.midi = BB_i_midE;
        Iflowedi.endi = BB_i_endE;
        FlowEdi_VI=AlignWithBreaths(FlowEdi_VI,Time,Iflow,Iflowedi);
    else
        FlowEdi_VI=NaN*VI';
    end
else
    FlowEdi_VI=NaN*VI';
end

%% FlowPes analysis (if it exists)
if sum(strcmp(ChannelsList,'FlowPes')==1)
    FlowPes = SigT.FlowPes;
    if sum(isnan(Pes))==0
        [~,~,BB_i_startP,BB_i_midP,BB_i_endP,~,FlowPes_VI,~,~,~,~,~,~,~,~,~,~,~] = ...
            VEfromFlow_sqrt_V16(Time,FlowPes);
        Iflowpes.starti = BB_i_startP;
        Iflowpes.midi = BB_i_midP;
        Iflowpes.endi = BB_i_endP;
        FlowPes_VI=AlignWithBreaths(FlowPes_VI,Time,Iflow,Iflowpes);
    else
        FlowPes_VI=NaN*VI';
    end
    FlowPes_VI(isnan(DeltaPes))=NaN;
else
    FlowPes_VI=NaN*VI';
end

%% GGpmax analysis (if it exists)
GGtonic = NaN*VI'; GGpeak = NaN*VI';
if sum(strcmp(ChannelsList,'GGpmax')==1)
    GGpmax = SigT.GGpmax;
    for i=1:length(Iflow.starti)
        GGtonic(i) = min(GGpmax(Iflow.starti(i):Iflow.endi(i)));
        GGpeak(i) = max(GGpmax(Iflow.starti(i):Iflow.endi(i)));
        if sum(isnan(GGpmax(Iflow.starti(i):Iflow.endi(i))))>0 %data output as unknown if any missing signal within a particular breath
            GGtonic(i)=NaN;
            GGpeak(i)=NaN;
        end
    end
end
%sum(GGpmax)


if settings.verbose==2
    delta_t = etime(clock, clocktimestart); % delta in seconds
    str = ['Analysis duration (PostPes): ', num2str(round(delta_t,1)), ' s']; disp(str);
end

%% FOT analysis (if it exists)
try
if sum(strcmp(ChannelsList,'YFOT')==1)
    
    FotT=table();
    
    figure(55);
    ax55(1)=subplot(6,1,1); plot(Time,SigT.YFOTclean);
    ax55(2)=subplot(6,1,2); plot(Time,SigT.YFOT);
    ax55(3)=subplot(6,1,3); plot(Time,SigT.CohFOT);
    ax55(4)=subplot(6,1,4); plot(Time,SigT.PflowFOT);
    ax55(5)=subplot(6,1,5); plot(Time,SigT.PpmaskFOT);
    ax55(6)=subplot(6,1,6); plot(Time,SigT.Flow);
    linkaxes(ax55,'x');
    
    blankcol = nan(length(BB_i_start),1);
    

    FotT.CohFOTmean=blankcol;
    FotT.CohFOTmin=blankcol;
    FotT.PpmaskFOTmean=blankcol;  
    FotT.FOTFnotisnan=blankcol;
    FotT.YFOTmean=blankcol;
    FotT.YFOTmeanI=blankcol;
    FotT.YFOTmeanIW=blankcol;
    FotT.YFOTminI=blankcol;
    FotT.YFOTmeanE=blankcol;
    FotT.YFOTmeanEW=blankcol;
    FotT.YFOTminE=blankcol;    
    FotT.YFOTendE=blankcol;
    FotT.YFOTmaxE=blankcol;
    
for i=1:length(BB_i_start)
    I = BB_i_start(i):BB_i_end(i);
    FotT.YFOTmean(i) = nanmean(SigT.YFOT(I));
    
    FotT.CohFOTmean(i) = nanmean(SigT.CohFOT(I));
    FotT.CohFOTmin(i) = min(SigT.CohFOT(I));
    FotT.PpmaskFOTmean(i) = nanmean(SigT.PpmaskFOT(I));
    FotT.PflowFOTmean(i) = nanmean(SigT.PflowFOT(I));
    FotT.FOTFnotisnan(i) = nanmean(~isnan(SigT.YFOT(I)));
    
    I = BB_i_start(i):BB_i_mid(i);
    FotT.YFOTmeanI(i) = nanmean(SigT.YFOT(I));
    FotT.YFOTmeanIW(i) = nanmean(SigT.YFOT(I).*SigT.Flow(I))./nanmean(SigT.Flow(I)+0*SigT.YFOT(I));
    FotT.YFOTminI(i) = min(SigT.YFOT(I));
    
    I = BB_i_mid(i):BB_i_end(i);
    FotT.YFOTmeanE(i) = nanmean(SigT.YFOT(I));
    FotT.YFOTmeanEW(i) = nanmean(SigT.YFOT(I).*SigT.Flow(I))./nanmean(SigT.Flow(I)+0*SigT.YFOT(I));
    FotT.YFOTminE(i) = min(SigT.YFOT(I));
    
    FotT.YFOTendE(i) = SigT.YFOT(I(end));
end
    FotT{FotT.PpmaskFOTmean<0.25,5:end}=NaN;

end

%% Pepi analysis (if it exists)
try
if sum(strcmp(ChannelsList,'Pepi')==1)

    Pepi = SigT.Pepi; % not indexed, 'find' stays
    %Filtered Pes (gentle low pass)
    [B_butter0,A_butter0] = butter(2,5/(1/dt/2),'low');
    Pepi = filtfilt(B_butter0,A_butter0,Pepi);
    
    dtW = 0.25;
    TimeW = [dtW:dtW:settings.windowlength-dtW]*60;
    clear Y Yt;
    for i=1:length(TimeW)
       I = Apnea_B==0 & (Time(BB_i_start)-Time(1)) > TimeW(i)-1*dtW*60 & (Time(BB_i_start)-Time(1)) < TimeW(i)+1*dtW*60;
       Y(i,1) = median([Pepi(BB_i_start(I));Pepi(BB_i_mid(I))]);
       Yt(i,1) = median([Time(BB_i_start(I));Time(BB_i_mid(I))]); 
    end
    PepiBaseline1 = median([Pepi(BB_i_start(Apnea_B==0));Pepi(BB_i_mid(Apnea_B==0))]);
    PepiBaseline = interp1(Yt,Y,Time,'linear','extrap');
        PepiBaseline(Time<Yt(1))=Y(1);
        PepiBaseline(Time>Yt(end))=Y(end);
        
    PepiBaseline0 = mean(PepiBaseline);    
    PepiBaseline2 = prctile(Pepi,100*median(Ti./Ttot_B));
    % what figure number should this be?
            if 1%settings.plotfigure&&0
                figure(218); clf(218);
                set(gcf,'color',[1 1 1]);
                
                ylabel('Pepi');
                ax218(1)=subplot(3,1,1);
                plot(Time,0*Pepi,'k:');
                hold('on')
                plot(Time,0*Pepi + PepiBaseline - PepiBaseline0,'k');
                
                plot(Time,Pepi-PepiBaseline0,'b');
                plot(Time(BB_i_start),Pepi(BB_i_start)-PepiBaseline0,'r.');
                ax218(2)=subplot(3,1,2);
                plot(Time,Vflow_out,'k'); %Flow
                linkaxes(ax218,'x');
                ax218(3)=subplot(3,1,3);
                plot(Time,Ab,'k');
                linkaxes(ax218,'x');
            end
            
            DeltaPepi = nan(length(BB_i_start),1);
            DeltaPepit = nan(length(BB_i_start),1);
            Sig = Pepi - PepiBaseline;
            Sig2 = Vflow_out;
            DeltaPepiE = nan(length(BB_i_start),1);
            DeltaPepiEt = nan(length(BB_i_start),1);
            Gepi200 = nan(length(BB_i_start),1);
	    Gepi100 = nan(length(BB_i_start),1);
            GepiPeakPepiI = nan(length(BB_i_start),1);
            GepiPeakPepiE = nan(length(BB_i_start),1);
            
for i=1:length(BB_i_start)
    I = BB_i_start(i):BB_i_mid(i);
    [val,I2] = min(Sig(I));
    maxPepi(i,1) = max(Sig(I));
    DeltaPepi(i) = -val;
    DeltaPepit(i) = (I2-1)*dt;
    flowatDeltaPepi = Sig2([BB_i_start(i)+I2-1]);
    GepiPeakPepiI(i) = flowatDeltaPepi/DeltaPepi(i); %assumes flow is relative to zero flow; "Pepi - PepiBaseline" is adequately corrected. 
    
    I = BB_i_mid(i):BB_i_end(i);
    [val,I2] = max(Sig(I));
    DeltaPepiE(i) = val;
    DeltaPepiEt(i) = (I2-1)*dt;
    flowatDeltaPepiE = -Sig2([BB_i_mid(i)+I2-1]);
    GepiPeakPepiE(i) = flowatDeltaPepiE/DeltaPepiE(i); %assumes flow is relative to zero flow; "Pepi - PepiBaseline" is adequately corrected. 
    
    
    I = [BB_i_start(i) BB_i_start(i)+round(0.2/dt)]; %anchored at same start insp time
    Gepi200(i) = -diff(Sig2(I))/diff(Sig(I));
    I = [BB_i_start(i) BB_i_start(i)+round(0.1/dt)]; %anchored at same start insp time
    Gepi100(i) = -diff(Sig2(I))/diff(Sig(I));
end

stairs(ax218(1),Time(BB_i_start),-DeltaPepi + PepiBaseline(BB_i_start + round(DeltaPepit/dt)) - PepiBaseline0); % 

    figure(219)
    ax219(1)=subplot(2,1,1);
    [xcf,lags,bounds] = crosscorr(detrend(Sig),detrend(Sig2),'NumLags',round(nanmedian(Ttot_B)/dt));
    plot(lags*dt,xcf); 
    xlabel('Flow lags Pepi by ? s')
    ax219(2)=subplot(2,1,2);
    plot(lags/round(nanmedian(Ttot_B)/dt),xcf); 
    xlabel('Flow lags Pepi by ? cycles')
    
    [val,I] = min(xcf);
    PepiFlowCorr = -val + zeros(length(BB_i_start),1);
    PepiFlowCorrLag = lags(I)/round(nanmedian(Ttot_B)/dt) + zeros(length(BB_i_start),1);
    
    GepiPeakPepiI(Apnea_B==1)=0;
    GepiPeakPepiE(Apnea_B==1)=0;
    Gepi200(Apnea_B==1)=0;
    Gepi100(Apnea_B==1)=0;
    
    
%% Pepi artifact removal
PepiArtifact = 1
if PepiArtifact == 1
    
PepiUpperLim = 4*abs(median(DeltaPepiE+DeltaPepi));
temp2 = maxPepi>PepiUpperLim; % changed from 0.5*abs(minPesB2); %test for artifact: upwards swing from "baseline" is larger than the downwards swing from the "baseline"
temp3 = DeltaPepiE>PepiUpperLim
temp=temp2|temp3

DeltaPepi(temp) = NaN;
    DeltaPepit(temp) = NaN;
    GepiPeakPepiI(temp) = NaN;
    DeltaPepiE(temp) = NaN;
    DeltaPepiEt(temp) = NaN;
    GepiPeakPepiE(temp) = NaN;
    Gepi200(temp) = NaN;
    Gepi100(temp) = NaN;
    PepiFlowCorr(temp) = NaN;
    PepiFlowCorrLag(temp) = NaN;
end


%     figure(219)
%     fakelagi = 250;
%     temp = [zeros(fakelagi,1);Sig2(1:end-fakelagi)];
%     [xcf,lags,bounds] = crosscorr(detrend(Sig2),detrend(temp),'NumLags',round(nanmedian(Ttot_B)/dt));
%     plot(lags,xcf);

else
    DeltaPepi = nan(length(BB_i_start),1);
    DeltaPepit = nan(length(BB_i_start),1);
    GepiPeakPepiI = nan(length(BB_i_start),1);
    DeltaPepiE = nan(length(BB_i_start),1);
    DeltaPepiEt = nan(length(BB_i_start),1);
    GepiPeakPepiE = nan(length(BB_i_start),1);
    Gepi200 = nan(length(BB_i_start),1);
    Gepi100 = nan(length(BB_i_start),1);
    PepiFlowCorr = nan(length(BB_i_start),1);
    PepiFlowCorrLag = nan(length(BB_i_start),1);
end
catch
    DeltaPepi = nan(length(BB_i_start),1);
    DeltaPepit = nan(length(BB_i_start),1);
    GepiPeakPepiI = nan(length(BB_i_start),1);
    DeltaPepiE = nan(length(BB_i_start),1);
    DeltaPepiEt = nan(length(BB_i_start),1);
    GepiPeakPepiE = nan(length(BB_i_start),1);
    Gepi200 = nan(length(BB_i_start),1);
    Gepi100 = nan(length(BB_i_start),1);
    PepiFlowCorr = nan(length(BB_i_start),1);
    PepiFlowCorrLag = nan(length(BB_i_start),1);
end



PepiBreathDataTable = table(DeltaPepi,DeltaPepit,GepiPeakPepiI,DeltaPepiE,DeltaPepiEt,GepiPeakPepiE,Gepi200,Gepi100,PepiFlowCorr,PepiFlowCorrLag);



%% Convert events to breath domain
% clippedinsp_B=zeros(length(BB_i_start),1);
% clippedexp_B=zeros(length(BB_i_start),1);
AR=zeros(1,length(BB_i_start));
ARei=zeros(1,length(BB_i_start));
AReif=zeros(1,length(BB_i_start));
ARf=zeros(1,length(BB_i_start));
ARonset=zeros(1,length(BB_i_start));
% E=zeros(1,length(BB_i_start)); %Events
% E_ScoredCentralApnea=zeros(1,length(BB_i_start)); %EventsScoredCentralApnea
% E_ScoredCentralHyp=zeros(1,length(BB_i_start)); %EventsScoredCentralHyp
Te_i1 = round(median(BB_i_end-BB_i_start));
for k=1:length(BB_i_start)
    %clippedinsp_B(k)=mean(clipped(BB_i_start(k):BB_i_mid(k)));
    %clippedexp_B(k)=mean(clipped(BB_i_mid(k):BB_i_end(k)));
    if 1
        if k==1
            I=max([BB_i_start(k)-Te_i1 1]):BB_i_mid(k);
        else
            I=BB_i_mid(k-1):BB_i_mid(k);
        end
        I2=BB_i_start(k):BB_i_end(k);
        AR(k)=max(Arousal(I2)); %if arousal is within breath
        ARf(k)=mean(Arousal(I2)); %if arousal is within breath
        ARei(k)=max(Arousal(I)); %if arousal is within breath
        AReif(k)=mean(Arousal(I)); %if arousal is within breath
        %E(k)=1-(max(Events(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        %         if settings.eventsarebreathsfullywithinmargins
        %             E(k)=1-(min(Events(BB_i_start(k):BB_i_end(k)))); %if breath is entirely inside an event
        %             E_ScoredCentralApnea(k)=(min(EventsScoredCentralApnea(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        %             E_ScoredCentralHyp(k)=(min(EventsScoredCentralHyp(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        %         else
        %             E(k)=1-(max(Events(BB_i_start(k):BB_i_end(k)))); %if breath is even slightly inside an event
        %             E_ScoredCentralApnea(k)=(max(EventsScoredCentralApnea(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        %             E_ScoredCentralHyp(k)=(max(EventsScoredCentralHyp(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        %         end
    else
        AR(k)=Arousal(BB_i_start(k)); %if start of breath is in arousal
        %         E(k)=1-(round(mean(Events(BB_i_start(k):BB_i_end(k))))); %if most of breath is within event
    end
end

%% Try phasic/tonic REM

REMphasic=nan(1,length(BB_i_start));
REMtonic=nan(1,length(BB_i_start));
REMconjugate=nan(1,length(BB_i_start));
REMconjugateFThres=nan(1,length(BB_i_start));
try
if ~isempty(find(strcmp(ChannelsList,'PhasicREM')==1))
PhasicREM=SigT.PhasicREM; %EEG selection uses scoring

TonicREM=SigT.TonicREM; %EEG selection uses scoring
EOGconjugateFThres=SigT.EOGconjugateFThres; %EEG selection uses scoring
EOGconjugate=SigT.EOGconjugate; %EEG selection uses scoring
for k=1:length(BB_i_start)
    %clippedinsp_B(k)=mean(clipped(BB_i_start(k):BB_i_mid(k)));
    %clippedexp_B(k)=mean(clipped(BB_i_mid(k):BB_i_end(k)));

        I2=BB_i_start(k):BB_i_end(k);
        REMphasic(k)=mean(PhasicREM(I2)); %if arousal is within breath
        REMtonic(k)=mean(TonicREM(I2)); %if arousal is within breath
        REMconjugateFThres(k)=median(EOGconjugateFThres(I2)); %if arousal is within breath
        REMconjugate(k)=median(EOGconjugate(I2)); %if arousal is within breath
end
end
catch me
   disp('failed Phasic REM code');
   disp(me.message);
end

%% EEG Power breath level
PdeltalogfiltB=NaN*zeros(1,length(BB_i_start));
PthetalogfiltB=NaN*zeros(1,length(BB_i_start));
PalphalogfiltB=NaN*zeros(1,length(BB_i_start));
%PsigmalogfiltB=NaN*zeros(1,length(BB_i_start));
PbetalogfiltB=NaN*zeros(1,length(BB_i_start));
WakeSleepB=NaN*zeros(1,length(BB_i_start));
if sum(strcmp(ChannelsList,'Pdeltalogfilt')==1)
    Pdeltalogfilt = SigT.Pdeltalogfilt;
    Pthetalogfilt = SigT.Pthetalogfilt;
    Palphalogfilt = SigT.Palphalogfilt;
    %Psigmalogfilt = DataEventHypnog_Mat(:,(find(strcmp(ChannelsList,'Psigmalogfilt')==1)));
    Pbetalogfilt = SigT.Pbetalogfilt;
    if sum(strcmp(ChannelsList,'WakeSleep')==1)
        WakeSleep = SigT.WakeSleep;
    end
    for k=1:length(BB_i_start)
        if 1
            if k==1
                I=max([BB_i_start(k)-Te_i1 1]):BB_i_mid(k);
            else
                I=BB_i_mid(k-1):BB_i_mid(k);
            end
        else
            I=BB_i_start(k):BB_i_end(k);
        end
        PdeltalogfiltB(k)=median(Pdeltalogfilt(I)); %if arousal is within breath
        PthetalogfiltB(k)=median(Pthetalogfilt(I)); %if arousal is within breath
        PalphalogfiltB(k)=median(Palphalogfilt(I)); %if arousal is within breath
        %        PsigmalogfiltB(k)=median(Psigmalogfilt(I)); %if arousal is within breath
        PbetalogfiltB(k)=median(Pbetalogfilt(I)); %if arousal is within breath
        if exist('WakeSleep')
            WakeSleepB(k)=median(WakeSleep(I)); %if arousal is within breath
        end
    end
end

%% breath analysis of events, clinical scoring
Etype = 0*BB_i_start;
for k=1:length(BB_i_start)
    % 3=CentralA 2=ObstructiveA 4=OHypopnea 5=MixedA 6=CentralHypopnea
    temp = SigT{BB_i_start(k):BB_i_end(k),find(strcmp(ChannelsList,'EventsResp')==1)};
    F = sum(temp~=0)/length(temp);
    if settings.eventsarebreathsfullywithinmargins==1
        if F==1
            Etype(k)=mode(temp(temp>0));
        end
    else
        if F>0
            Etype(k)=mode(temp(temp>0));
        end
    end
end

%make this go somewhere later with new scoring of events
%Etype = nan(length(BB_i_start),1);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%%
if sum(strcmp(ChannelsList,'SnoreDB')==1)
    SnoreDB = SigT.SnoreDB;
    SnoreRMSsq = 10.^((SnoreDB/10))*(0.00002^2);
    keepcols = [1 2 3 7 8];
    SnoreOut = SignaltoBreathStats(SnoreRMSsq,BB_i_start,BB_i_mid,BB_i_end,'SnoreDB',keepcols);
    SnoreOut{:,:}=10*(log10(SnoreOut{:,:}/(0.00002^2))); %convert back to dB after averaging
    %10*log10((Smpa/1000).^2./(0.00002^2)) %10*(log10(mean_SnoreRMSsq_sleep/(0.00002^2)));
end
if sum(strcmp(ChannelsList,'NoxAudio')==1)||sum(strcmp(ChannelsList,'Audio_Volume')==1)
    NoxAudio = SigT.NoxAudio; %mPa %to do, add option for when signal is called Audio_Volume
    SnoreDBNox = 10*log10((NoxAudio).^2./(0.00002^2)); %10*(log10(mean_SnoreRMSsq_sleep/(0.00002^2)));
    SnoreRMSsqNox = 10.^((SnoreDBNox/10))*(0.00002^2);
    keepcols = [1 2 3 7 8];
    SnoreOutNox = SignaltoBreathStats(SnoreRMSsqNox,BB_i_start,BB_i_mid,BB_i_end,'SnoreDBNox',keepcols);
    SnoreOutNox{:,:}=10*(log10(SnoreOutNox{:,:}/(0.00002^2))); %convert back to dB after averaging
end
if 1 && exist('SnoreDB')||exist('SnoreDBNox')
    figure(89); clf(89); set(gcf,'color',[1 1 1])
     ax(1)=subplot(3,1,[1:2]);
    if exist('SnoreDBNox')
        plot(Time,SnoreDBNox);
        hold on
    end
    if exist('SnoreDB')
        plot(Time,SnoreDB);
    end
    hold off
    box off
    set(gca,'xtick',[],'xcolor',[1 1 1]);
    ylabel('Snore dB');
    ax(2)=subplot(3,1,3);
    plot(Time,Vflow_out);
    box off
    linkaxes(ax,'x');
end

%% RIP features, TBD
if sum(strcmp(ChannelsList,'kFlow')==1)
    kFlow = SigT.kFlow;
end
if exist('Ab')&&exist('Th')
RIPcorr=nan(length(BB_i_start),1);
for i=1:length(BB_i_start)
    I = BB_i_start(i):BB_i_end(i);
    %Thvar = Th(I) - detrend(Th(I))
    RIPcorr(i) = corr(Th(I),Ab(I));
end
RIPtable = table(RIPcorr);
end

%% Derive Ventilation from thorax and abdomen RIP bands:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
RIP_Thorax=Th;
RIP_Abdo=Ab;
if 0
    RIP_Thorax = filtfilt(B_butter1,A_butter1,RIP_Thorax);
    RIP_Abdo = filtfilt(B_butter1,A_butter1,RIP_Abdo);
end

%check for poor quality RIP trace
RIPinfo='';
if abs(corr(RIP_Thorax,RIP_Abdo))<0.3
    %there is a low quality trace; duplicate the good one.
    if abs(corr(RIP_Abdo,Vflow_out))>abs(corr(RIP_Thorax,Vflow_out))
        RIP_Thorax=RIP_Abdo;
        RIPinfo='Poor Th';
    else
        RIP_Abdo=RIP_Thorax;
        RIPinfo='Poor Ab';
    end
end

VEThorax=zeros(1,length(BB_i_start));
VEAbdo=zeros(1,length(BB_i_start));
for k=1:length(BB_i_start)
    VEThorax(k)=max(RIP_Thorax((BB_i_start(k):BB_i_end(k))))-min(RIP_Thorax((BB_i_start(k):BB_i_end(k))));
    VEAbdo(k)=max(RIP_Abdo((BB_i_start(k):BB_i_end(k))))-min(RIP_Abdo((BB_i_start(k):BB_i_end(k))));
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Take ventilation as the average of inspired/expired levels,
% and find RIP thorax to abdo ratios
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
VI=(VI+VE)/2; % Mean of inspiration and expiration volumes
%note that currently VI and VE are already the same and equal VT/Ttot
meanVIbeforenormalizing = mean(VI);
VI=VI/mean(VI); %normalize
VEThorax=VEThorax/mean(VEThorax); %normalize
VEAbdo=VEAbdo/mean(VEAbdo); %normalize
VEThorax(VEThorax>3)=3; %limit extremes
VEAbdo(VEAbdo>3)=3; %limit extremes
VEThorax=VEThorax/mean(VEThorax); %normalize thorax RIP again
VEAbdo=VEAbdo/mean(VEAbdo); %normalize abdo RIP again

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% Attempt to get the Stable Breathng data
StableBB_Channel = find(strcmp(ChannelsList,'StableBreathing')==1);
if isfinite(StableBB_Channel)
    try
        StableBreathingFull = SigT.StableBreathing;
        StableBreathingBB = StableBreathingFull(BB_i_mid);
    catch
        StableBreathingBB = NaN(length(BB_i_mid),1);
    end
else
    StableBreathingBB = NaN(length(BB_i_mid),1);
end

if 0
    figure(1001); clf(figure(1001)); fig = gcf;
    fig.Color = [ 1 1 1 ]; fig.Units = 'Inches';
    fig.Position = [27 1.5 13 3];
    plot(StableBreathing); hold on;    
    ylim([-1 max(StableBreathing)+1]);
end

%% Additional breath info
%pos_B = round(positionfactor*DataEventHypnog_Mat(BB_i_start(:),ColumnHeads(10)));
pos_B = PositionSig(BB_i_start);
pos_B = pos_B(:);
hypnog_B = SigT{BB_i_start(:),find(strcmp(ChannelsList,'Epochs')==1)};
hypnog_B = hypnog_B(:);
Time0 = Time(1)+0*BB_i_start;
IEratio = IEratio+0*BB_i_start;
IEratioEstimated = IEratioEstimated+0*BB_i_start;
leak_A = leak+0*BB_i_start;
FlowDriftAmplitudeA_ = FlowDriftAmplitudeA+0*BB_i_start;
FlowDriftAmplitudeB_ = FlowDriftAmplitudeB+0*BB_i_start;

BreathDataList = {'Time0','Time_start','Time_mid','Time_end','BB_i_start','BB_i_mid','BB_i_end',...
    'VI','AR','ARF','ARei','AReiF','spo2','pos_B','hypnog_B','Etype','DeltaPes','DeltaPmus','DeltaEdi','PeakEdi','VIpes','VIedi',...
    'GGpeak','GGtonic','FlowPes_VI','FlowEdi_VI','VE','Veup',...
    'Pdelta','Ptheta','Palpha','Pbeta','WakeSleep','ApneaB',...
    'FlowDriftAmplitudeA','FlowDriftAmplitudeB',...
    'FclipIB', 'FclipEB','VIOrig','REMphasic','REMtonic','REMconjugate','REMconjugateFThres','StableBreathing','Ti','Te','Ttot','IEratio','IEratioEstimated','leak','leak2'};
BreathDataTable=table(Time0, Time(BB_i_start), Time(BB_i_mid), Time(BB_i_end), BB_i_start, BB_i_mid, BB_i_end, ...
    VI(:), AR(:), ARf(:),ARei(:),AReif(:),spo2(BB_i_start), pos_B(:), hypnog_B(:), Etype(:), DeltaPes(:), DeltaPmus(:), DeltaEdi(:), PeakEdi(:), VIpes(:), VIedi(:), ...
    GGpeak, GGtonic, FlowPes_VI, FlowEdi_VI, meanVIbeforenormalizing*VI(:), meanVIbeforenormalizing+0*AR(:),...
    PdeltalogfiltB(:),PthetalogfiltB(:),PalphalogfiltB(:),PbetalogfiltB(:),WakeSleepB(:),Apnea_B,...
    FlowDriftAmplitudeA_(:),FlowDriftAmplitudeB_(:),...
    FclipIB,FclipEB,VIOrig(:),REMphasic(:),REMtonic(:),REMconjugate(:),REMconjugateFThres(:),StableBreathingBB(:),...
    Ti(:),Te(:),Ti(:)+Te(:),IEratio(:),IEratioEstimated(:),leak_A(:),leak_B(:),...
    'VariableNames', BreathDataList);

if isfield(settings,'intentionalclipping')==1 && settings.intentionalclipping==1
    BreathDataTable = [BreathDataTable BreathDataTableClipped];
end

if exist('FlowDriveOut')
    BreathDataTable = [BreathDataTable FlowDriveOut];
end
if exist('SnoreOut')
    BreathDataTable = [BreathDataTable SnoreOut];
end
if exist('SnoreOutNox')
    BreathDataTable = [BreathDataTable SnoreOutNox];
end
if exist('RIPtable')
    BreathDataTable = [BreathDataTable RIPtable];
end
if exist('PepiBreathDataTable')
    BreathDataTable = [BreathDataTable PepiBreathDataTable];
end
if exist('FlowCompartments')
    BreathDataTable = [BreathDataTable FlowCompartments];
end
if exist('FotT')
    BreathDataTable = [BreathDataTable FotT];
end

if exist('StimON')
    BreathDataTable = [BreathDataTable StimONtbl];
end

if exist('HGNSon')
    BreathDataTable = [BreathDataTable HGNSonTbl];
end

if exist('StimAmpl')
    BreathDataTable = [BreathDataTable StimAmplTbl];
end

if exist('StimSelect')
    BreathDataTable = [BreathDataTable StimSelectTbl];
end
%% Noise Analysis for flow shape
%
if 1
    BreathDataTable = SignalToNoiseWin(BreathDataTable,Flow-leak, Time);
end

%% Inverted flow reporting
% EAS added on 2021-12-01 to report inverted flow windows
if (isfield(settings,'reportInvertedFlowWindows') && settings.reportInvertedFlowWindows==1)
    listNames = {'VTi','VTe','VT','InvertedFlow'};
    data=table(VTi(:),VTe(:),VT(:),InvertedFlow(:),'VariableNames', listNames);
    BreathDataTable = [BreathDataTable data];
end

%% Export local Flow and Edi signals
switch settings.GetLocalSignals
    case 'FlowEdi'
        LocalSignals=[Vflow_out,Edi];
    case 'FlowOnly'
        LocalSignals=[Vflow_out];
    case 'BreathSignature'
        %Create BreathMatrixResampledFlow
        BreathSignature = nan(length(BB_i_start),500);
        for i = 1:length(BB_i_start)
            inspTemp = Vflow_out(BB_i_start(i):BB_i_mid(i));
            expTemp = Vflow_out(BB_i_mid(i):BB_i_end(i));
            offset = inspTemp(1);
            inspZero = (inspTemp - offset);
            expZero = (expTemp - offset);
            inspNorm = inspZero/max(inspZero);
            expNorm = expZero/max(inspZero);
            inspLen = length(inspNorm);
            expLen = length(expNorm);
            inspRS = resample(inspNorm, 250, inspLen);
            expRS = resample(expNorm, 250, expLen);
            BreathSignature(i,1:250) = inspRS';
            BreathSignature(i,251:500) = expRS';
        end
        LocalSignals=[BreathSignature];
    case 'BreathSignature2'
        %Create BreathMatrixResampledFlow; not normalized yet, whole breath signal is here
        BreathSignature = nan(length(BB_i_start),500);
        for i = 1:length(BB_i_start)
            inspTemp = Vflow_out(BB_i_start(i):BB_i_end(i));
            offset = inspTemp(1);
            inspTemp2 = (inspTemp - offset);
            inspLen = length(inspTemp2);
            inspRS = resample(inspTemp2, 500, inspLen);
            BreathSignature(i,1:500) = inspRS';
        end
        LocalSignals=[BreathSignature];        
    otherwise
        LocalSignals=NaN;
end

%% Break if only running FlowShapes
if settings.flowshapesonly==1||settings.flowshapesonly==2
    LoopGainParameters(1:17)=NaN;
    EventsInfo(1:10)=NaN;
    LG_QualityInfo(1:13)=NaN;
    DataOut=NaN;
    ArousalDat=NaN;
    FitQual=NaN;
    MiscEpochData=NaN;
    disp('exiting 3: skipping loop gain analysis');
    return
end

if settings.verbose
    delta_t = etime(clock, clocktimestart); % delta in seconds
    D_win = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
    str = ['Analysis of window complete: ', char(D_win), ' (hh:mm:ss)']; disp(str);
end

% 1 Time(BB_i_start)
% 2 Time(BB_i_mid)
% 3 Time(BB_i_end)
% 4 VI'
% 5 Vdr_est
% 6 E1'
% 7 E_recover'
% 8 E_Terminate'
% 9 Error
% 10 AR
% 11 meanVIbeforenormalizing+0*AR
% 12 VAr_est
% 13 pos_B
% 14 hypnog_B
% 15 meanVIbeforenormalizing*VI'
% 16 DeltaPes
% 17 DeltaPmus
% 18 VIpes
% 19 DeltaEdi
% 20 VIedi
% 21 GGpeak
% 22 GGtonic
% 23 FlowPes_VI
% 24 FlowEdi_VI];
end

