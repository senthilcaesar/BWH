  %% MIFL - Method for Identifiying Flow Limitation
% This function is a composite of many methods identifying flow limitation.
% It is designed to work on a window of data, stepping through each breath.
%
% Inputs to this function:
%  Time, the time series of the window data
%  Flow, the flow data for the window
%  Breath Timing, the start, mid and end points of each breath
%
% Outputs from this function:
%  BreathDataTable, each line is a breath, and each column is a feature
%
% The objective is that these features should provide some direction on the
% severity of pharyngeal airway obstruction
%
% The basic process is to split up the data series, and work on each breath
% Once a breath is segmented, we perform a few basic tests on the original
% data, and then normalise the breath (max flow = 1), and perform a range
% of additional tests that require normalised data. While the orignal data
% provides information of breath amplitude etc, this may be misleading
% when comparing shapes, and as such, the normalised signal should remove
% bias that may arise due to signal variation (e.g. if you multiply the
% incoming signal by two, that shouldn't change the resutls from any of the
% shape based features).
%
% v1 - original method
% v2 - designed to be run with either orignal or modified timing
% v3 - added VTe/VT "feature"
%    - added breath offset options
%    - moved flutter to work on window, not breath
% v4 - RemoveFlutter turned off (simplify)
%    - Excluding flutter features (output for these ftrs set to NaN)
% v5 - dedicated tertile peak detection for SS_Area, no longer uses AApeaks
%    - scaled Ellipse, Hypcos and Parab area to match other area measures
%    - 
% v6 - explicit versioning removed when added to git repository
%    - updates to Ali peaks features (where often highly NaN'd)
%
% search for "ToDo" to find unresolved issues.
%
%% MIFL - Method for Identifiying Flow Limitation
function [BreathDataTable] = MIFL(Time, Flow, BB_original, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B, MIFLsettings)
%% options
RemoveFlutter = 0; % set to 1 to remove flutter, 0 for unadjusted flow
ShowFigures = 0; % set as 1 to enable figures for each breath, 0 for quiet operation
ShowSummaryFigures = 0; % set as 1 to enable summary figures, 0 for quiet operation
PublicationFigure = 0; % set as 1 to enable publication style figure (showing key features), 0 for quiet operation
Verbose = 0; % set as 1 to get debug info during processing, 0 for quiet
ShowErrors = 1; % set as 1 to get error info during processing, 0 for quiet
downsampledFs = MIFLsettings(1); % set as 0 to not change Flow, or set as downsample Hz
useOriginalTiming = MIFLsettings(2); % set to 1 to use original timing, 0 for mod timing
if length(MIFLsettings)>2enough
    newFsrate = MIFLsettings(3); % change the Fs to artificially shorten/lengthen breaths, keeping dt
else
    newFsrate = 0;
end
IgnoreApneaBreaths = 1; % set as 1 to skip apnea breaths, 0 to attempt analysis
signal_adjustment =  'OverallMinMax';%  'SingleOffset';% 'None';% 'BaselineChord';%
IncludeFlutterFtrs = 1; % set to 1 to include flutter based features, 0 to exclude these ftrs

% list of figures used in this function
% Summary Figures (window length)
% 50 : window length flow and downsampled flow and BB timing marks
% 51 : 
% 52 : 
% Individual Figures (per breath)
% 53 : split insp & exp flow, and points of interest, incl peaks, sine fit, 
%       morgensterns, Insp flow and vol, VPEF / VTi, and ATS exp flow limit (VTe), 
%       Exp flow and vol, SS area
% 54 : teschler
% 55 : plot for publication, shows complete breath (insp+exp) flow and vol
%       with featues overlay
% 56 : other curve fitting methods (Elipse, Inv Parabola, HypCos)
% 57 : histogram with data representing insp flow
% 58 : histogram with data representing exp flow
% 59 : mansour, flow top panel, derivative of polyfit bottom panel
% 60 : wiriyaporn expiratory decy model
% 61 : series, incl flow volume loop
% 62 : 
% 63 : SSarea

% close all the figures this function generates
if 0
    for j = 53:63
        close(figure(j));
    end
end

%% signal pre-processing
%transpose flow signal if needed
if size(Flow,2)<size(Flow,1); Flow=Flow'; end
%transpose time signal if needed:
if size(Time,2)<size(Time,1); Time=Time'; end

dt=(Time(end)-Time(1))/(length(Time)-1);
Fs = round(1/dt);

%% resample (currently unused here because this is done in Analysis.m)
% may be coming back to here, as we want Analysis.m to use full sample rate
Flow_original = Flow; % for plot
if downsampledFs~=0
    originalFs = Fs;
    % downsample
    [~,q] = rat(downsampledFs / originalFs);
    Flow_ds = resample(Flow,1,q);
    % upsample
    [p,~] = rat(originalFs / downsampledFs);
    Flow = resample(Flow_ds,p,1);
end
% p is used later, so clear to be sure
clear p q

%% testing the effects of artificially longer/shorter breaths
% this resets the Fs, but keeps the dt. So the interval is the same, but it
% appears as though the breath is longer, or shorter. 
% BB times (indexes) are adjusted accordingly.
% This will break (at very least cause different results) in FFT features
if newFsrate~=0
    FsFactor = newFsrate / Fs;
    Flow=interp1(linspace(1,length(Flow)*FsFactor, length(Flow)),Flow,1:length(Flow)*FsFactor,'spline'); 
    Flow_original = Flow; 
    Time = 0:dt:length(Flow)/Fs;            % make new "Time"
    BB_s = BB_original; % this is for plots below
    BB_original = BB_original .* FsFactor;  % realign the breath timing  ToDo: may need 'round' if not integer adjustment
    BB_Ttrans = BB_Ttrans .* FsFactor;
    TiTrans = TiTrans .* FsFactor; %these are calculated before this fn, so scale to keep magnitude
    TeTrans = TeTrans .* FsFactor;
end

if 0 % looking at the testing of the duration change
SetLength = 1000;
Flow_set = interp1(linspace(1,SetLength, length(Flow)),Flow,1:SetLength,'spline');
BB_flow_set = round(BB_original .* (SetLength/length(Flow)));
Flow_original_set = interp1(linspace(1,SetLength, length(Flow_original)),Flow_original,1:SetLength,'spline');
BB_flow_original_set = round(BB_s .* (SetLength/length(Flow_original)));
isequal(BB_flow_set, BB_flow_original_set)

figure(1); clf(figure(1));
plot(1:1:SetLength, Flow_set, 'r'); hold on;
plot(1:1:SetLength, Flow_original_set, 'g');
plot((BB_flow_set(:,1)), Flow_set(BB_flow_set(:,1)), 'ro');%'markersize',10);
plot((BB_flow_set(:,2)), Flow_set(BB_flow_set(:,2)), 'rx');%'markersize',10);
plot((BB_flow_set(:,3)), Flow_set(BB_flow_set(:,3)), 'r*');
end

%% Remove flutter from full window 
% works less well on indiv breaths, doesn't work at all on insp or exp only
Flow_DSwithFlutter = Flow; % keep a backup of the original flow signal
if RemoveFlutter
    Flutter = FindWaveLetCoefs(Flow,'db4',Fs,'Flow Noise');
    Flow = Flow - Flutter;
end

%% Set breath timing
% The idea is that MIFL function can work on either original timing, or the 
% Ttran modified timing, depending on what version of breath timing is 
% passed in. TiTrans and TeTrans are always attempted, using whichever
% timing is currently used. 
if useOriginalTiming
    BB_i_start = BB_original(:,1);
    BB_i_end = BB_original(:,2);
    BB_e_start = BB_original(:,2);
    BB_e_end = BB_original(:,3);   
else 
    BB_i_start = BB_Ttrans(:,1);
    BB_i_end = BB_Ttrans(:,2);
    BB_e_start = BB_Ttrans(:,3);
    BB_e_end = BB_Ttrans(:,4);
end

%% Make the Features list and all of the associated variables
BBs=length(BB_i_start);
FeatureNames = MakeFeaturesList(BBs);
% Ignore messages about variables changing size on every loop iteration.
% The variables are intialised at run time, by the function above, 
% and as such, will have been pre-allocated when they are used below.

if ShowSummaryFigures
    figure(50); clf(figure(50));
    plot(Time, Flow_original, 'g'); hold on;  % original flow
    plot(Time, Flow_DSwithFlutter,'r');       % down sampled flow
    plot(Time, Flow, 'k');                    % flow w flutter removed
    plot(Time(BB_Ttrans(:,1)), Flow(BB_Ttrans(:,1)), 'm^');
    plot(Time(BB_Ttrans(:,2)), Flow(BB_Ttrans(:,2)), 'm^');
    plot(Time(BB_Ttrans(:,3)), Flow(BB_Ttrans(:,3)), 'mv');
    plot(Time(BB_Ttrans(:,4)), Flow(BB_Ttrans(:,4)), 'mv');
    plot(Time(BB_original(:,1)), Flow(BB_original(:,1)), 'ro');%'markersize',10);
    plot(Time(BB_original(:,2)), Flow(BB_original(:,2)), 'rx');%'markersize',10);
    plot(Time(BB_original(:,3)), Flow(BB_original(:,3)), 'r*');
    refline(0,0);   
%     plot(Time(BB_i_start), Flow(BB_i_start), 'go');
%     plot(Time(BB_i_end), Flow(BB_i_end), 'g^');
%     plot(Time(BB_e_start), Flow(BB_e_start), 'rv');
%     plot(Time(BB_e_end), Flow(BB_e_end), 'rv');
    for i=1:length(BB_original)
        str = num2str(i);
        text(Time(BB_original(i,1)), nanmean(Flow(Flow>0)), str ,'FontSize',8, 'Color','b'); % print text on plot at loc
    end
    legend('Original','Downsampled', 'DeFluttered');
    title('Flows and Timings');    
end

warning('off','all');

settingstext = ['[',num2str(MIFLsettings),'] ']; % this is just used in logging and error reporting

%% Step through each breath, and determine values for all of the features
% this is usually run as a parfor loop, however, currently when run as
% parfor, it returns the following (hard to resolve) error:
% % Warning: X is rank deficient to within machine precision.
% % > In regress (line 84)
% %   In parallel_function>make_general_channel/channel_general (line 914)
% %   In remoteParallelFunction (line 38)
for i = 1:BBs  
    if Verbose
        displaytext=['Analyzing breath ' num2str(i) '/' num2str(BBs)];
        disp(displaytext);
    end
    
    if isnan(BB_i_start(i))
        % if we don't have timing for a breath, we can't analyse it.
        if ShowErrors
            displaytext=[settingstext 'No timing data for breath ' num2str(i) '/' num2str(BBs)];
            disp(displaytext);
        end
        continue
    end
   
    if Apnea_B(i)
        % this breath has been determined in LGfromFlow as apneic
        % the results of MIFL for this breath should be treated carefully
        if IgnoreApneaBreaths
            if Verbose
                displaytext=[settingstext 'Ignoring Apnea breath ' num2str(i) '/' num2str(BBs)];
                disp(displaytext);
            end
            continue
        end
    end
   
    %% set the breath indices for the current breath
    BB_I_start = 1;
    BB_I_end = BB_i_end(i) - BB_i_start(i);
    BB_E_start = BB_e_start(i) - BB_i_start(i);
    BB_E_end = BB_e_end(i) - BB_i_start(i);
    
    if BB_I_end-BB_I_start < (Fs/5) % at least 0.2 of a second, probably could increase this... 
        if ShowErrors
            displaytext=[settingstext 'Short Inspiration, skipping breath ' num2str(i) '/' num2str(BBs)];
            disp(displaytext);
        end
        continue
    end
    
    if BB_E_end-BB_E_start < (Fs/5)
        if ShowErrors
            displaytext=[settingstext 'Short Expiration, skipping breath ' num2str(i) '/' num2str(BBs)];
            disp(displaytext);
        end
        continue
    end
    
    %% set up the flow variables to use throughout most processing
    flow = Flow(BB_i_start(i):BB_e_end(i));% get the full flow data for the breath (insp and exp)
    flow_original = Flow_original(BB_i_start(i):BB_e_end(i));% get the full flow data for the breath (insp and exp)
    flow_DSwithFlutter = Flow_DSwithFlutter(BB_i_start(i):BB_e_end(i));% get the full flow data for the breath (insp and exp)
    %time_actual = Time(BB_i_start(i):BB_i_end(i));
    time =(0:length(flow)-1)/Fs; % make a time vector starting from zero for this breath
   
    % split in to Insp and Exp portions (and also get the original signal)
    InspFlow_original = flow_original(BB_I_start:BB_I_end); % original inspiratory flow
    InspFlow_DSwithFlutter = flow_DSwithFlutter(BB_I_start:BB_I_end); % original inspiratory flow
    InspFlow = flow(BB_I_start:BB_I_end); % flutter removed inspiratory flow
    InspTime = time(BB_I_start:BB_I_end); % inspiratory time
    
    ExpFlow_original = flow_original(BB_E_start:BB_E_end); % original expiratory flow
    ExpFlow_DSwithFlutter = flow_DSwithFlutter(BB_E_start:BB_E_end); % original expiratory flow
    ExpFlow = flow(BB_E_start:BB_E_end); % flutter removed expiratory flow
    ExpTime = time(BB_E_start:BB_E_end); % expiratory time, note in Ttran, will not be continuous with Insp time
    
    %% baseline adjustment
    % 'SingleOffset'  - use the min of (insp start or insp end points),
    %                   and max of (exp start or exp end points).
    % 'BaselineChord' - subtracts a chord between start and end insp, and
    %                   then same for exp. works well in clean breaths with 
    %                   start and end timing already near zero crosses
    % 'OverallMinMax' - use the overall min of insp and max of exp,
    %                   good if breath scoops are > than start/end markers,
    %                   this can happen with scoopy breaths and start/end 
    %                   times are high on the flow (worse in Ttran timing)  
    % 'None'          - do no offsets of adjustments to signal
    switch signal_adjustment
        case 'SingleOffset' % offset by the min(insp) or max(exp)
            % this is problematic if the start and end times are a little
            % way up to flow signal, and the flow dips during the breath,
            % as this will cause the dipped flow to be negative.
            
            % inspiration
            InspOffset = min(InspFlow(1), InspFlow(end));
            InspFlow = InspFlow - InspOffset;
            % expiration
            ExpOffset = max(ExpFlow(1), ExpFlow(end));
            ExpFlow = ExpFlow - ExpOffset;
            
        case 'BaselineChord' % offset by a chord drawn between the start 
                             % and end of insp, and then again for exp
            % this sounds pretty neat, and in an ideal world, it would be
            % good, but it can cut through the breath if it dips.
            
            % inspiration
            x = [InspTime(1) InspTime(end)];
            y = [InspFlow(1) InspFlow(end)];
            c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector       
            chord = zeros(length(InspTime),1);
            for k=1:length(InspTime)
                chord(k) = (c(2)*(InspTime(k)))+c(1);
            end
            InspFlow = InspFlow - chord'; 
            % expiration
            x = [ExpTime(1) ExpTime(end)];
            y = [ExpFlow(1) ExpFlow(end)];
            c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector       
            chord = zeros(length(ExpTime),1);
            for k=1:length(ExpTime)
                chord(k) = (c(2)*(ExpTime(k)))+c(1);
            end
            ExpFlow = ExpFlow - chord'; 
        case 'OverallMinMax' 
            % this seems most realiable. no frills, nothing exciting.
            InspOffset = min(InspFlow);
            InspFlow = InspFlow - InspOffset;
            % expiration
            ExpOffset = max(ExpFlow);
            ExpFlow = ExpFlow - ExpOffset;
        otherwise
            % do nothing to adjust the signal
            % use it as it is.
    end  
        
    %% normalise to max inspiratory flow
    FlowNorm = max(InspFlow);
    InspFlow = InspFlow/FlowNorm;
    ExpFlow = ExpFlow/FlowNorm;
    %flow = flow/FlowNorm;
    flow = flow/max(flow);
    flow_reconstituted = [InspFlow,ExpFlow];
    
    %% Force inspiration as positive, and expiration as negative
    %InspFlow(InspFlow<0)=0.000001; % force postive, or  near zero
    InspFlow = max(InspFlow, 0); % force postive or zero
    %ExpFlow(ExpFlow>0)=-0.0000001; % force negative, or near zero
    ExpFlow = min(ExpFlow, 0); % force postive or zero

    %% normalize duration
    % if we want to ensemble average breaths, we can scale the length of the breath to a fixed value.
    % we also scale the length to normalize duration, as this effects some features
    SetInspLength = 200; % 
    SetExpLength = 250;  %
    InspFlow_SetLength=interp1(linspace(1,SetInspLength, length(InspFlow)),InspFlow,1:SetInspLength,'spline'); 
    ExpFlow_SetLength=interp1(linspace(1,SetExpLength, length(ExpFlow)),ExpFlow,1:SetExpLength,'spline');
    
    %testing
    if 0
        M1(i) = mean(InspFlow_SetLength);
        M2(i) = median(InspFlow_SetLength);
        M3(i) = std(InspFlow_SetLength);
        M4(i) = skewness(InspFlow_SetLength);
        M5(i) = kurtosis(InspFlow_SetLength);
    
        M1(i) = mean(InspFlow);
        M2(i) = median(InspFlow);
        M3(i) = std(InspFlow);
        M4(i) = skewness(InspFlow);
        M5(i) = kurtosis(InspFlow);
    end
    
    if 0
        figure();
        subplot(2,1,1);
        plot(InspTime, InspFlow);
        subplot(2,1,2);
        plot(1:SetInspLength, InspFlow_SetLength);
    end
    
    if ShowFigures
        figure(53);
        clf(figure(53));
        subplot(2,1,1);
        plot(InspTime, InspFlow, 'k'); hold on;
        plot(InspTime, InspFlow_original, 'b-.');
        refline(0,0);
        legend('Adjusted','Original'); box off;
        ylabel('InspFlow (normalised to max flow)'); xlabel('Time (seconds)');
        ax = gca;
        ax.TickDir = 'out';
        subplot(2,1,2);
        plot(ExpTime, ExpFlow, 'k'); hold on;
        plot(ExpTime, ExpFlow_original, 'b-.');
        refline(0,0);
        ylabel('ExpFlow (normalised with InspFlow)'); xlabel('Time (seconds)');
        legend('Adjusted','Original', 'Location','Southeast'); box off;
        fig = gcf;
        fig.Color = [1 1 1]; % set background colour to white
        fig.Units = 'inches';
        fig.Position = [1.5 1.5 6 8];
        ax = gca;
        ax.TickDir = 'out';
    end
    if PublicationFigure 
        TickFntSz = 12;
        
        figure(55);
        clf(figure(55));
        
        fig = gcf;
        fig.Color = [1 1 1]; fig.Units = 'inches';       
        fig.Position = [1.5 1.5 10 8];

        ax(1)=subplot(4,1,[1:2]);
        plot(time, flow, 'k'); hold on;
        rline1 = refline(0,0);
        rline1.Color = [0.8 0.8 0.8];
        set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
        set(gca,'xtick',[],'box','off');
        set(gca,'xcolor',[1 1 1])
        ylim([-1.45 1.45]);
        ylabel('Flow (normalised to max flow)'); xlabel('Time (seconds)');

        ax(4)=subplot(4,1,4);
        VolForFig = cumsum(flow).*dt;
        VolForFig = VolForFig/max(VolForFig);
        plot(time, VolForFig, 'k'); hold on;
        rline2 = refline(0,0);
        rline2.Color = [0.8 0.8 0.8];
        set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
        set(gca,'xtick',[],'box','off');
        set(gca,'xcolor',[1 1 1])
        ylim([-0.1 1.1]);
        ylabel('Volume (normalised to 1)'); xlabel('Time (seconds)');
        
    end
    
    %% Ali CalcAsymmetry
    AsymIndex(i)=CalcAsymmetry(flow,ExpFlow,Fs);
        
    %% >>>>>>>>>>>>>> Set up some timings and flows <<<<<<<<<<<<<<<<<<<
    %vol=cumtrapz(time,flow_original); %figure();plot(vol);
    % we do the insp and exp volumes separately. 
    % Ttran timing has gaps between insp and exp segments.
    
    InspVol = cumtrapz(InspTime, InspFlow);
    ExpVol = cumtrapz(ExpTime, ExpFlow);
    
    %VTi = vol(BB_I_end) - vol(BB_I_start);
    VTi = max(InspVol(end), max(InspVol));
    
    %VTe = vol(BB_E_start) - vol(BB_E_end);
    VTe = -1*(min(ExpVol(end), min(ExpVol)));
    
    %Ti = (BB_I_end-1)/Fs;
    Ti = length(InspTime)/Fs;
    
    %Te = (length(vol)-BB_E_start)/Fs;
    Te = length(ExpTime)/Fs;
    
    VT = ((VTi*Te)+(VTe*Ti))/(Ti + Te);  
    
    if useOriginalTiming
        % Ttot for original timing is just Ti + Te
        Ttot = Ti + Te;
    else
        % Ttot for Ttran timing is Ti+TiTran+Te+TeTran
        if i==BBs % last breath
            % we don't have TeTrans for the last breath (because we don't
            % know when the next one starts). BUT, we don't want to add lots
            % of NaN's to the data, so we will set Ttot for the last breath
            % in each window as Ti+TiTran+Te+median(TeTran)
            Ttot = Ti + TiTrans(i) + Te + (nanmedian(TeTrans));
        else
            % normal breath in middle of analysis window
            Ttot = Ti + TiTrans(i) + Te + TeTrans(i);
        end
    end
    
    % peak flows and time to peak flows
    [PIF,TPIF_i] = max(InspFlow); % PIF should be one... because scaled above
    TPIF = TPIF_i*dt;
    [PEF,TPEF_i] = min(ExpFlow);
    TPEF = TPEF_i*dt;
    
    % average flows
    MIF=mean(InspFlow);
    MEF=mean(ExpFlow);
       
    % indices for proportions through breath time
    % the index location of 10,25,50,75,90 and 95% through insp
    % Ti05 is only used in Discontinuity measures
    Ti05 = round(0.05*(BB_I_end),0); Ti05_set = round(0.05*(SetInspLength),0);    
    if Ti05==0; Ti05=1; end; if Ti05_set==0; Ti05_set=1; end 
    Ti10 = round(0.10*(BB_I_end),0); Ti10_set = round(0.10*(SetInspLength),0);
    Ti25 = round(0.25*(BB_I_end),0); Ti25_set = round(0.25*(SetInspLength),0);
    Ti30 = round(0.30*(BB_I_end),0); Ti30_set = round(0.30*(SetInspLength),0);
    Ti50 = round(0.50*(BB_I_end),0); Ti50_set = round(0.50*(SetInspLength),0);
    Ti75 = round(0.75*(BB_I_end),0); Ti75_set = round(0.75*(SetInspLength),0);
    Ti90 = round(0.90*(BB_I_end),0); Ti90_set = round(0.90*(SetInspLength),0);
    Ti95 = round(0.95*(BB_I_end),0); Ti95_set = round(0.95*(SetInspLength),0);
    % the index location of 20,25,50,75 and 80% through exp  
    Te20 = round(0.20*(BB_E_end-BB_E_start),0); Te20_set = round(0.20*(SetExpLength),0);
    Te25 = round(0.25*(BB_E_end-BB_E_start),0); Te25_set = round(0.25*(SetExpLength),0);
    Te50 = round(0.50*(BB_E_end-BB_E_start),0); Te50_set = round(0.50*(SetExpLength),0);
    Te75 = round(0.75*(BB_E_end-BB_E_start),0); Te75_set = round(0.75*(SetExpLength),0);
    Te80 = round(0.80*(BB_E_end-BB_E_start),0); Te80_set = round(0.80*(SetExpLength),0);

    if ShowFigures
        figure(53); 
        subplot(2,1,1); hold on;
        plot(InspTime(TPIF_i),PIF, 'r^');
        plot(InspTime, InspVol, 'r-.');
        legend 'off';
        subplot(2,1,2); hold on;
        plot(ExpTime(TPEF_i),PEF, 'rv');
        plot(ExpTime, ExpVol, 'r-.');
        legend 'off';
    end

    %% >>>>>>>>>>>>>>>>>>>>>>  RATIOS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Ti_Ttot(i) = Ti/Ttot;
    VTi_VTe(i) = VTi/VTe;
    VTi_VT(i) = VTi/VT;
    VTe_VT(i) = VTe/VT;
    MIF_PIF(i)=MIF/PIF; % PIF should be one, therefore = MIF
    PIF_MIF(i)=PIF/MIF; % PIF should be one, therefore = 1 / MIF
    %MIF_VT(i)=MIF/VT;
    %PIF_VT(i)=PIF/VT;
    %PEF_VT(i)=PEF/VT;    
    
    %% >>>>>>>>>>>>>>>>>>>>  FLATTENING  <<<<<<<<<<<<<<<<<<<<<<<<
    % This collection of functions is aimed at directly measuring flatness
    % or the degree of flattening in inspiration and expiration
    
    % Ali's methods (these are quite similar in principle to Teschler)
    InspFlow_Ali=InspFlow(Ti10:Ti90);
    %InspFlatnessIndex_ = trapz(abs(InspFlow_Ali-MIF))*dt;
    InspFlatnessIndex_ = trapz(InspTime(Ti10:Ti90),abs(InspFlow_Ali-MIF));
    Ali_InspFlat(i)=InspFlatnessIndex_/Ti; 

    ExpFlow_Ali=ExpFlow(Te20:Te80);
    %ExpFlatnessIndex_ = trapz(abs(ExpFlow_Ali-MEF))*dt;
    ExpFlatnessIndex_ = trapz(ExpTime(Te20:Te80),abs(ExpFlow_Ali-MEF));
    Ali_ExpFlat(i)=ExpFlatnessIndex_/Te; % 
    
    Ali_InspExpFlat(i) = Ali_InspFlat(i) / Ali_ExpFlat(i);
       
    if ShowFigures
        figure(53);
        subplot(2,1,1); hold on;
        plot(InspTime(Ti10:Ti90), MIF*ones(length(InspFlow(Ti10:Ti90)),1),'k--'); hold on;
        plot(InspTime(Ti10:Ti90), abs(InspFlow_Ali-MIF), 'c--');
        subplot(2,1,2); hold on;
        plot(ExpTime(Te20:Te80), MEF*ones(length(ExpFlow(Te20:Te80)),1),'k--'); hold on;
        plot(ExpTime(Te20:Te80), -(abs(ExpFlow_Ali-MEF)), 'c--');
    end
    
    if PublicationFigure 
        figure(55);
        subplot(4,1,[1:2]); hold on;
        plot(time(1:length(InspFlow)), MIF*ones(length(flow(1:length(InspFlow))),1),'k:'); hold on;
        text(time(length(InspFlow)+1), MIF, 'MIF');
    end
    
    % AA's methods
    AA_InspFlat9020(i)=sum(InspFlow(InspFlow>0.9))/sum(InspFlow(InspFlow>0.2));
    AA_InspFlat8020(i)=sum(InspFlow(InspFlow>0.8))/sum(InspFlow(InspFlow>0.2));
    AA_InspFlat7020(i)=sum(InspFlow(InspFlow>0.7))/sum(InspFlow(InspFlow>0.2));
    AA_InspFlat6020(i)=sum(InspFlow(InspFlow>0.6))/sum(InspFlow(InspFlow>0.2));
    AA_InspFlat5020(i)=sum(InspFlow(InspFlow>0.5))/sum(InspFlow(InspFlow>0.2));
    AA_InspFlat90Ti(i)=sum(InspFlow(InspFlow>0.9))/sum(InspFlow);
    AA_InspFlat80Ti(i)=sum(InspFlow(InspFlow>0.8))/sum(InspFlow);
    AA_InspFlat70Ti(i)=sum(InspFlow(InspFlow>0.7))/sum(InspFlow);
    AA_InspFlat60Ti(i)=sum(InspFlow(InspFlow>0.6))/sum(InspFlow);
    AA_InspFlat50Ti(i)=sum(InspFlow(InspFlow>0.5))/sum(InspFlow);
    
    % same as insp, but relative to PEF (because PIF was scaled to one)
    AA_ExpFlat9020(i)=sum(-ExpFlow(-ExpFlow>(0.9*-PEF)))/sum(-ExpFlow(-ExpFlow>(0.2*-PEF)));
    AA_ExpFlat8020(i)=sum(-ExpFlow(-ExpFlow>(0.8*-PEF)))/sum(-ExpFlow(-ExpFlow>(0.2*-PEF)));
    AA_ExpFlat7020(i)=sum(-ExpFlow(-ExpFlow>(0.7*-PEF)))/sum(-ExpFlow(-ExpFlow>(0.2*-PEF)));
    AA_ExpFlat6020(i)=sum(-ExpFlow(-ExpFlow>(0.6*-PEF)))/sum(-ExpFlow(-ExpFlow>(0.2*-PEF)));
    AA_ExpFlat5020(i)=sum(-ExpFlow(-ExpFlow>(0.5*-PEF)))/sum(-ExpFlow(-ExpFlow>(0.2*-PEF)));
    AA_ExpFlat90Te(i)=sum(-ExpFlow(-ExpFlow>(0.9*-PEF)))/sum(-ExpFlow);
    AA_ExpFlat80Te(i)=sum(-ExpFlow(-ExpFlow>(0.8*-PEF)))/sum(-ExpFlow);
    AA_ExpFlat70Te(i)=sum(-ExpFlow(-ExpFlow>(0.7*-PEF)))/sum(-ExpFlow);
    AA_ExpFlat60Te(i)=sum(-ExpFlow(-ExpFlow>(0.6*-PEF)))/sum(-ExpFlow);
    AA_ExpFlat50Te(i)=sum(-ExpFlow(-ExpFlow>(0.5*-PEF)))/sum(-ExpFlow);
    
    %AA_ExpFlat(i)=sum(-ExpFlow>0.85)/sum(-ExpFlow>0.4);
    
    % mean flow over the middle 50% of insp and exp
    MIF50(i) = mean(InspFlow(Ti25:Ti75));
    MEF50(i) = mean(ExpFlow(Te25:Te75));
     
    % Teschler(1996) flattening/curvature index
    % the old version was not normalised for duration
    % this version is normalised by Ti
    % Inspiratory airflow is scaled by mean inspiratory flow to one unit
    inspflow_forTeschler = InspFlow/mean(InspFlow); %normalised to mean flow
    % The curvature index is a measure of the deviation from unit scaled flow over the 
    % middle 50% of inspiratory time, and time is expressed as a fraction of inspiratory duration
    %Teschler(i) = (trapz(abs(inspflow_forTeschler(Ti25:Ti75)' - ones(length(InspTime(Ti25:Ti75)),1)))*dt)/Ti;
    Teschler(i) = (trapz(InspTime(Ti25:Ti75), abs(inspflow_forTeschler(Ti25:Ti75)' - ones(length(InspTime(Ti25:Ti75)),1))))/Ti;
  
    if ShowFigures && 0 % figure uses old (non-normalized) flow
        figure(54);
        clf(figure(54));
        plot(InspTime, inspflow_forTeschler,'k');hold on;
        refline(0,1);
        plot(InspTime(Ti25:Ti75), ones(length(InspTime(Ti25:Ti75))), 'r');
        %should show vert lines indicating the area being analysed
        plot(InspTime(Ti25:Ti75),...
        abs(inspflow_forTeschler(Ti25:Ti75)' - ones(length(InspTime(Ti25:Ti75)),1)),'r--');
        title('Teschler inspiration and reference line');
    end
   
    % SS_Area is included in the list as an assessment of flatness, however
    % it uses peaks as determined by AA peaks, and therefore runs later
    
    %% >>>>>>>>>>>>>>>>>>>>>>>  TEMPLATE   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    % "sine fit" using simple polynomial curve fitting
    % fit 2nd degree polynomial to three set points, i.e start, max, & end
    % expiration
    threeX = [ExpTime(1); ExpTime(round(length(ExpTime)/2,0)); ExpTime(end)];
    threeY = [ExpFlow(1); min(ExpFlow); ExpFlow(end)];
    [pp,~,mu] = polyfit(threeX, threeY, 2);
    y1 = polyval(pp, ExpTime,[],mu);
    % find area difference between flow and curve
    %SinE(i)=abs(trapz(ExpTime,ExpFlow) - trapz(ExpTime,y1));
    %SinE_ = (trapz(abs(ExpFlow-y1))*dt)/Te;
    SinE(i) = trapz(ExpTime, abs(ExpFlow-y1))/Te; 
    
    % inspiration 
    threeX = [InspTime(1); InspTime(round(length(InspTime)/2,0)); InspTime(end)];
    threeY = [InspFlow(1); max(InspFlow); InspFlow(end)];
    [pp,~,mu] = polyfit(threeX, threeY, 2);
    y2 = polyval(pp, InspTime,[],mu);
    % find area difference between flow and curve
    %SinI(i)=abs(trapz(InspTime,InspFlow) - trapz(InspTime,y2));
    %SinI(i) = (trapz(abs(InspFlow-y2))*dt)/Ti;
    SinI(i) = trapz(InspTime, abs(InspFlow-y2))/Ti; 
    SinI50(i) = trapz(InspTime(Ti25:Ti75), abs(InspFlow(Ti25:Ti75)-y2(Ti25:Ti75)))/Ti;   
    
%     figure() % just to make sure the Sin50 is working correctly
%     plot(InspTime, InspFlow); hold on;
%     plot(InspTime, y2);
%     plot(InspTime, abs(InspFlow-y2));
%     plot(InspTime(Ti25:Ti75), abs(InspFlow(Ti25:Ti75)-y2(Ti25:Ti75)));
     
    % options and starting points for next few curve fitting options
    % because we make a new x for these, we need to *dt to get back to
    % scale, and then we normalize the result by /Ti.
    
    c=max(InspFlow)*1;
    N=length(InspFlow);
    x=(1:N)-round(N/2);   %x=(1:N); 
    
    % inverted parabola
    y3=c*((1-(x.^2)/(N/2)^2));
    %InvParabI(i)=abs(trapz(x,InspFlow) - trapz(x,y3));
    InvParabI(i)=(trapz(x,abs(InspFlow-y3))*dt)/Ti; % added *dt)/Ti to scale this result
    
    % ellipseI
    y4=c*(((1-(x.^2)/(N/2)^2)).^(1/2));
    %EllipseI(i)=abs(trapz(x,InspFlow) - trapz(x,y4));   
    EllipseI(i)=(trapz(x, abs(InspFlow-y4))*dt)/Ti; % added *dt)/Ti to scale this result
    
    % hypcosI, hyperbolic cosine arch 0<b<Inf
    b=10;
    y5=c*((1-cosh((x/(N/2))*log(1+b+sqrt(2+b)*sqrt(b)))+b)/b);
    %HypcosI(i)=abs(trapz(x,InspFlow) - trapz(x,y5));
    HypcosI(i)=(trapz(x, abs(InspFlow-y5))*dt)/Ti; % added *dt)/Ti to scale this result
    
    if ShowFigures
        figure(53); hold on; 
        subplot(2,1,1);
        %plot(InspTime,InspFlow,'k');hold on
        plot(InspTime, y2, 'r');
        plot(InspTime, abs(InspFlow-y2),'r--');
        %refline(0,0);
        %title('Polyfit to inspiration');
        subplot(2,1,2);
        %plot(ExpTime,ExpFlow,'k');hold on
        plot(ExpTime, y1, 'r');
        plot(ExpTime, (-abs(ExpFlow-y1)),'r--');
        %refline(0,0);
        %title('Polyfit to expiration');
        hold off
    
        figure(56);
        clf(figure(56));
        subplot(3,1,[1:2]);
        plot(x, InspFlow,'k'); hold on;
        plot(x, y2,'c');
        plot(x, y3,'r');
        plot(x, y4,'g');
        plot(x, y5,'b');
        legend('Flow','Sine','InvP','Ellipse','HypCos');
        title('Various curve fits to inspiration');
        subplot(3,1,3);
        plot(x, abs(InspFlow-y2),'c--'); hold on;
        plot(x, abs(InspFlow-y3),'r--'); 
        plot(x, abs(InspFlow-y4),'g--');
        plot(x, abs(InspFlow-y5),'b--');
        ylim([0 1.5]);
        title('Abs diff btw curve fits and inspiration');
    end  
    
    if PublicationFigure 
        figure(55); hold on; 
        subplot(4,1,[1:2]); hold on;
        plot(InspTime, y5, 'r--');
        plot(ExpTime, y1, 'm--');
        
        ax(3)=subplot(4,1,3); hold on;
        plot(InspTime, abs(InspFlow-y5),'r');
        plot(ExpTime, abs(ExpFlow-y1),'m');
        rline3 = refline(0,0);
        rline3.Color = [0.8 0.8 0.8];
        set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
        set(gca,'xtick',[],'box','off');
        set(gca,'xcolor',[1 1 1])
        ylabel('Flow \Delta');
    end
    

    %% >>>>>>>>>>>>>>>>>>>> ASYMMETRY <<<<<<<<<<<<<<<<<<<<<<<<
    % Skewness is a measure of the asymmetry of the data around the sample mean.
    % If skewness is negative, the data are spread out more to the left of
    % the mean than to the right. If skewness is positive, the data are
    % spread out more to the right. The skewness of the normal distribution
    % (or any perfectly symmetric distribution) is zero.
    
    % Kurtosis is a measure of how outlier-prone a distribution is.
    % The kurtosis of the normal distribution is 3.
    % Distributions that are more outlier-prone than the normal distribution
    % have kurtosis greater than 3;
    % distributions that are less outlier-prone have kurtosis less than 3.
    
    %% Converting the flow shape to represent a distribution
    % and then use standard skewness and kurtosis measures
    % Inspiratory Skewness and Kurtosis
    %clear distdataInsp % can't do clear in parfor
    distdataInsp = [];
    for jj = 1:length(InspFlow_SetLength)-1
        distdataInsp = cat(1,distdataInsp,ones(round(InspFlow_SetLength(jj)*100),1)*jj);
    end
    SkewDistInsp(i)=skewness(distdataInsp);
    KurtDistInsp(i)=kurtosis(distdataInsp,0);
    
    % Expiratory Skewness and Kurtosis
    %clear distdataExp % can't do clear in parfor
    distdataExp = [];
    for jj = 1:length(ExpFlow_SetLength)-1
        distdataExp = cat(1,distdataExp,ones(round(abs(ExpFlow_SetLength(jj))*100),1)*jj);
    end
    SkewDistExp(i)=skewness(distdataExp);
    KurtDistExp(i)=kurtosis(distdataExp,0);
    
    if ShowFigures
        figure(57); clf(figure(57));
        subplot(2,1,1);
        plot(InspFlow_SetLength);
        title('Normalised inspiratory flow');
        subplot(2,1,2);
        histogram(distdataInsp,1:length(InspFlow_SetLength), 'facealpha',0.2, 'edgealpha', 0.2);
        title('Histogram of data representing insp flow');
        figure(58); clf(figure(58));
        subplot(2,1,1);
        %plot(InspTime, InspFlow);
        plot(-ExpFlow_SetLength);
        title('Normalised expiratory flow (inverted)');
        subplot(2,1,2);
        histogram(distdataExp,1:length(ExpFlow_SetLength), 'facealpha',0.2, 'edgealpha', 0.2);
        title('Histogram of data representing exp flow');
    end
    
    %% Using the flow waveform, and looking at asymmetry around mid time
    % and then standard kurtosis measure, and a home-made skewness  
    % Inspiratory and Expiratory Skewness 
    %InspFlowArea=trapz(InspFlow)*dt;
    %areaL=trapz(InspFlow(1:Ti50))*dt;
    %areaR=trapz(InspFlow(Ti50:end))*dt;
    %SkewDataInsp(i)=(areaL/InspFlowArea)-(areaR/InspFlowArea);
    try
    InspFlowArea=trapz(InspFlow_SetLength);
    areaLi=trapz(InspFlow_SetLength(1:Ti50_set));
    areaRi=trapz(InspFlow_SetLength(Ti50_set:end));
    AsymmetryInsp(i)=(areaLi/InspFlowArea)-(areaRi/InspFlowArea);
    % values near zero indicate symmetry in flow
    % positive values indicate right skew (i.e. more early insp)
    % negative values indicate left skew (i.e. more late insp)
    
    %ExpFlowArea=trapz(ExpFlow)*dt;
    %areaL=trapz(ExpFlow(1:Te50))*dt;
    %areaR=trapz(ExpFlow(Te50:end))*dt;
    %SkewDataExp(i)=(areaL/ExpFlowArea)-(areaR/ExpFlowArea);
    
    ExpFlowArea=trapz(ExpFlow_SetLength);
    areaLe=trapz(ExpFlow_SetLength(1:Te50_set));
    areaRe=trapz(ExpFlow_SetLength(Te50_set:end));
    AsymmetryExp(i)=(areaLe/ExpFlowArea)-(areaRe/ExpFlowArea);
    catch Asymmerror
        disp(Asymmerror.message); Asymmerror.getReport
    end
    % Inspiratory and Expiratory Kurtosis
    KurtDataInsp(i)=kurtosis(InspFlow_SetLength,0);
    KurtDataExp(i)=kurtosis(ExpFlow_SetLength,0);
    
    %% >>>>>>>>>>>>>>>>>>>> CURVE FITTING <<<<<<<<<<<<<<<<<<<<<<<<
    % here we attempt to fit a curve to the breath.
    % as such, these are sublty different to the template fits above
    
    %% mansourFL,
    if 0
    % fit 3rd order polynomial to insp, then find the
    % slope of the polyfit at max flow (+ve = NIFL, -ve = IFL)
    PolyOrder = 3;
    range_end = TPIF_i+0.1*(TPIF_i); % go to peak, then 10% more
    if range_end >= BB_I_end
        range_end = BB_I_end - 1;
    end    
    range=(1:range_end); 
    %range=(1:TPIF_i);
    %range=(1:BB_I_endi(i)-BB_I_starti(i));
    warning('off','MATLAB:polyfit:PolyNotUnique');
    [PolyfitValues,~,mu] = polyfit(InspTime(range),InspFlow(range),PolyOrder);
    PolyLine = polyval(PolyfitValues,InspTime(range),[],mu); % Produce a line that represents the polyfit
    PolyLine_slope = gradient(PolyLine);% find the gradient of the PolyLine
    % special case where TPIF_i > range_end
    Peaki_mod = min([TPIF_i range_end]);
    MansourFL(i) = PolyLine_slope(Peaki_mod); % Slope at Max Flow or end
    if ShowFigures
        figure(59)
        clf(figure(59))
        subplot(2,1,1)
        plot(InspTime(range),InspFlow(range)); hold on;
        plot(InspTime(range), PolyLine);
        plot(InspTime(Peaki_mod),PolyLine(Peaki_mod),'ro');
        xlim([InspTime(range(1)),InspTime(range(end))]);
        title('3^{rd} order polyfit to inspiration');
        subplot(2,1,2)
        plot(InspTime(range), PolyLine_slope); hold on;
        plot(InspTime(Peaki_mod),PolyLine_slope(Peaki_mod),'ro');
        xlim([InspTime(range(1)),InspTime(range(end))]);
        refline(0,0);
        title('Slope of polyfit');
    end
    end
    
    %% mansourRUA,
    if 0
    % ToDo: work out why this gives bad results
    % flow at min/max slope, and constant C from polynomial
    % Next, let's find Rua, as R=C/(F^2*t) where:
    %   C is the C from PolyValues above,
    %   F is the flow at the end of linear portion, as determined by
    %       min of derivative plot, if breath is unobstructed, or
    %       max of derivative plot, if breath is obstructed
    %   t is time
    if ShowFigures
        if MansourFL(i)>0 % unobstructed, therefore looking for min
            [SlopeVal,SlopeIndexOffset] = min(PolyLine_slope);
        elseif MansourFL(i)<=0 % obstructed, therefore looking for max
            [SlopeVal,SlopeIndexOffset] = max(PolyLine_slope);
        end
    else
        if MansourFL(i)>0 % unobstructed, therefore looking for min
            [~,SlopeIndexOffset] = min(PolyLine_slope);
        elseif MansourFL(i)<=0 % obstructed, therefore looking for max
            [~,SlopeIndexOffset] = max(PolyLine_slope);
        end
    end
    C = PolyfitValues(3);
    t = InspTime(SlopeIndexOffset)-InspTime(1);
    F = InspFlow(SlopeIndexOffset);
    MansourRUA(i)= abs(C /((F^2)*t)); 
    if ShowFigures
        figure(59)
        subplot(2,1,1); hold on;
        plot(InspTime(SlopeIndexOffset), PolyLine(SlopeIndexOffset), 'rx');
        xlim([InspTime(range(1)),InspTime(range(end))]);
        subplot(2,1,2); hold on;
        plot(InspTime(SlopeIndexOffset),SlopeVal,'rx'); 
        xlim([InspTime(range(1)),InspTime(range(end))]);
    end
    end 
    
    %% Wiriyaporn exponential decay fit to end expiratory flow
    if 0 % Expiratory decay constant 
        % currently switched off, as 'fit' function below is fussy...
        %fit decay model from peak exp flow to end expiration
        % using 6 co-ords at 0.2 sec intervals after PEF, before end exp
        FlowPts=ExpFlow(TPEF_i:end);
        TimePts=ExpTime(TPEF_i:end); 
        if numel(FlowPts) > 20
            
            % needs a better fit function, maybe try lsqcurvefit
            ft = fittype('exp1');
            [f3, gof3, ~] = fit_new(TimePts', FlowPts', ft); %
            
            cftool(TimePts', FlowPts');
        
            if gof3.rsquare > 0.8 % was 0.95
                Wiriyaporn(i) = -1/(f3.b);
                if ShowFigures
                    figure(60); clf(figure(60));
                    plot(TimePts, FlowPts,'k'); hold on;
                    fitdata = f3(TimePts);
                    plot(TimePts, fitdata,'r');
                    legend('Flow','Fit');
                    title('Exponential decay from PEF to end of expiration');
                    hold off;
                end
            end
        else
            % Wiriyaporn(i) = 0; % not enough points to calculate decay
        end
    end
    
    %% >>>>>>>>>>>>>>>>>>>>  I:E MEASURES   <<<<<<<<<<<<<<<<<<<<<<<<
    Ti_Te(i) = Ti/Te;
    PIF_PEF(i)=PIF/PEF; %PIF should be one
    MIF_MEF(i)=MIF/MEF;
    MIF50_MEF50(i) = MIF50(i) / MEF50(i);
    
    %% SeriesIEflow50, ratio of insp flow to exp flow at mid volumes, (<0.97)
    % mid-tidal volume-flow ratio
    halfinspvol = 0.5*VTi; % find 1/2 insp vol
    Inspind = find(InspVol>halfinspvol,1,'first'); % find when 1/2 insp vol occurs in time
    inspflowatmidTV = InspFlow(Inspind); % find flow at this time
    % repeat for expiratory
    halfexpvol = -0.5*VTe;
    Expind = find(ExpVol<halfexpvol,1,'first');
    expflowatmidTV = abs(ExpFlow(Expind));
    if ~isempty(Inspind) && ~isempty(Expind)
        SeriesIEflow(i) = inspflowatmidTV/expflowatmidTV;
    end
    
    %% modifications on Series
    % ratio of insp time to exp time at mid volumes
    SeriesIEtime(i) = (InspTime(Inspind)-InspTime(1)) / (ExpTime(Expind)-ExpTime(1));
    
    % ratio of insp vol to exp volume at mid time (this would be same as Kaplan...)
    %SeriesIEvol(i) = InspVol(round(length(InspTime)/2)) / ExpVol(round(length(ExpTime)/2));
    
    if ShowFigures
        figure(61); clf(figure(61));
        subplot(4,1,1);
        plot(InspTime, InspFlow, 'k'); hold on;
        plot(InspTime, InspVol, 'r');
        plot(InspTime(Inspind), InspVol(Inspind), 'ro');
        subplot(4,1,2);
        plot(ExpTime, ExpFlow, 'k'); hold on;
        plot(ExpTime, ExpVol, 'r');
        plot(ExpTime(Expind), ExpVol(Expind), 'ro');
        subplot(4,1,[3:4]);
        %[val1,ind1] = min(ExpFlow);
        %x=[expvol(ind1),expvol(end)];
        %y=[val1,ExpFlow(end)];
        plot(InspVol,InspFlow);
        hold on;
        plot(fliplr(-ExpVol),ExpFlow);
        %plot(x,y);
        ylabel('Flow (-exp, +insp)')
        xlabel('Volume')
        axis('tight')
    end
    
    %% KaplanIEvol50, ratio of insp vol to exp vol at mid times
    MidInspTime = round((BB_I_end-BB_I_start)/2);
    MidExpTime = round((BB_E_end-BB_E_start)/2);
    KaplanIEvol(i) = InspVol(MidInspTime) / ExpVol(MidExpTime);

    %% >>>>>>>>>>>>>>>>>>>>  OTHERS   <<<<<<<<<<<<<<<<<<<<<<<<
    %TAA(i) = nanmedian(TAA_(BB_I_start:BB_E_end));

    % Ttran values
    TTran_i_Ti(i) = TiTrans(i)/Ti;
    TTran_i_Ttot(i) = TiTrans(i)/Ttot;
    TTran_e_Te(i) = TeTrans(i)/Te;
    TTran_e_Ttot(i) = TeTrans(i)/Ttot;
    
    % morgensterns measures
    MorgensVPIFVTi(i) = InspVol(TPIF_i) / VTi;  % volume at time of PIF, normalized by VTi
    %MorgensV05(i) = VTiMod(Ti50)/VTi;
    MorgensV03(i) = InspVol(Ti30) / VTi;        % volume at 1/3 Ti, normalized by VTi
    
    if ShowFigures
        figure(53); hold on; subplot(2,1,1);
        %plot(InspTime, InspFlow, 'k'); hold on;
        %plot(InspTime, InspVol, 'r');
        plot(InspTime(TPIF_i), InspVol(TPIF_i), 'b^');
        plot(InspTime(Ti50), InspVol(Ti50), 'bo');
        plot(InspTime(Ti30), InspVol(Ti30), 'bo');
        %title('Inspiratory flow and volume');
        %plot(InspTime(Ti25:Ti75), 0.5*ones(length(InspTime(Ti25:Ti75)),1), 'm--');
        %plot(InspTime(Ti25:Ti75), MorgensVPIFVTi(i)*ones(length(InspTime(Ti25:Ti75)),1), 'r--');
        % 0.5 is "normal". away in either direction is abnormal. WHY ?
        % should analyse feature as abs(delta(0.5, breath_value))
    end
    
    if PublicationFigure
        figure(55); hold on; subplot(4,1,4);
        plot(InspTime(Ti30), MorgensV03(i), 'rs');
    end
    
    
    %% >>>>>>>>>>>>>>>>  RISE, DWELL AND FALL TIMES  <<<<<<<<<<<<<<<<<<<<<<
    % Determine insp and exp rise times, fall times, dwell times.
    PIF10i = find(InspFlow>(0.1*PIF), 1, 'first');
    PIF90i = find(InspFlow>(0.9*PIF), 1, 'first');
    PIF90ii = find(InspFlow>(0.9*PIF), 1, 'last');
    PIF10ii = find(InspFlow>(0.1*PIF), 1, 'last');
    if ~isempty(PIF10i) && ~isempty(PIF90i)
        RTi(i) = ((PIF90i-PIF10i)*dt) / Ti;  % insp rise time
    end
    if ~isempty(PIF10ii) && ~isempty(PIF90ii)
        FTi(i) = ((PIF10ii-PIF90ii)*dt) / Ti; % insp fall time
    end
    if ~isempty(PIF90ii) && ~isempty(PIF90i)
        DTi(i) = ((PIF90ii-PIF90i)*dt) / Ti; % insp dwell time
    end
    
    PEF10i = find(ExpFlow<(0.1*PEF), 1, 'first');
    PEF90i = find(ExpFlow<(0.9*PEF), 1, 'first');
    PEF90ii = find(ExpFlow<(0.9*PEF), 1, 'last');
    PEF10ii = find(ExpFlow<(0.1*PEF), 1, 'last');
    if ~isempty(PEF10i) && ~isempty(PEF90i)
        RTe(i) = ((PEF90i-PEF10i)*dt) / Te;  % exp rise time
    end
    if ~isempty(PEF10ii) && ~isempty(PEF90ii)
        FTe(i) = ((PEF10ii-PEF90ii)*dt) / Te; % exp fall time
    end
    if ~isempty(PEF90ii) && ~isempty(PEF90i)
        DTe(i) = ((PEF90ii-PEF90i)*dt) / Te; % exp dwell time
    end
    
    if ShowFigures
        figure(53); hold on; 
        subplot(2,1,1); hold on; legend off;
        plot(InspTime(PIF10i), InspFlow(PIF10i), 'rx');
        plot(InspTime(PIF90i), InspFlow(PIF90i), 'rx');
        plot(InspTime(PIF90ii), InspFlow(PIF90ii), 'rx');
        plot(InspTime(PIF10ii), InspFlow(PIF10ii), 'rx');
        subplot(2,1,2); hold on;
        plot(ExpTime(PEF10i), ExpFlow(PEF10i), 'rx');
        plot(ExpTime(PEF90i), ExpFlow(PEF90i), 'rx');
        plot(ExpTime(PEF90ii), ExpFlow(PEF90ii), 'rx');
        plot(ExpTime(PEF10ii), ExpFlow(PEF10ii), 'rx');
    end
    
    %% >>>>>>>>>>>>>>>>  cf REFERENCE BREATH  <<<<<<<<<<<<<<<<<<<<<<
    % Mooney
    if 0 % removed by DLM, uses reference breath
        MooneyTi(i) = Ti;
    end
    % Morris, TPEF/Te
    if 0    % removed by DLM, same as TpeakE_Te(i)
        MorrisTPEF_Te(i) = TPEF/Te;
    end
    
    % ATS, VPEF/VTe (expiratory equivalent to morgens insp Vol at PIF)  
    ATS_VPEF_VTe(i) = abs(ExpVol(TPEF_i) / VTe);
   
    if ShowFigures
        figure(53); hold on; subplot(2,1,2);
        %plot(ExpTime, ExpFlow,'k'); hold on;
        %plot(ExpTime, ExpVol,'r');
        plot(ExpTime(TPEF_i), ExpVol(TPEF_i), 'bv');
        %title('Expiratory flow and volume');
        % plot(ExpTime(Te25:Te75), -1*(ATS_VPEF_VTe(i)*ones(length(ExpTime(Te25:Te75)),1)), 'r--');
        % plot(ExpTime(Te25:Te75), -0.5*ones(length(ExpTime(Te25:Te75)),1), 'm--');
        % -0.5 is "normal". away in either direction is abnormal. WHY?
        % should analyse feature as abs(delta(0.5, breath_value))
    end
    
    % Clark template - can not use as requires reference breath
    
    %% >>>>>>>>>>>>>>>>>>>>>>  FLUTTER   <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if IncludeFlutterFtrs   
        FlowToUse = 'SetLengthFlow';% 'common'; % 'original'; % 'DSwithFlutter'; %
        % 'common' is the flow that is used everywhere else in this code,
        % which is downsampled, and normalised, with optional de-fluttering
        % 'original' is the flow passed into this function, no modification
        % 'DSwithFlutter' is downsampled, but not flutter removed.
        % note that if RemoveFlutter is switched off, then 'DSwithFlutter'
        % will be a non-normalised version of the 'common' signal.
        % 'SetLengthFlow' is flow forced to 200 samples for insp, 250 exp.
        % 
        switch FlowToUse
            case 'original' % does not include baseline adjustment
                InspFlow_forFlutter = InspFlow_original;
                ExpFlow_forFlutter = ExpFlow_original;
            case 'DSwithFlutter' % does not include baseline adjustment
                InspFlow_forFlutter = InspFlow_DSwithFlutter;
                ExpFlow_forFlutter = ExpFlow_DSwithFlutter;
            case 'SetLengthFlow' %
                InspFlow_forFlutter = InspFlow_SetLength;
                ExpFlow_forFlutter = ExpFlow_SetLength;
            otherwise
                InspFlow_forFlutter = InspFlow;
                ExpFlow_forFlutter = ExpFlow;
        end

        %% the original method
        [Pxx1,f1] = periodogram(InspFlow_forFlutter,[],[],Fs);
        [Pxx2,f2] = periodogram(ExpFlow_forFlutter,[],[],Fs);
        % [Pxx1w,f1]=pwelch(InspFlow_forFlutter,[],[],[],Fs);
        % [Pxx2w,f2]=pwelch(ExpFlow_forFlutter,[],[],[],Fs);
        Pow1_Orig=bandpower(Pxx1,f1,[5,round(Fs/2)-1],'psd');
        Pow2_Orig=bandpower(Pxx2,f2,[5,round(Fs/2)-1],'psd');
        InspFlutPowOrig(i)=Pow1_Orig/PIF^2;
        ExpFlutPowOrig(i)=Pow2_Orig/PEF^2;
        InspExpFlutPowOrig(i)=Pow1_Orig/Pow2_Orig;
        InspExpFlutPowOrig_Sum(i)=Pow1_Orig/(Pow1_Orig+Pow2_Orig);
        %InspFlutPow_VT2(i)=Pow1/VT^2; % !! units of 1/ times squared
        %ExpFlutPow_VT2(i)=Pow2/VT^2;

        %% Scotty's alternate (and usually quicker) method    
        Frange_low = [4 7];
        Frange_high = [8 12];
        
        % inspiration
        N1=length(InspFlow_forFlutter);
        % w1 = hann(N1); w1 = w1/rms(w1);
        % X1 = fft((InspFlow_original-mean(InspFlow_original)).*w)/N*2;
        T01 = N1*dt; % DLM. query use of set length, dt could be different...
        df1 = 1/T01;
        X1 = fft(InspFlow_forFlutter)/N1*2;
        Pxx1 = conj(X1).*X1; % equivalent to Pxx1 = abs(X1.^2);
        
        if 0
            Pxx1 = conj(X1).*X1;
            Pxx1_ = abs(X1);
            Pxx1__ = abs(X1.^2);
            ax(1)=subplot(3,1,1); plot(Pxx1);
            ax(2)=subplot(3,1,2); plot(Pxx1_);
            ax(3)=subplot(3,1,3); plot(Pxx1__);
            linkaxes(ax);
        end
               
        % 4 to 7
        Frangei1 = round(Frange_low/df1); %frequency bins of interest based on frequency range
        Pow1 = sum(Pxx1(Frangei1(1):Frangei1(2)))*df1; %Pow1_ = sum(Pxx1(Frange))*df1;
        % 8 to 12
        Frangei3 = round(Frange_high/df1); %frequency bins of interest based on frequency range
        Pow3 = sum(Pxx1(Frangei3(1):Frangei3(2)))*df1; %Pow3_ = sum(Pxx1(Frange3))*df1;

        % expiration
        N2=length(ExpFlow_forFlutter);
        % w2 = hann(N2); w2 = w2/rms(w2);
        % X2 = fft((ExpFlow_original-mean(ExpFlow_original)).*w)/N*2;
        T02 = N2*dt; df2 = 1/T02;
        X2 = fft(ExpFlow_forFlutter)/N2*2;
        Pxx2 = conj(X2).*X2;
        % 4 to 7
        Frangei2 = round(Frange_low/df2); %frequency bins of interest based on frequency range
        Pow2 = sum(Pxx2(Frangei2(1):Frangei2(2)))*df2;
        % 8 to 12
        Frangei4 = round(Frange_high/df2); %frequency bins of interest based on frequency range
        Pow4 = sum(Pxx2(Frangei4(1):Frangei4(2)))*df2;
        
        InspFlutPow4to7(i)=Pow1/PIF^2; % recall that PIF = 1.
        ExpFlutPow4to7(i)=Pow2/PEF^2;
        InspExpFlutPow4to7(i)=Pow1/Pow2;
        InspExpFlutPow4to7_Sum(i)=Pow1/(Pow1+Pow2);
        InspFlutPow8to12(i)=Pow3/PIF^2;
        ExpFlutPow8to12(i)=Pow4/PEF^2;
        InspExpFlutPow8to12(i)=Pow3/Pow4;
        InspExpFlutPow8to12_Sum(i)=Pow3/(Pow3+Pow4);
        
        %% dlm other alternate method (faster than built-in, 0.002sec slower than Scotty's)
        % same ranges as per Scotty's method above
        if 0
        x = InspFlow_forFlutter; x = x - nanmean(x);
        x(isnan(x))=0; % set zero for NaNs
        N=length(x); 
        nfft = 2^nextpow2(N); % next larger power of 2
        y = fft(x,nfft); % Fast Fourier Transform
        y = abs(y.^2); % raw power spectrum density
        %y = abs(y);
        y = y(1:nfft/2); % half-spectrum, plot(y)
        f_scale = (1:nfft/2)* Fs/nfft; % frequency scale,  plot(f_scale, y)
        Pow1_dm = sum(y(f_scale>Frange_low(1)&f_scale<Frange_low(2)))*dt; % dt is set on original breath timing/duration/spacing
        Pow3_dm = sum(y(f_scale>Frange_high(1)&f_scale<Frange_high(2)))*dt;
        
        x = ExpFlow_forFlutter; x = x - nanmean(x);
        x(isnan(x))=0; % set zero for NaNs
        N=length(x); nfft = 2^nextpow2(N); % next larger power of 2
        y = fft(x,nfft); % Fast Fourier Transform
        y = abs(y.^2); % raw power spectrum density
        y = y(1:1+nfft/2); % half-spectrum
        f_scale = (0:nfft/2)* Fs/nfft; % frequency scale
        Pow2_dm = sum(y(f_scale>Frange_low(1)&f_scale<Frange_low(2)))*dt;
        Pow4_dm = sum(y(f_scale>Frange_high(1)&f_scale<Frange_high(2)))*dt;            
        
        InspFlutPow4to7DLM(i)=Pow1_dm/PIF^2; % recall that PIF = 1.
        ExpFlutPow4to7DLM(i)=Pow2_dm/PEF^2;
        InspExpFlutPow4to7DLM(i)=Pow1_dm/Pow2_dm;
        InspFlutPow8to12DLM(i)=Pow3_dm/PIF^2;
        ExpFlutPow8to12DLM(i)=Pow4_dm/PEF^2;
        InspExpFlutPow8to12DLM(i)=Pow3_dm/Pow4_dm;
        end
    end
    
    %% >>>>>>>>>>>>>>>>>>>>  Ali Peaks   <<<<<<<<<<<<<<<<<<<<<<<<
    % This set of functions has been taken from Ali's code, which
    % identifies peaks in the inspiratory flow, and then performs a
    % range of analysis upon those peaks.
    [pks,locs,w,p]=findpeaks(InspFlow,'MinPeakProminence',0.1,'Annotate','extents');
    if ShowFigures
        figure(53); hold on;
        subplot(2,1,1); hold on;
        legend('off');
        plot(InspTime(locs),pks,'gd');
    end
    if ~isempty(locs)
        if locs(1)<0.5*BB_I_end
            I_01=locs(1);
        else
            I_01=round(0.5*BB_I_end);
        end
        if locs(end)>0.8*BB_I_end
            I_09=locs(end);
        else
            I_09=round(0.8*BB_I_end);
        end
    else
        I_01=Ti10;
        I_09=Ti90;
    end
    if ~isempty(p)
        % modified by DLM to handle multiple values returned by max
        pp = find(p==max(p), 1, 'first');
        WidthMostProminentPeak=w(pp)/Fs/Ti; % normalized to Ti
        PromMostProminentPeak=p(pp);
    else
        WidthMostProminentPeak=NaN;
        PromMostProminentPeak=NaN;
    end
    
    try
    if length(locs)>=3
        % modified by DLM to handle more than three peaks
        % if there are more than three peaks, then need to assign which one
        % is to be peak 1, which is to be peak 2 and which is to be peak 3.
        % Do this by finding the three biggest p (prominence) values, and
        % using them in the order they appear through the breath
        [~, a] = sort(p, 'descend');
        [a, ~] = sort(a(1:3), 'ascend');
        fp1_idx=locs(a(1))/Fs;
        fp2_idx=locs(a(2))/Fs;
        fp3_idx=locs(a(3))/Fs;
        fp1=pks(a(1)); fp2=pks(a(2));  fp3=pks(a(3));    
        pkProm1=p(a(1))/PromMostProminentPeak;
        pkProm2=p(a(2))/PromMostProminentPeak;
        pkProm3=p(a(3))/PromMostProminentPeak;       
    elseif length(locs)==2
        fp1_idx=locs(1)/Fs; fp2_idx=locs(2)/Fs; fp3_idx=NaN;
        fp1=pks(1); fp2=pks(2);  fp3=NaN;      
        pkProm1=p(1)/PromMostProminentPeak;
        pkProm2=p(2)/PromMostProminentPeak;
        pkProm3=NaN;    
    elseif length(locs)==1
        fp1_idx=locs(1)/Fs; fp2_idx=NaN; fp3_idx=NaN;
        fp1=pks(1); fp2=NaN; fp3=NaN;    
        pkProm1=p(1)/PromMostProminentPeak;
        pkProm2=NaN;  pkProm3=NaN;    
    else
        fp1_idx=NaN; fp2_idx=NaN; fp3_idx=NaN;
        fp1=NaN; fp2=NaN; fp3=NaN;     
        pkProm1=NaN; pkProm2=NaN; pkProm3=NaN;
    end 
    catch me
        disp(me.message); me.getReport
    end
    
    try    
    %Vpeak1_Vpeak(i)=fp1/PIF;
    %Vpeak2_Vpeak(i)=fp2/PIF;
    %Vpeak3_Vpeak(i)=fp3/PIF;
    TpeakI_Ti(i)=TPIF/Ti;
    %Tpeak1_Ti(i)=fp1_idx/Ti;
    %Tpeak2_Ti(i)=fp2_idx/Ti;
    %Tpeak3_Ti(i)=fp3_idx/Ti;
    TpeakE_Te(i)=TPEF/Te;               % equivalent to MorrisTPEF_Te
    TpeakI_TpeakE(i)=TPIF/TPEF;
    %MostPromPeakW(i)=WidthMostProminentPeak;
    MPPW_Ti(i)=WidthMostProminentPeak; % already normalised to Ti, don't do again here.
    %PkProm1_n(i)=pkProm1;
    %PkProm2_n(i)=pkProm2;
    %PkProm3_n(i)=pkProm3;
    catch me
        disp(me.message); me.getReport
    end
    if ShowFigures
        figure(53); hold on;
        subplot(2,1,1); hold on;
        if ~isnan(fp1_idx)
        plot(InspTime(locs(1)),pks(1),'g^');
        end
        if ~isnan(fp2_idx)
        plot(InspTime(locs(2)),pks(2),'g^');
        end
        if ~isnan(fp3_idx)
        plot(InspTime(locs(3)),pks(3),'g^');
        end  
    end
    
    %% >>>>>>>>>>>>>>>>>>>>  ALI NED   <<<<<<<<<<<<<<<<<<<<<<<<
    % these NED methods require the respective peaks as determined above
    % this one is excluded because it really requires the peak to be in the
    % first third of the breath, and the 'min' should be an average of the
    % second two thirds. The Span is a good way to go, but sometimes it is
    % not possible after the TPIF_i.
    if 0
    NEDSpan=round(0.25*(BB_I_end-1));
    if ~isnan(pks)
        NEDMinValIdx=TPIF_i:round(BB_I_end/2+NEDSpan);
        if ~isempty(NEDMinValIdx)
            Ali_NED(i)=(PIF)-min(InspFlow(NEDMinValIdx))/PIF;
        end
    end  
    end
    
    %% >>>>>>>>>>>>>>>>>>>>  AA Peaks   <<<<<<<<<<<<<<<<<<<<<<<<
    % This set of functions has been taken from Alessandra and Angelos code
    % where peaks are identified in the inspiratory flow, however, these
    % peaks are distributed across tertiles of the inspiratory flow. 
    % After detecting up to three possible peaks, some analysis then occurs
    [AA_NumOfPeaks(i),fp1, fp2, fp3, AA_IsTerminalPeak(i),Flow_peaks, Flow_troughs]= ...
        AApeaks(InspFlow,BB_I_end,PIF);
    
    %AA_Pf_1_3(i)=fp1;
    %AA_Pf_2_3(i)=fp2;
    %AA_Pf_3_3(i)=fp3;
    
    % AA PeaksRatio
    if 0 %removed on 20171005, SS. Highly non-monotonic relationship to FL
        P=[fp1 fp2 fp3];
        P(isnan(P))=0; % replace NaN's with zeros
        AA_PeaksRatio(i)=sum(P./max(P))/3;
    end
    
    if ShowFigures 
        figure(53); hold on;
        subplot(2,1,1); hold on;
        plot(InspTime(Flow_peaks(:,1)),Flow_peaks(:,2),'ks');
        plot(InspTime(Flow_troughs(:,1)),Flow_troughs(:,2),'g.');
    end
    
    %% >>>>>>>>>>>>>>>>>>>>  AA NED   <<<<<<<<<<<<<<<<<<<<<<<<
%     % need to confirm use of normalised flow signal
%     if ~isnan(fp1)
%         AA_NED(i)=100*(fp1-InspFlow(round(BB_I_end/2)))/fp1;
%     elseif ~isnan(fp2)
%         AA_NED(i)=100*(fp2-InspFlow(round(BB_I_end/2)))/fp2;
%     elseif ~isnan(fp3)
%         AA_NED(i)=100*(fp3-InspFlow(round(BB_I_end/2)))/fp3;
%     end
    % reworked by DLM - what proportion of PIF is flow at mid inspiration
    % PIF is always 1, so this is 1 - Flow@mid-inspiration
    AA_NED(i) = (PIF - InspFlow(Ti50) / PIF); 
    if AA_NED(i)<0; AA_NED(i)=0; end  
    if AA_NED(i)>1; AA_NED(i)=1; end
    
    %% >>>>>>>>>>>>>>  AA Expiratory Flow Limitation   <<<<<<<<<<<<<<<<<<<
    % Find Expiratory Flow Limitation Index by AA method
    if 0
        AA_EFLI(i)=IsExpFlowLimited(flow,BB_I_start,BB_E_end,Fs);
    end
    
    %% >>>>>>>>>>>>>>>>>>>>  DISCONTINUITY   <<<<<<<<<<<<<<<<<<<<<<<<
    
    %% Ali
    % these use SetLength flow
    if 0 % option to not run "discontinuity", if process taking too long
        % Calculating smooth derivative
        if 1  % using SetLength InspFlow
            FlowI1_I2=InspFlow_SetLength(Ti10_set:Ti90_set);
            TimeI1_I2 = 0:dt:length(FlowI1_I2)/Fs; % make new "Time"
        else % using normal InspFlow
            FlowI1_I2=InspFlow(I_01:I_09);
            TimeI1_I2=InspTime(I_01:I_09)';
        end
        
        %tt=-(I_01_chnge):length(Flow)-(I_01_chnge)-1;
        ipt=findchangepts(FlowI1_I2,'MinThreshold' ,0.005,'Statistic','linear');
        
        if 0
            figure(1); clf(figure(1));
            %plot(FlowI1_I2);
            plot(TimeI1_I2, FlowI1_I2); hold on;
            plot(TimeI1_I2(ipt), FlowI1_I2(ipt), 'ro');         
        end
        
        if ~isempty(ipt)
            ipt2=[1,ipt,length(FlowI1_I2)]; % add start and end points
            slp = nan(length(ipt2)-1,1);
            stFlow = nan(length(ipt2)-1,1);
            endFlow = nan(length(ipt2)-1,1);
            for ii=2:length(ipt2)
                X = [ones(size(TimeI1_I2(ipt2(ii-1):ipt2(ii)-1),2),1), TimeI1_I2(ipt2(ii-1):ipt2(ii)-1)'];
                % b = regress(y,X)
                % b is a vector of coefficient estimates for a multilinear 
                % regression of the responses in y on the predictors in X
                warning('off','stats:regress:RankDefDesignMat');
                b = regress(FlowI1_I2(ipt2(ii-1):ipt2(ii)-1)',X); % ignores nans, and removes them
                slp(ii-1)=b(2);
                intercept=b(1);
                stFlow(ii-1)=intercept+TimeI1_I2(ipt2(ii-1))*slp(ii-1);
                endFlow(ii-1)=intercept+TimeI1_I2(ipt2(ii)-1)*slp(ii-1);
            end 
            rise=endFlow-stFlow;
            Disc=[0; stFlow(2:end)-endFlow(1:end-1)]; 
            riseSlp=abs(rise.*slp);
            rise=abs(rise);
            Disc=abs(Disc);
            FlowChangeSrtd=NaN(5,1);
            riseSrtd=NaN(5,1);
            DiscSrtd=NaN(5,1);
            RiseSlope=NaN(5,1);
            
            if ~isempty(slp)
                [FlowChangeSrtd_temp,srtd_idx]=sort(abs(slp),'descend');
                rise_tmp=rise(srtd_idx);
                Disc_tmp=Disc(srtd_idx);
                riseSlp_tmp=sort(riseSlp,'descend');
                if length(FlowChangeSrtd_temp)<5
                    FlowChangeSrtd_temp(length(FlowChangeSrtd_temp):5)=min(FlowChangeSrtd_temp);
                    rise_tmp(length(rise_tmp):5)=rise_tmp(end);
                    Disc_tmp(length(Disc_tmp):5)=Disc_tmp(end);
                    riseSlp_tmp(length(riseSlp_tmp):5)=riseSlp_tmp(end);
                end
                
                for ii=1:length(FlowChangeSrtd_temp)
                    if ii<=5
                        FlowChangeSrtd(ii)=FlowChangeSrtd_temp(ii);
                        riseSrtd(ii)=rise_tmp(ii);
                        DiscSrtd(ii)=Disc_tmp(ii);
                        RiseSlope(ii)=riseSlp_tmp(ii);
                    else
                        break;
                    end
                end
                if 1 % using SetLength InspFlow, do not normalise again
                    D1(i)=FlowChangeSrtd(1);
                    D2(i)=riseSrtd(1);
                    D3(i)=DiscSrtd(1);
                    D4(i)=RiseSlope(1);
                else 
                    D1(i)=FlowChangeSrtd(1)/VT;
                    D2(i)=riseSrtd(1)/VT;
                    D3(i)=DiscSrtd(1)/VT;
                    D4(i)=RiseSlope(1)/VT^2;
                end
            end
        end
    end
    
    %% AA
    if 0 % option to not run "rate of change", if process taking too long
        % Calculating smooth rate of change
        if 1 % using SetLength InspFlow
            [~,~,~,FlowDiff]=Vbox(InspFlow_SetLength);
            flowDiffSmooth=FlowDiff(Ti05_set:Ti95_set);
        else % using normal InspFlow
            [~,~,~,FlowDiff]=Vbox(InspFlow);
            flowDiffSmooth=FlowDiff(Ti05:Ti95);
        end
        
        %figure(); plot(Time(Ti05:Ti95), flowDiffSmooth); hold on;
        %plot(Time(Ti05:Ti95), InspFlow(Ti05:Ti95));
        
        MinFlowChange(i)=nanmin(flowDiffSmooth);
        MedianFlowChange(i)=nanmedian(flowDiffSmooth);
        STDFlowChange(i)=nanstd(flowDiffSmooth); 
    end
    
    
    %% >>>>>>>>>>>>>>>>>>>>>>>>  SS AREA   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % In v4 and before, this used the peaks and troughs returned from the
    % AA peaks method above. 
    % From v5 onwards, this uses it's own peakdet function, that finds 
    % peaks located in each tertile of the inspiratory flow.
    % SS Area finds the area between the Insp flow signal and a chord drawn
    % between the peaks.

    [ThePeaks] = FindPeaksInInspTertiles(InspFlow);
    
    %[ThePeaks_] = FindPeaksInInspTertiles(flow(1:116));
    
%     % for debugging, compare Flow_peaks and ThePeaks
%     if ~isequal(ThreePeaks, Flow_peaks)
%         disp('Not equal peaks');
%         Flow_peaks_bckup = Flow_peaks;
%         Flow_peaks = ThreePeaks;
%     end

    Flow_peaks = ThePeaks; % overwrite Flow_peaks, not used hereafter (but doing this will break the parfor loop)
    
    % could also use coefficients = polyfit([x(1), x(2)], [y(1), y(2)], 1);
    % to get a = coefficients (1); and b = coefficients (2);, but
    % polyfit is much slower than mldivide (\) used below
    switch size(Flow_peaks,1) %AA_NumOfPeaks(i)
        case 1 % if AA_NumOfPeaks(i) = 1, return zero
            SS_Area_ = 0; % set to zero, 
        case 2 % if AA_NumOfPeaks(i) = 2, simple area
            % find the area between flow signal, and a chord between peaks
            x = [InspTime(Flow_peaks(1,1)) InspTime(Flow_peaks(2,1))];
            y = [Flow_peaks(1,2) Flow_peaks(2,2)];
            c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
            % alternative method
            %coefficients = polyfit([x(1), x(2)], [y(1), y(2)], 1);
            %m = coefficients(1); c = coefficients(2);
            %X = [x(1),y(1);x(2),y(2)];
            %ln = round(pdist(X,'euclidean'))
            Xrange = InspTime(Flow_peaks(1,1):Flow_peaks(2,1));
            ln = length(Xrange);        
            chordA = NaN(ln,1);
            for k=1:ln
                chordA(k) = (c(2)*(Xrange(k)))+c(1);
            end
            % find abs area between the chord (between peaks) and flow
            %SS_Area(i)=abs(trapz(Xrange,InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))) ...
            %    - trapz(Xrange,chordA));
            diff_line = chordA' - InspFlow(Flow_peaks(1,1):Flow_peaks(2,1));
            diff_line = max(diff_line,0); % only want the positive values
            SS_Area_ =trapz(Xrange,diff_line); % did have /Ti here, but it does this at the end
            
            if ShowFigures
                %figure(53); hold on; legend off, subplot(2,1,1); hold on;
                figure(63); plot(InspTime, InspFlow); hold on;
                plot(Xrange, chordA, 'r--');
                plot(Xrange, diff_line,'r--');
            end
            
             if PublicationFigure
                figure(55); hold on; subplot(4,1,[1:2]);
                plot(Xrange, chordA, 'b--');
                
                subplot(4,1,3);
                plot(Xrange, diff_line,'b');
            end
            
        case 3 % if AA_NumOfPeaks(i) = 3, max ( A-B + B-C ; A-C )
            % AB chord (= between peaks 1 and 2)
            x = [InspTime(Flow_peaks(1,1)) InspTime(Flow_peaks(2,1))];
            y = [Flow_peaks(1,2) Flow_peaks(2,2)];
            c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
            XrangeAB = InspTime(Flow_peaks(1,1):Flow_peaks(2,1));
            ln = length(XrangeAB);        
            chordAB = NaN(ln,1);
            for k=1:ln
                chordAB(k) = (c(2)*(XrangeAB(k)))+c(1);
            end
            % AB area
            %Area_AB=abs(trapz(XrangeAB,InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))) ...
            %    - trapz(XrangeAB,chordAB)); 
            diff_lineAB = chordAB' - InspFlow(Flow_peaks(1,1):Flow_peaks(2,1));
            diff_lineAB = max(diff_lineAB,0); % only want the postive values
            Area_AB = trapz(XrangeAB,diff_lineAB);
            
            % BC chord (= between peaks 2 and 3)
            x = [InspTime(Flow_peaks(2,1)) InspTime(Flow_peaks(3,1))];
            y = [Flow_peaks(2,2) Flow_peaks(3,2)];
            c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
            XrangeBC = InspTime(Flow_peaks(2,1):Flow_peaks(3,1));
            ln = length(XrangeBC);        
            chordBC = NaN(ln,1);
            for k=1:ln
                chordBC(k) = (c(2)*(XrangeBC(k)))+c(1);
            end           
            % BC area
            %Area_BC=abs(trapz(XrangeBC,InspFlow(Flow_peaks(2,1):Flow_peaks(3,1))) ...
            %    - trapz(XrangeBC,chordBC));
            diff_lineBC = chordBC' - InspFlow(Flow_peaks(2,1):Flow_peaks(3,1));
            diff_lineBC = max(diff_lineBC,0); % only want the postive values
            Area_BC = trapz(XrangeBC,diff_lineBC);
            
            % AC chord (= between peaks 1 and 3)
            x = [InspTime(Flow_peaks(1,1)) InspTime(Flow_peaks(3,1))];
            y = [Flow_peaks(1,2) Flow_peaks(3,2)];
            c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
            XrangeAC = InspTime(Flow_peaks(1,1):Flow_peaks(3,1));
            ln = length(XrangeAC);        
            chordAC = NaN(ln,1);
            for k=1:ln
                chordAC(k) = (c(2)*(XrangeAC(k)))+c(1);
            end           
            % AC area
            %Area_AC=abs(trapz(XrangeAC,InspFlow(Flow_peaks(1,1):Flow_peaks(3,1))) ...
            %    - trapz(XrangeAC,chordAC));
            diff_lineAC = chordAC' - InspFlow(Flow_peaks(1,1):Flow_peaks(3,1));
            diff_lineAC = max(diff_lineAC,0); % only want the postive values
            Area_AC = trapz(XrangeAC,diff_lineAC);
            
            % use the greatest value
            [SS_Area_, indx] = max([(Area_AB + Area_BC) Area_AC]);
            
            if ShowFigures
                %figure(53); hold on; subplot(2,1,1); 
                figure(63); plot(InspTime, InspFlow); hold on;
                if indx == 1 % we used chords AB and BC
                    plot(XrangeAB, chordAB, 'r--');
                    plot(XrangeBC, chordBC, 'r--');
                    plot(XrangeAB,diff_lineAB,'r--');
                    plot(XrangeBC,diff_lineBC,'r--');
                else % we used chordAC
                    plot(XrangeAB, chordAB, 'g:');
                    plot(XrangeBC, chordBC, 'g:');
                    plot(XrangeAC, chordAC, 'r--');
                    plot(XrangeAC,diff_lineAC,'r--');
                end
            end 
            
            if PublicationFigure
                figure(55); hold on; 
                if indx == 1 % we used chords AB and BC
                    subplot(4,1,[1:2]); hold on;
                    plot(XrangeAB, chordAB, 'b--');
                    plot(XrangeBC, chordBC, 'b--');
                    subplot(4,1,3);
                    plot(XrangeAB, diff_lineAB,'b');
                    plot(XrangeBC, diff_lineBC,'b');
                else  % we used chordAC
                    subplot(4,1,[1:2]); hold on;
                     plot(XrangeAC, chordAC, 'b--');
                     subplot(4,1,3);
                     plot(XrangeAC,diff_lineAC,'b');
                end
            end
        otherwise % error case
            SS_Area_ = NaN; % set to NaN
    end
    SS_Area(i) = SS_Area_/Ti; 
    %pause(3);
    if PublicationFigure
        linkaxes(ax,'x');
        
        % do FFT image
        [fft_x, fft_y] = CheckSignalForNoisewFig(InspTime,InspFlow, Fs);
        figure(55); hold on;
        if 0
            axes('Position',[0.17 0.48 0.25 0.25]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]
        else
            axes('Position',[0.65 0.72 0.22 0.22]);
        end
        box on
        plot(fft_x, fft_y, 'k-'); yticks([]);
        xlim([-1,21]); xlabel('Hz'); ylabel('Amplitude');
        title('FFT');
        
        savefigs = 0;
        if savefigs
            str = ['..\Figures\Figure5\KeyFeatures_Pt2_BB_', num2str(i),'_wFFT_altPos'];
            saveas(fig, str, 'png'); %savefig(str);
        end %
    end
    
    if ShowFigures
        figure(53); fig = gcf; % bring this one to the front
        fig.Position = [ -6.5 -1 6 8]; % and move it somewhere
    end
end % Finished processing each Breath

warning('on', 'all');

%% combine into one big table
FT = FeatureNames';
BreathDataTable = table(...
    ... % M1,M2,M3,M4,M5,...  % used in testing only
    ... %% Flattening N=28
    MIF_PIF, Ali_InspFlat, Ali_ExpFlat, Ali_InspExpFlat, ...    
    AA_InspFlat9020, AA_InspFlat8020, AA_InspFlat7020, AA_InspFlat6020, AA_InspFlat5020, ...
    AA_InspFlat90Ti, AA_InspFlat80Ti, AA_InspFlat70Ti, AA_InspFlat60Ti, AA_InspFlat50Ti, ...
    AA_ExpFlat9020, AA_ExpFlat8020, AA_ExpFlat7020, AA_ExpFlat6020, AA_ExpFlat5020, ...
    AA_ExpFlat90Te, AA_ExpFlat80Te, AA_ExpFlat70Te, AA_ExpFlat60Te, AA_ExpFlat50Te, ...
    MIF50, MEF50, Teschler, MPPW_Ti, ...
    ... %% Scooping N=8
    AA_NED, SS_Area, SinI, SinI50, SinE, InvParabI, EllipseI, HypcosI, ...
    ... %% Asymmetry N=9
    AsymIndex, ... 
    SkewDistInsp, KurtDistInsp, SkewDistExp, KurtDistExp, ...
    AsymmetryInsp, AsymmetryExp, KurtDataInsp, KurtDataExp, ...
    ... %%  Ratios N=28
    Ti_Ttot, Ti_Te, TTran_i_Ti, TTran_i_Ttot, TTran_e_Te, TTran_e_Ttot, ...
    VTi_VTe, VTi_VT, VTe_VT, ATS_VPEF_VTe, MorgensVPIFVTi, MorgensV03, ...
    PIF_MIF, FTi, RTi, DTi, FTe, RTe, DTe, ...
    PIF_PEF, MIF_MEF, MIF50_MEF50, SeriesIEflow, SeriesIEtime, KaplanIEvol, ...
    TpeakI_Ti, TpeakE_Te, TpeakI_TpeakE, ...
    ... %% Flutter N=12
    InspFlutPowOrig, ExpFlutPowOrig, InspExpFlutPowOrig, InspExpFlutPowOrig_Sum, ...
    InspFlutPow4to7, ExpFlutPow4to7, InspExpFlutPow4to7, InspExpFlutPow4to7_Sum, ...
    InspFlutPow8to12, ExpFlutPow8to12, InspExpFlutPow8to12, InspExpFlutPow8to12_Sum, ...
    'VariableNames', FT);

end % end of MIFL function
