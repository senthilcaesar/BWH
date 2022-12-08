%% DeriveBreathData
% This code measures 60+ different things about each breath
%
% Input to this function:
%  Time
%  Flow
%  Clinical Scoring (Hypnog, arousals, resp events)
%  - and optionally
%  ReferenceBreathData
%
% Output from this function:
%  BreathData, each line is one breath, and each column is as follows: 
%
%ClinicalScoring
% 01 - SS - Sleep stage (0=wake, 1=N1 ... 3=N3, 5=REM, 8=undefined)
% 02 - AS - Arousal scoring (1 = wakefulness, 0 = sleep)
% 03 - CS - Central event scoring (1 = event, 0 = no event), incl Ap and Hypop
% 04 - OS - Obstructive event scoring (1 = event, 0 = no event), Ap and Hypop
%
%Times
% 05 - Time(BB_I_starti) - breath inspiratory start time
% 06 - Time(BB_I_endi) - breath inspiratory end time
% 07 - Time(BB_E_starti) - breath expiratory start time
% 08 - Time(BB_E_endi) - breath expiratory end time
% 09 - BB_I_starti - breath inspiratory start index
% 10 - BB_I_endi - breath inspiratory end index
% 11 - BB_E_starti - breath expiratory start index
% 12 - BB_E_endi - breath expiratory end index
% 13 - Ti - inspiratory time (seconds)
% 14 - Te - expiratory time (seconds)
% 15 - TiTran - inspiratory transition
% 16 - TeTran - expiratory transition
% 17 - Ttot - total breath time
% 18 - TPIF - time to peak inspiratory flow, from insp start
% 19 - TPEF - time to peak expiratory flow, from exp start
%
%Flows
% 20 - MIF - mean inspiratory flow
% 21 - MEF - mean expiratory flow
% 22 - PIF - peak inspiratory flow
% 23 - PEF - peak expiratory flow
% 24 - MIF50 - mean inspiratory flow, for middle 50% of breath by time
% 25 - MEF50 - mean expiratory flow, for middle 50% of breath by time
%
%Volumes
% 26 - VT - tidal volume
% 27 - VTi - inspiratory tidal volume
% 28 - VTe - expiratory tidal volume 
%       note, in this data, these three should essentially be the same
% 29 - Vdot - minute ventilation based on VT
% 30 - VI - VTi./Ttot, minute inspiratory volume
% 31 - VE - VTe./Ttot, minute expiratory volume
%       note, in this data, these three should essentially be the same
%
%BreathRates
% 32 - fR - Respiratory rate or frequency as 1./Ttot.
% 33 - FTi - inspiratory fall time, 90% to 10% (of peak flow)
% 34 - RTi - inspiratory rise time, 10% to 90% (of peak flow)
% 35 - DTi - inspiratory dwell time, flow >= 90% (of peak flow)
% 36 - FTe - expiratory fall time, 90% to 10% (of peak flow)
% 37 - RTe - expiratory rise time, 10% to 90% (of peak flow)
% 38 - DTe - expiratory dwell time, flow >= 90% (of peak flow)
%   Rise, dwell, and fall time values are also used in miller, 1998.
%
%BreathAreas, to standard curves
% 39 - SinI - Sine fit to inspiration
% 40 - SinE - Sine fit to expiration
%   an alternative to Clark area index, instead of using reference breath, 
%   could use difference to sinusoidal centroid shape, as above
% 41 - teschler
% 42 - invParabI
% 43 - ellipseI
% 44 - hypcosI
%
%BreathScores
% 45 - mansourFL, slope of polyfit at max flow (+ve = NIFL, -ve = IFL)
% 46 - mansourRUA, flow at min/max slope, and constant C from polynomial
% 47 - seriesIE50, ratio of insp flow to exp flow at mid volumes, (<0.97)
% 48 - morgensPIF, peak insp flow, normalized to pt peak flow
% 49 - morgensVPIFVTi, VTi at PIF / VTi, with flow norm'd to peak flow
% 50 - morgensV05, VTi at one half Ti, with flow normalized to peak flow
% 51 - morgensV03, VTi at one third Ti, with flow normalized to peak flow
% 52 - TAA, Thoracoabdominal asynchrony
%
%RefData
% 53 - clark, area index??
% 54 - schneider, Ti/Ttot
% 55 - mooney, Ti
% 56 - wiriyaporn, Expiratory decay constant,
% 57 - morris, TPEF/Te
% 58 - ats, VPEF/VTe
% 59 - kaplan, PTIF/MIF
%
% Potentially also outputs (when Reference data is included in input, i.e. varargin)
%Relative scores
% 60 - clarkRef, area index, difference btw ref breath and test breath as %age
%  of total area under test breath
% 61 - schneiderRef, Ti/Ttot relative to ref breath
% 62 - mooneyRef, Ti relative to ref breath
% 63 - wiriyapornRef, Expiratory decay constant, if bigger than ref = IFL
% 64 - morrisRef, TPEF/Te, relative to ref
% 65 - atsRef, VPEF/VTe, relative to ref
% 66 - kaplanRef, PTIF/MIF, relative to ref
%
%
%Possible expansions to the data configuration
%DAT{1,PT_Number}{1,ReferenceNumber} = "Breath Data", for the test breaths
%DAT{1,PT_Number}{2,ReferenceNumber} = "Breath Data", for Reference breaths
%
% but... do we just have one (possibly) averaged reference breath, or do we
% grab reference periods throughout the recording, and each block relates
% to it's coresponding reference period?
% and then do we use an average from a few breaths, as the reference data?

%% DeriveBreathData function
function [BreathData] = DeriveBreathData(time, flow, clinical_scoring, I, Ttot, VTi, VTe, varargin)
%% options and preprocessing

ShowFigures = 0; % set as 1 to enable figures, 0 for quiet operation
range = 0.95; % the range for AUC timing measurements

% figures used so far
% 9 : flow and timing marks
% 10 : whole breath with 10 and 90 markers
% 11 : teschler
% 12 : insp top panel, exp bottom panel
% 13 : insp only, zero centered
% 14 : flow top panel, derivative of polyfit bottom panel

%transpose if needed:
if size(flow,2)<size(flow,1)
    flow=flow';
end
%transpose if needed:
if size(time,2)<size(time,1)
    time=time';
end

%% Set up breath times - using previously found timing
BB_I_starti = I.starti;
BB_I_endi = I.midi;
BB_E_starti = I.midi;
BB_E_endi = I.endi;

if ShowFigures    
    figure(9)
    clf(figure(9));
    plot(time, flow); hold on
    plot(time(BB_I_starti), flow(BB_I_starti), 'go');
    plot(time(BB_I_endi), flow(BB_I_endi), 'g^');
    plot(time(BB_E_starti), flow(BB_E_starti), 'rv');
    plot(time(BB_E_endi), flow(BB_E_endi), 'r.');
    for i=1:length(BB_I_starti)
        str = num2str(i);
        text(time(BB_I_starti(i)), 0.4, str ,'FontSize',8, 'Color','b'); % print text onplot at loc
    end
    title('Breaths in current window');
end

%% Variables
dt=(time(end)-time(1))/(length(time)-1);

SS = nan(size(BB_I_starti,1),1);
AS = nan(size(BB_I_starti,1),1);
CS = nan(size(BB_I_starti,1),1);
OS = nan(size(BB_I_starti,1),1);

PIF = nan(size(BB_I_starti,1),1);
PEF = nan(size(BB_I_starti,1),1);
TPIF = nan(size(BB_I_starti,1),1);
TPEF = nan(size(BB_I_starti,1),1);
TPIFi = nan(size(BB_I_starti,1),1);
TPEFi = nan(size(BB_I_starti,1),1);
MIF = nan(size(BB_I_starti,1),1);
MEF = nan(size(BB_I_starti,1),1);
MIF50 = nan(size(BB_I_starti,1),1);
MEF50 = nan(size(BB_I_starti,1),1);

Ti25 = nan(size(BB_I_starti,1),1); % 1/4 insp time in breath, for entire series
TiI25 = nan(size(BB_I_starti,1),1); % 1/4 insp time in breath, for indiv breaths
Ti75 = nan(size(BB_I_starti,1),1);% 3/4 insp time in breath, for entire series
TiI75 = nan(size(BB_I_starti,1),1);% 3/4 insp time in breath, for indiv breaths
Te25 = nan(size(BB_I_starti,1),1);
Te75 = nan(size(BB_I_starti,1),1);

fR = nan(size(BB_I_starti,1),1);
FTi = nan(size(BB_I_starti,1),1);
RTi = nan(size(BB_I_starti,1),1);
DTi = nan(size(BB_I_starti,1),1);
FTe = nan(size(BB_I_starti,1),1);
RTe = nan(size(BB_I_starti,1),1);
DTe = nan(size(BB_I_starti,1),1);

SinI = nan(size(BB_I_starti,1),1); 
SinE = nan(size(BB_I_starti,1),1);
teschler = nan(size(BB_I_starti,1),1);
invParabI = nan(size(BB_I_starti,1),1);
ellipseI = nan(size(BB_I_starti,1),1);
hypcosI = nan(size(BB_I_starti,1),1);

mansourFL = nan(size(BB_I_starti,1),1);
mansourRUA = nan(size(BB_I_starti,1),1);
seriesIE50 = nan(size(BB_I_starti,1),1);

morgensPIF = nan(size(BB_I_starti,1),1);
morgensVPIFVTi = nan(size(BB_I_starti,1),1);
morgensV05 = nan(size(BB_I_starti,1),1);
morgensV03 = nan(size(BB_I_starti,1),1);

TAA = nan(size(BB_I_starti,1),1);

clark = nan(size(BB_I_starti,1),1);
schneider = nan(size(BB_I_starti,1),1);
mooney = nan(size(BB_I_starti,1),1);
wiriyaporn = nan(size(BB_I_starti,1),1);
morris = nan(size(BB_I_starti,1),1);
ats = nan(size(BB_I_starti,1),1);
kaplan = nan(size(BB_I_starti,1),1);

%% Step through each breath, and determine values for all features
for i = 1:size(BB_I_starti,1) 
    %% ClinicalScoring, assigning the predominant state/condition to each breath
    SS(i) = mode(clinical_scoring(BB_I_starti(i):BB_E_endi(i),1),1);
    AS(i) = mode(clinical_scoring(BB_I_starti(i):BB_E_endi(i),2),1);
    CS(i) = mode(clinical_scoring(BB_I_starti(i):BB_E_endi(i),3),1);
    OS(i) = mode(clinical_scoring(BB_I_starti(i):BB_E_endi(i),4),1);
    
    %% set up some useful variables to use throughout processing
    InspFlow = flow(BB_I_starti(i):BB_I_endi(i));
    InspTime = time(BB_I_starti(i):BB_I_endi(i));
    ExpFlow = flow(BB_E_starti(i):BB_E_endi(i));
    ExpTime = time(BB_E_starti(i):BB_E_endi(i));
    
    %% Flows and Volumes, not already derived above
    [PIF(i),TPIFi(i)] = max(InspFlow);
    TPIF(i) = TPIFi(i)*dt;   
    
    [PEF(i),TPEFi(i)] = min(ExpFlow);
    TPEF(i) = TPEFi(i)*dt;
    
    MIF(i) = mean(InspFlow); % mean or median
    MEF(i) = mean(ExpFlow); % mean or median
    
    Ti25(i) = round(BB_I_starti(i)+(0.25*(BB_I_endi(i)-BB_I_starti(i))),0);
    TiI25(i) = Ti25(i)-BB_I_starti(i);
    Ti75(i) = round(BB_I_starti(i)+(0.75*(BB_I_endi(i)-BB_I_starti(i))),0);
    TiI75(i) = Ti75(i)-BB_I_starti(i);
    Te25(i) = round(BB_E_starti(i)+(0.25*(BB_E_endi(i)-BB_E_starti(i))),0);
    Te75(i) = round(BB_E_starti(i)+(0.75*(BB_E_endi(i)-BB_E_starti(i))),0);
    
    MIF50(i) = mean(flow(Ti25(i):Ti75(i)));
    MEF50(i) = mean(flow(Te25(i):Te75(i)));
    
    %% Determine breath rates, ie. fR, dwell times etc.
    fR(i) = 1./Ttot(i); %frequency
    
    PIF10i = find(InspFlow>(0.1*PIF(i)), 1, 'first');
    PIF90i = find(InspFlow>(0.9*PIF(i)), 1, 'first');
    PIF90ii = find(InspFlow>(0.9*PIF(i)), 1, 'last');
    PIF10ii = find(InspFlow>(0.1*PIF(i)), 1, 'last');
    if ~isempty(PIF10i) && ~isempty(PIF90i)
    RTi(i) = (PIF90i-PIF10i)*dt;  % insp rise time
    end
    if ~isempty(PIF10ii) && ~isempty(PIF90ii)
    FTi(i) = (PIF10ii-PIF90ii)*dt; % insp fall time
    end
    if ~isempty(PIF90ii) && ~isempty(PIF90i)
    DTi(i) = (PIF90ii-PIF90i)*dt; % insp dwell time
    end
    
    PEF10i = find(ExpFlow<(0.1*PEF(i)), 1, 'first');
    PEF90i = find(ExpFlow<(0.9*PEF(i)), 1, 'first');
    PEF90ii = find(ExpFlow<(0.9*PEF(i)), 1, 'last');
    PEF10ii = find(ExpFlow<(0.1*PEF(i)), 1, 'last');
    if ~isempty(PEF10i) && ~isempty(PEF90i)
    RTe(i) = (PEF90i-PEF10i)*dt;  % exp rise time
    end
    if ~isempty(PEF10ii) && ~isempty(PEF90ii)
    FTe(i) = (PEF10ii-PEF90ii)*dt; % exp fall time
    end
    if ~isempty(PEF90ii) && ~isempty(PEF90i)
    DTe(i) = (PEF90ii-PEF90i)*dt; % exp dwell time
    end
    if ShowFigures
        figure(10);% show complete breath, with things found so far
        clf(figure(10));
        plot(time(BB_I_starti(i):BB_E_endi(i)), flow(BB_I_starti(i):BB_E_endi(i))); hold on;
        plot(InspTime(TPIFi(i)), PIF(i), 'r^');
        plot(InspTime(PIF10i), InspFlow(PIF10i), 'rx');
        plot(InspTime(PIF90i), InspFlow(PIF90i), 'rx');
        plot(InspTime(PIF90ii), InspFlow(PIF90ii), 'rx');
        plot(InspTime(PIF10ii), InspFlow(PIF10ii), 'rx');
        plot(ExpTime(TPEFi(i)), PEF(i), 'rv');
        plot(ExpTime(PEF10i), ExpFlow(PEF10i), 'rx');
        plot(ExpTime(PEF90i), ExpFlow(PEF90i), 'rx');
        plot(ExpTime(PEF90ii), ExpFlow(PEF90ii), 'rx');
        plot(ExpTime(PEF10ii), ExpFlow(PEF10ii), 'rx');
        title('Breath cycle with 10 and 90 marks');
    end
    
    %% teschler flattening/curvature index
    % Inspiratory airflow is scaled by mean inspiratory flow to one unit
    inspflow_forTeschler = InspFlow;
    inspflow_forTeschler = inspflow_forTeschler/mean(inspflow_forTeschler); %normalised to mean flow
    % The curvature index is a measure of the deviation from unit scaled
    % flow over the middle 50% of inspiratory time
    teschler(i) = trapz(time(TiI25(i):TiI75(i)),inspflow_forTeschler(TiI25(i):TiI75(i))) - ...
        trapz(time(TiI25(i):TiI75(i)),ones(length(time(TiI25(i):TiI75(i))),1));
    
    if ShowFigures
        figure(11);
        clf(figure(11));
        plot(InspTime, inspflow_forTeschler);hold on;
        refline(0,1);
        %should show vert lines indicating the area being analysed
        title('Teschler inspiration and reference line');
    end
    
    %% "sine fit" using simple polynomial curve fitting
    % fit 2nd degree polynomial to three set points, i.e start, max, & end
    % do insp first
    threeX = [InspTime(1); InspTime(round(length(InspTime)/2,0)); InspTime(end)];
    threeY = [InspFlow(1); max(InspFlow); InspFlow(end)];
    [pp,~,mu] = polyfit(threeX, threeY, 2);
    y2 = polyval(pp, InspTime,[],mu);
    if ShowFigures
        figure(12)
        clf(figure(12));
        subplot(2,1,1);
        plot(InspTime,InspFlow);hold on
        plot(InspTime, y2, 'r');
        title('Polyfit to inspiration');
        hold off
    end
    
    % find area difference between flow and curve
    SinI(i)=trapz(InspTime,InspFlow) - trapz(InspTime,y2);
    
    % repeat for the process for expiration
    threeX = [ExpTime(1); ExpTime(round(length(ExpTime)/2,0)); ExpTime(end)];
    threeY = [ExpFlow(1); min(ExpFlow); ExpFlow(end)];
    [pp,~,mu] = polyfit(threeX, threeY, 2);
    y2 = polyval(pp, ExpTime,[],mu);
    if ShowFigures
        figure(12)
        subplot(2,1,2);
        plot(ExpTime,ExpFlow);hold on
        plot(ExpTime, y2, 'r');
        title('Polyfit to expiration');
        hold off
    end
    % find area difference between flow and curve
    SinE(i)=trapz(ExpTime,ExpFlow) - trapz(ExpTime,y2);
 

    %% options and starting points for next few curve fitting options
    c=max(InspFlow)*1.1;
    N=length(InspFlow);
    x=(1:N)-round(N/2);

    %% inverted parabola
    y3=c*((1-(x.^2)/(N/2)^2));    
    invParabI(i)=trapz(x,InspFlow) - trapz(x,y3);    
    
    %% ellipseI
    y4=c*(((1-(x.^2)/(N/2)^2)).^(1/2)); 
    ellipseI(i)=trapz(x,InspFlow) - trapz(x,y4); 
    
    %% hypcosI, hyperbolic cosine arch 0<b<Inf
    b=10; 
    y5=c*((1-cosh((x/(N/2))*log(1+b+sqrt(2+b)*sqrt(b)))+b)/b);
    hypcosI(i)=trapz(x,InspFlow) - trapz(x,y5); 
     
    if ShowFigures
        figure(13);
        clf(figure(13));
        hold on;
        plot(x, InspFlow); hold on;
        plot(x, y3, x, y4, x, y4, x, y5);
        legend('Flow','InvP','Ellipse','HypCos');
        title('Various curve fits to inspiration');
    end
 
    
    %% mansourFL, 
    % fit 3rd order polynomial to insp, then find the 
    % slope of the polyfit at max flow (+ve = NIFL, -ve = IFL)
    PolyOrder = 3;
    range=(1:TPIFi(i));
    %range=(1:BB_I_endi(i)-BB_I_starti(i));
    [PolyfitValues,~,mu] = polyfit(InspTime(range),InspFlow(range),PolyOrder); 
    PolyLine = polyval(PolyfitValues,InspTime(range),[],mu); % Produce a line that represents the polyfit
    PolyLine_slope = gradient(PolyLine);% find the gradient of the PolyLine
    mansourFL(i) = PolyLine_slope(TPIFi(i)); % Slope at Max Flow
    if ShowFigures
        figure(14)
        clf(figure(14))
        subplot(2,1,1)
        plot(InspTime(range),InspFlow(range)); hold on;
        plot(InspTime(range), PolyLine);
        plot(InspTime(TPIFi(i)),PolyLine(TPIFi(i)),'ro');
        xlim([InspTime(range(1)),InspTime(range(end))]);
        title('3^{rd} order polyfit to inspiration');
        subplot(2,1,2)
        plot(InspTime(range), PolyLine_slope); hold on;
        plot(InspTime(TPIFi(i)),PolyLine_slope(TPIFi(i)),'ro');
        xlim([InspTime(range(1)),InspTime(range(end))]);
        refline(0,0);
        title('Slope of polyfit');
    end
    
    %% mansourRUA,
    % flow at min/max slope, and constant C from polynomial
    % Next, let's find Rua, as R=C/(F^2*t) where:
    %   C is the C from PolyValues above,
    %   F is the flow at the end of linear portion, as determined by
    %       min of derivative plot, if breath is unobstructed, or
    %       max of derivative plot, if breath is obstructed
    %   t is time
    if mansourFL(i)>0; % unobstructed, therefore looking for min
        [SlopeVal,SlopeIndexOffset] = min(PolyLine_slope);
    elseif mansourFL(i)<=0; % obstructed, therefore looking for max
        [SlopeVal,SlopeIndexOffset] = max(PolyLine_slope);
    end
    C = PolyfitValues(3);
    t = InspTime(SlopeIndexOffset)-InspTime(1);
    F = InspFlow(SlopeIndexOffset);
    mansourRUA(i)= abs(C /((F^2)*t)); %to do - work out why this gives bad results
    if ShowFigures
        figure(14)
        subplot(2,1,1); hold on;
        plot(InspTime(SlopeIndexOffset), PolyLine(SlopeIndexOffset), 'rx');
        xlim([InspTime(range(1)),InspTime(range(end))]);
        subplot(2,1,2); hold on;
        plot(InspTime(SlopeIndexOffset),SlopeVal,'rx');
        xlim([InspTime(range(1)),InspTime(range(end))]);
    end
    
    %% seriesIE50, ratio of insp flow to exp flow at mid volumes, (<0.97)
    inspcumsum = cumsum(InspFlow)*dt;
    halfinspvol = 0.5*(inspcumsum(end));
    Inspind = find(inspcumsum>halfinspvol,1,'first');
    inspflowatmidTV = flow(BB_I_starti(i)+Inspind);
    
    expcumsum = cumsum(ExpFlow);  
    halfexpvol = 0.5*(expcumsum(end));
    Expind = find(expcumsum<halfexpvol,1,'first');
    expflowatmidTV = flow(BB_E_starti(i)+Expind);
    
    if ~isempty(Inspind) && ~isempty(Expind)
        seriesIE50(i) = inspflowatmidTV/abs(expflowatmidTV);
    end
    
    %% morgensterns measures
    [~, ind] = max(InspFlow);
    inspflow_forMorgens = InspFlow/max(InspFlow); % normalised to max flow
    
    %figure();
    %plot(InspTime, InspFlow, 'b'); hold on
    %plot(InspTime, inspflow_forMorgens, 'r');
    
    VTiMod = cumsum(inspflow_forMorgens)*dt;
    morgensVPIFVTi(i) = VTiMod(ind);
    morgensV05(i) = VTiMod(round(length(VTiMod)*0.5));
    morgensV03(i) = VTiMod(round(length(VTiMod)*0.3));
    
%     if 0
%         % exp_fit from Morgenstern
%         [x, y] = prepareCurveData( x, y );
%         if 0
%             ExpValues = fit(x,y,'exp2'); % use standard two term exponential
%         else
%             % Use Morgenstern exponential formula
%             fit = @(b,x)  b(1).*(x.^b(2)).*(exp(b(3)*x)); % Function to fit
%             fcn = @(b) sum((fit(b,x) - y).^2);    % Least-Squares cost function
%             s = fminsearch(fcn, x);   % Minimise Least-Squares
%             ExpLine=fit(s,x);
%         end
%     end  
    
    % to do - morgensPIF

    morgens = [morgensPIF morgensVPIFVTi morgensV05 morgensV03];
    
    %% do the base measure for the relative measures
    % clark, area index, difference btw ref breath and test breath as a
    % %age of total area under test breath
    clark(i) = NaN; % to do - ref breath average shape???
    
    Ti(i) = time(BB_I_endi(i))-time(BB_I_starti(i));
    Te(i) = time(BB_E_endi(i))-time(BB_E_starti(i));
    
    % schneider, Ti/Ttot 
    schneider(i) = Ti(i)./Ttot(i);
    
    % mooney, Ti relative to ref breath
    mooney(i) = Ti(i);
    
    % wiriyaporn, Expiratory decay constant 
    if 0
        %fit decay model from peak exp flow to end expiration
        % 6 co-ords at 0.2 sec intervals after PEF, before end exp
        step = 10; % set to 20 to copy paper
        TimePts=ExpTime(TPEFi(i):step:end); 
        TimePtsLin=linspace(0,1,numel(TimePts));
        FlowPts=-1.*ExpFlow(TPEFi(i):step:end);

        if numel(FlowPts)>4
        [f3, gof3, ~] = fit(TimePtsLin', FlowPts', 'exp1');
            if gof3.rsquare > 0.95
                wiriyaporn(i) = -1/(f3.b);
                if show_figs
                    figure(12) % should centre and scale to display correctly
                    subplot(2,1,2); hold on;
                    %plot(f3, ExpTime(TPEFi(i):end)', ExpFlow(TPEFi(i):end)');
                    hold off;
                end
            end
        end
    end
    
    % morris, TPEF/Te
    morris(i) = TPEF(i)./Te(i);
    
    % ats, VPEF/VTe
    [~, ind]=min(ExpFlow);
    Exp_sum = cumsum(ExpFlow);
    VPEF = abs(Exp_sum(ind));
    ats(i) = VPEF./VTe(i);
    
    % kaplan, PIF/MIF
    kaplan(i) = PIF(i)./MIF(i);  
    
    % to do - TAA 
    % RIPS are currently not passed to this function, therefore, no can do.
    %[ phasedifference ] = ImmanuelTAA( , x2 )
    
end

%% BreathScores - Relative to reference breath
% only performed during OSA segment testing, 
% when Reference data is included in the function call
if nargin==8
clarkRef = nan(size(BB_I_starti,1),1);
schneiderRef = nan(size(BB_I_starti,1),1);
mooneyRef = nan(size(BB_I_starti,1),1);
wiriyapornRef = nan(size(BB_I_starti,1),1);
morrisRef = nan(size(BB_I_starti,1),1);
atsRef = nan(size(BB_I_starti,1),1);
kaplanRef = nan(size(BB_I_starti,1),1);
% Set up the reference breath from nargin(4)
% the required info for the breath being analysed is at nameRef(i), as
% calculatied above, and the relative comparison is at varargin(1).nameRef
RelativeValues = varargin(1);
% 53 - clarkRef, area index??
% 54 - schneiderRef, Ti/Ttot
% 55 - mooneyRef, Ti
% 56 - wiriyapornRef, Expiratory decay constant,
% 57 - morrisRef, TPEF/Te
% 58 - atsRef, VPEF/VTe
% 59 - kaplanRef, PTIF/MIF

for i = 1:size(BB_I_starti,1)  
    
    % clark, area index, difference btw ref breath and test breath as
    % age of total area under test breath
    clarkRef(i) = clark(i) / RelativeValues{1,1}(53);
    
    % schneider, Ti/Ttot relative to ref breath 
    schneiderRef(i) = schneider(i) / RelativeValues{1,1}(54);
    
    
    % mooney, Ti relative to ref breath
    mooneyRef(i) = mooney(i) / RelativeValues{1,1}(55);
    
    % wiriyaporn, Expiratory decay constant, if bigger than ref = IFL
    wiriyapornRef(i) = wiriyaporn(i) / RelativeValues{1,1}(56);
    
    % morris, TPEF/Te, relative to ref
    morrisRef(i) = morris(i) / RelativeValues{1,1}(57);
    
    % ats, VPEF/VTe, relative to ref
    atsRef(i) = ats(i) / RelativeValues{1,1}(58);
    
    % kaplan, PIF/MIF, relative to ref
    kaplanRef(i) = kaplan(i) /RelativeValues{1,1}(59);
end
RelScores = [clarkRef schneiderRef mooneyRef wiriyapornRef morrisRef atsRef kaplanRef];
end

%% Compile all the BreathData together
ClinicalScoring = [SS AS CS OS];

[TF, ModeSize, FaultList] = CheckDimensions( ...
          BB_I_starti,BB_I_endi,BB_E_starti,BB_E_endi, ...
          Ti', Te', TiTran', TeTran', Ttot', TPIF, TPEF);
if ~TF
    % report to user the at fault locations
    disp('Dimension mismatch for '); FaultList
    % then set variables in FaultList to ModeSize
    for i=1:size(FaultList,1)
        VarName=char(FaultList(i,:));
        VarData=eval(VarName);            
        VarData(end+1:ModeSize(1,1))=NaN;
        VarData=transpose(VarData);
        feval(@()assignin('caller',VarName, VarData)); 
    end
end

Times = [ time(BB_I_starti)', time(BB_I_endi)', time(BB_E_starti)', time(BB_E_endi)', ...
          BB_I_starti, BB_I_endi, BB_E_starti, BB_E_endi, ...
          Ti', Te', TiTran', TeTran', Ttot', TPIF, TPEF];
      
Flows = [MIF, MEF, PIF, PEF, MIF50, MEF50];

Volumes = [VT', VTi', VTe', Vdot', VI', VE'];

BreathRates = [ fR, FTi, RTi, DTi, FTe, RTe, DTe ];

BreathAreas = [ SinI, SinE, teschler, invParabI, ellipseI, hypcosI];

RefData = [clark, schneider, mooney, wiriyaporn, morris, ats, kaplan];
% this is perhaps a little misleading. this RefData is used for two purposes:
%(1) when analysing reference segments these are used to store averages of those values;
%(2) when analysing OSA segments, these are the base values which are then compared
%to the reference breath averages from use (1), and the result is stored in RelScores.

if nargin==8
    BreathScores = [mansourFL mansourRUA seriesIE50 morgens TAA RefData RelScores];
else
    BreathScores = [mansourFL mansourRUA seriesIE50 morgens TAA RefData];
end

BreathData = [ClinicalScoring Times Flows Volumes BreathRates BreathAreas BreathScores];

end

function [TF, varargout] = CheckDimensions(varargin)
Dims = NaN(nargin,2);
range = 1:nargin;
for i = range
    Dims(i,:) = size(varargin{1,i});
end
TF = (all(Dims(:,1) == Dims(1,1))) && (all(Dims(:,2) == Dims(1,2)));
varargout{1} = [];
varargout{2} = [];
if ~TF
    varargout{1} = [mode(Dims(:,1)), mode(Dims(:,2))];
    D1 = (Dims(:,1)~=Dims(1,1));
    D2 = (Dims(:,2)~=Dims(1,2));
    D3 = D1|D2;
    AtFault=range(D3);
    for i = AtFault
        varargout{2} = cat(1,varargout{2}, inputname(i));
    end
end
end