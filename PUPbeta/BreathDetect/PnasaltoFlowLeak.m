function [F,Parameters_out,FlowSignal]=PnasaltoFlowLeak(IEratio,Ydata1,Time,dt,exponent,plotfigs)
% global settings
leakdelta=0; %is the zero flow baseline: initial value is zero, will be updated below

%%
for w=1:2 %repeat in loop once to get closer estimate of zero flow baseline
FlowSignal = Ydata1;
%FlowSignal(5000:25000) = -0.03;

% exponent=settings.scalingexponent; %0.67 is default for nasal pressure


FlowSignal=FlowSignal - leakdelta; %subtract current guess zeroflow baseline

%linearize the signal
FlowSignal(FlowSignal>0)=(FlowSignal(FlowSignal>0).^(exponent))/(IEratio^0.5);
FlowSignal(FlowSignal<0)=(-(-FlowSignal(FlowSignal<0)).^(exponent))*(IEratio^0.5);

%perform breath detection, which yields a new guess at the zeroflowbaseline
    [VT,VTi,VTe,leaktotal]=VflowanalysisFastLeak(FlowSignal,Time,dt,1-plotfigs); 
   
%update zero flow baseline    
    if leaktotal>0
        leakdelta = leakdelta + (leaktotal*(IEratio^0.5))^(1/exponent); %undo transofmration for the zero flow baseline.
    else
        leakdelta = leakdelta - (-leaktotal/(IEratio^0.5))^(1/exponent);
    end
end
%leakdelta is found.

%% For current estimate of IEratio (and accompnying guess of leakdelta), find estimate of accuracy (how much leak is there for large breaths, y2)
    %toc
    FVTi = VTi-VTe; %is the "leak" parameter. Ideally is zero. 
    Nbreaths = length(VTe);
    minbreaths = 5;
    prclow = 10;
        if Nbreaths*prclow/100<minbreaths
            prclow = minbreaths/Nbreaths;
        end
    VTXtile = prctile(VT,prclow);
    prclow2 = 67;
    VTXtile2L = prctile(VT,prclow2);
    if VTXtile==0
        VTXtile=min(VT(VT>0));
    end
    if VTXtile2L==0
        VTXtile2L=min(VT(VT>0));
    end
    VTXtile2 = max(VT); %75
    VTbelowVTXtile=VT<=VTXtile; %Index of small breaths (lowest 10%)
    VTwithinVTXtile2=VT>=VTXtile2L&VT<=VTXtile2; %Index of large breaths (top third largest)
    
    FVTilowestXtile = median(FVTi(VTbelowVTXtile==1)); %leak for small breaths (y1)
    VTlowestXtile = median(VT(VTbelowVTXtile==1)); %size of small breaths (x1)
    FVTisecondlowestXtile = median(FVTi(VTwithinVTXtile2==1)); %leak for large breaths (y2)
    VTsecondlowestXtile = median(VT(VTwithinVTXtile2==1));  %size of large breaths (x2)
    
    
    if plotfigs||0
        figure(99); plot(VT,FVTi,'.',[VTlowestXtile VTsecondlowestXtile],[FVTilowestXtile FVTisecondlowestXtile],'ko-');
    end
    projectedFVTizeroVT = FVTilowestXtile-VTXtile*(FVTisecondlowestXtile-FVTilowestXtile)/(VTXtile2-VTXtile); %y intercept of plot
    if 0
        F = abs(projectedFVTizeroVT)+abs(FVTisecondlowestXtile); %
    else
        F = FVTisecondlowestXtile; %is the error term we are aiming to minimize
    end
    Parameters_out = [leakdelta FVTilowestXtile FVTisecondlowestXtile]; %leak for large breaths (y2) is kept as a measure of how accurate IEratio guess is
    %toc




