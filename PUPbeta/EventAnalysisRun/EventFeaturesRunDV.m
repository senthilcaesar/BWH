% Get published versions of VE features
function [EvtFtrs]=EventFeaturesRunDV(Boxes,Ensembles,RespT)
EvtFtrs = struct();

%% calculate event depth
EventDepth = struct();
LastEvtBr = 101;

% Event depth boundaries at 90
EventDepth.EventDepth_All90 = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'All','All','mean');
EventDepth.EventDepth_All90_NREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'All','NREM','mean');
EventDepth.EventDepth_All90_REM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'All','REM','mean');
EventDepth.EventDepth_All90_Supine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Supine','All','mean');
EventDepth.EventDepth_All90_NonSupine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Nonsupine','All','mean');
EventDepth.EventDepth_All90_SupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Supine','NREM','mean');
EventDepth.EventDepth_All90_SupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Supine','REM','mean');
EventDepth.EventDepth_All90_NonSupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Nonsupine','NREM','mean');
EventDepth.EventDepth_All90_NonSupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Nonsupine','REM','mean');                         

% Event depth boundaries at 80
EventDepth.EventDepth_All80 = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'All','All','mean');
EventDepth.EventDepth_All80_NREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'All','NREM','mean');
EventDepth.EventDepth_All80_REM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'All','REM','mean');
EventDepth.EventDepth_All80_Supine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Supine','All','mean');
EventDepth.EventDepth_All80_NonSupine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Nonsupine','All','mean');
EventDepth.EventDepth_All80_SupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Supine','NREM','mean');
EventDepth.EventDepth_All80_SupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Supine','REM','mean');
EventDepth.EventDepth_All80_NonSupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Nonsupine','NREM','mean');
EventDepth.EventDepth_All80_NonSupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Nonsupine','REM','mean');   

% Event depth boundaries at 70
EventDepth.EventDepth_All70 = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'All','All','mean');
EventDepth.EventDepth_All70_NREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'All','NREM','mean');
EventDepth.EventDepth_All70_REM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'All','REM','mean');
EventDepth.EventDepth_All70_Supine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Supine','All','mean');
EventDepth.EventDepth_All70_NonSupine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Nonsupine','All','mean');
EventDepth.EventDepth_All70_SupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Supine','NREM','mean');
EventDepth.EventDepth_All70_SupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Supine','REM','mean');
EventDepth.EventDepth_All70_NonSupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Nonsupine','NREM','mean');
EventDepth.EventDepth_All70_NonSupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Nonsupine','REM','mean'); 

% Event depth boundaries at 90
EventDepth.EventDepth_Med90 = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'All','All','median');
EventDepth.EventDepth_Med90_NREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'All','NREM','median');
EventDepth.EventDepth_Med90_REM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'All','REM','median');
EventDepth.EventDepth_Med90_Supine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Supine','All','median');
EventDepth.EventDepth_Med90_NonSupine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Nonsupine','All','median');
EventDepth.EventDepth_Med90_SupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Supine','NREM','median');
EventDepth.EventDepth_Med90_SupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Supine','REM','median');
EventDepth.EventDepth_Med90_NonSupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Nonsupine','NREM','median');
EventDepth.EventDepth_Med90_NonSupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.9,0.9,'Nonsupine','REM','median');                         

% Event depth boundaries at 80
EventDepth.EventDepth_Med80 = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'All','All','median');
EventDepth.EventDepth_Med80_NREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'All','NREM','median');
EventDepth.EventDepth_Med80_REM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'All','REM','median');
EventDepth.EventDepth_Med80_Supine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Supine','All','median');
EventDepth.EventDepth_Med80_NonSupine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Nonsupine','All','median');
EventDepth.EventDepth_Med80_SupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Supine','NREM','median');
EventDepth.EventDepth_Med80_SupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Supine','REM','median');
EventDepth.EventDepth_Med80_NonSupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Nonsupine','NREM','median');
EventDepth.EventDepth_Med80_NonSupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.8,0.8,'Nonsupine','REM','median');   

% Event depth boundaries at 70
EventDepth.EventDepth_Med70 = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'All','All','median');
EventDepth.EventDepth_Med70_NREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'All','NREM','median');
EventDepth.EventDepth_Med70_REM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'All','REM','median');
EventDepth.EventDepth_Med70_Supine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Supine','All','median');
EventDepth.EventDepth_Med70_NonSupine = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Nonsupine','All','median');
EventDepth.EventDepth_Med70_SupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Supine','NREM','median');
EventDepth.EventDepth_Med70_SupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Supine','REM','median');
EventDepth.EventDepth_Med70_NonSupineNREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Nonsupine','NREM','median');
EventDepth.EventDepth_Med70_NonSupineREM = CalculateEventDepth(Boxes,RespT,LastEvtBr,0.7,0.7,'Nonsupine','REM','median'); 

% %% Various conidition (Sleep states/position) - Mean method only
% % NREM all positions
% VIidx = RespT.Epochs ~=3;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_NREMAllPos = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_NREMAllPos = nan;
% end
% 
% % REM all position
% VIidx = RespT.Epochs ==3;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_REMAllPos = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_REMAllPos = nan;
% end
% 
% % Supine all sleep
% VIidx = RespT.Position==1;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_SupineAllSleep = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_SupineAllSleep = nan;
% end
% 
% % Nonsupine all sleep
% VIidx = RespT.Position~=1;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_NonSupineAllSleep = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_NonSupineAllSleep = nan;
% end
% 
% % Supine NREM
% VIidx = RespT.Position==1 & RespT.Epochs~=3;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_SupineNREM = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_SupineNREM = nan;
% end
% 
% % Supine REM
% VIidx = RespT.Position==1 & RespT.Epochs==3;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx)> 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_SupineREM = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_SupineREM = nan;
% end
% 
% % NonSupine NREM
% VIidx = RespT.Position~=1 & RespT.Epochs~=3;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_NonSupineNREM = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_NonSupineNREM = nan;
% end
% 
% % NonSupine REM
% VIidx = RespT.Position~=1 & RespT.Epochs==3;
% VI = nanmean(Boxes.VI(VIidx,:));
% if sum(VIidx) > 3 && sum(isnan(VI)) < length(VI)
%     EventDepth.EventDepth_NonSupineREM = CalculateEventDepth(VI,LastEvtBr,0.9,0.9);
% else
%     EventDepth.EventDepth_NonSupineREM = nan;
% end

EvtFtrs.EventDepth = EventDepth;

%% calculate desaturation slope
% find max pre event
meanSpO2 = Ensembles.SpO2;

% refine pre event baseline
ii = 100;
while meanSpO2(ii) < meanSpO2(ii-1)
    ii = ii-1;
    if ii == 92  
        break
    end
end
PreEvtBslnIdx = ii;
PreEvtBsln = meanSpO2(PreEvtBslnIdx);

[~, minSpO2idxTemp] = min(meanSpO2(PreEvtBslnIdx:PreEvtBslnIdx+50));
minSpO2idx = PreEvtBslnIdx + minSpO2idxTemp-1;

MeanEvSat = meanSpO2(PreEvtBslnIdx:minSpO2idx);
Time = (PreEvtBslnIdx:minSpO2idx)';
DesatSlope = CalculateDesatSlope(MeanEvSat,Time,1);
DesatSlope.SatAtBsln = PreEvtBsln;

EvtFtrs.DesatSlope = DesatSlope;
end

function EventDepth = CalculateEventDepth(Boxes,RespT,LastEvtBr,RVentLim,LVentLim,pos,sleep,method)
if strcmp(pos,'All') && strcmp(sleep,'All')
    VIidx = true(size(Boxes.VI,1),1);
elseif strcmp(pos,'All') && strcmp(sleep,'NREM')
    VIidx = RespT.Epochs~=3;
elseif strcmp(pos,'All') && strcmp(sleep,'REM')
    VIidx = RespT.Epochs==3;
elseif strcmp(pos,'Supine') && strcmp(sleep,'All')
    VIidx = RespT.Position==1;
elseif strcmp(pos,'Nonsupine') && strcmp(sleep,'All')
    VIidx = RespT.Position~=1;
elseif strcmp(pos,'Supine') && strcmp(sleep,'NREM')
    VIidx = RespT.Position==1 & RespT.Epochs~=3;    
elseif strcmp(pos,'Supine') && strcmp(sleep,'REM')
    VIidx = RespT.Position==1 & RespT.Epochs==3;        
elseif strcmp(pos,'Nonsupine') && strcmp(sleep,'NREM')
    VIidx = RespT.Position~=1 & RespT.Epochs~=3;        
elseif strcmp(pos,'Nonsupine') && strcmp(sleep,'REM')
    VIidx = RespT.Position~=1 & RespT.Epochs==3;    
end

if strcmp(method,'mean')
    meanVE = nanmean(Boxes.VI(VIidx,:),1);
elseif strcmp(method,'median')
    meanVE = nanmedian(Boxes.VI(VIidx,:),1);
end

if sum(isnan(Boxes.VI)) == length(Boxes.VI)
    EventDepth = nan;
    return
end

EupBreathTemp = find(meanVE(LastEvtBr:-1:1) <= RVentLim, 1, 'first');
EupBreathR_mean = LastEvtBr - EupBreathTemp+1;

EupBreathTemp = find(meanVE(EupBreathR_mean:-1:1) >= LVentLim, 1, 'first');
EupBreathL_mean = EupBreathR_mean - EupBreathTemp+1;

if ~isempty(EupBreathR_mean) && ~isempty(EupBreathL_mean) 
    EventDepth = 1 - nanmean(meanVE(EupBreathL_mean:EupBreathR_mean));
    
    % check that minimum value is within event indices
    [~,minIdx] = min(meanVE);
    if minIdx < EupBreathR_mean && minIdx > EupBreathL_mean
    else
        disp('WARNING: Minimum ventilation is outside the bounds of the event in event depth calculation')
    end
else
    EventDepth = nan;
end


        
end

function DesatSlope = CalculateDesatSlope(MeanEvSat,Time,plotfig)

DesatSlope = struct();
%%% Fit a line from 25th to 75th percentile
try
    Sat25to75 = prctile(MeanEvSat,[25 75]);
    MeanEvSat2 = MeanEvSat(MeanEvSat>=Sat25to75(1) & MeanEvSat<=Sat25to75(2));
    Time2 = Time(MeanEvSat>=Sat25to75(1) & MeanEvSat<=Sat25to75(2));

    % Fit a line
    p = polyfit(Time2,MeanEvSat2,1);
    DesatSlope.Lin2575Slp = p(1); % linear slope
%     R = corrcoef(Time2,MeanEvSat2);
%     RsqLinear2575 = R(1,2).^2; % Rsquared of fit
%     SatAt5025to75 = 100-prctile(MeanEvSat2,50);

%     if plotfig
%         lineToPlot=polyval(p,Time2);
%         figure(23), plot(Time2,lineToPlot,'b--','LineWidth',2)
%     end
catch
    DesatSlope.Lin2575Slp = nan;
%     RsqLinear2575 = nan;
%     SatAt5025to75 = nan;
end
%% Fit a line from 10th to 90th percentile
try
    Sat10to90 = prctile(MeanEvSat,[10 90]);
    MeanEvSat2 = MeanEvSat(MeanEvSat>=Sat10to90(1) & MeanEvSat<=Sat10to90(2));
    Time2 = Time(MeanEvSat>=Sat10to90(1) & MeanEvSat<=Sat10to90(2));

    % Fit a line
    p = polyfit(Time2,MeanEvSat2,1);
    DesatSlope.Lin1090Slp = p(1); % linear slope
%     R = corrcoef(Time2,MeanEvSat2);
%     RsqLinear1090 = R(1,2).^2; % Rsquared of fit
%     SatAt5010to90 = 100-prctile(MeanEvSat2,50);

%     if plotfig
%         lineToPlot=polyval(p,Time2);
%         figure(23), plot(Time2,lineToPlot,'-.','Color',[0.5 0.5 0.5],'LineWidth',2)
%     end
catch
    DesatSlope.Lin1090Slp = nan;
%     RsqLinear1090 = nan;
%     SatAt5010to90 = nan;
end
end
