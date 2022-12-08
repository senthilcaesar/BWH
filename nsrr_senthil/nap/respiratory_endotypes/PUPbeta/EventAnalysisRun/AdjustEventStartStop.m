function EvtsOut = AdjustEventStartStop(EvtsIn,BreathDataTable2,Fs)

try
    VI = BreathDataTable2.VI;
catch
    disp('There is no ventilation data')
    return
end


EvtStart = EvtsIn.RespT.EventStart;
EvtEnd = EvtsIn.RespT.EventStart + EvtsIn.RespT.EventDuration;
EvStIdx = nan(length(EvtStart),1);
EvEndIdx = nan(length(EvtStart),1);
rowsoutside = zeros(length(EvtStart),1);

for ii=1:length(EvtStart)
    if EvtStart(ii) < BreathDataTable2.Time_start(1) || EvtStart(ii) > BreathDataTable2.Time_start(end)
       rowsoutside(ii) = 1;
       continue
    end

    % find event start index and time
    TimeDiff = BreathDataTable2.Time_start - EvtStart(ii);
    TimeDiff(TimeDiff>0 | BreathDataTable2.FDuplicated2~=0) = nan;
    [~,EvStIdx(ii,1)] = max(TimeDiff);

    Windowtemp = BreathDataTable2.Time0(EvStIdx(ii,1));

    % find event end index and time
    TimeDiff = BreathDataTable2.Time_end - EvtEnd(ii);
    TimeDiff(TimeDiff<0 | BreathDataTable2.Time0 ~= Windowtemp) = nan;
    [~,EvEndIdx(ii,1)] = min(TimeDiff);
end

% remove nans
EvStIdx(isnan(EvStIdx)) = [];
EvEndIdx(isnan(EvEndIdx)) = [];
EvStIdx(isnan(EvStIdx)) = [];

% store original start times for later difference comparisons
EventStartOrig = BreathDataTable2.Time_start(EvStIdx);
EventEndOrig = BreathDataTable2.Time_end(EvEndIdx);

% % Get associated time variables (e.g. event start time)
% EvtStart2 = BreathDataTable2.Time_start(EvStIdx2);
% EvtEnd2 = BreathDataTable2.Time_end(EvEndIdx2);
% Window = BreathDataTable2.Time0(EvStIdx2);

BrSrchWin = 5; % size of the breath search window
EvStIdxNew = nan(length(EvEndIdx),1);
EvEndIdxNew = nan(length(EvEndIdx),1);

%%% Thresholds for deciding to move events %%%
% IF event end breath is VI > thresh1 then we're probably in
% the arousal so move it back. If breath after event end breath is VI < 0.7
% then event should be longer, move it forward. 
% There is a zone of do nothing between thresh1 and thresh2, wherein if
% event end is between (e.g 0.85), don't move it. 
thresh1 = 0.95; % WHAT SHOULD THESE BE??
thresh2 = 0.8;
thresh3 = 0.9;
thresh4 = 0.7;

format short
for evtnum = 1:length(EvStIdx)
    try
        if isnan(EvStIdx(evtnum))
            continue
        end
        
        %%% Adjust event end locations %%%
        if VI(EvEndIdx(evtnum)) > thresh1 % IF event end breath is in arousal - bring it back inside the event
            diffIdx = find(VI(EvEndIdx(evtnum):-1:EvStIdx(evtnum)) < thresh1, 1, 'first');
            if ~isempty(diffIdx)
                EvEndIdxNew(evtnum) = EvEndIdx(evtnum)-diffIdx+1;
            end       
        elseif  VI(EvEndIdx(evtnum)+1) < thresh2 % iimplies Ap/Hyp breaths after event end - so move it forward
            tempIdx = EvEndIdx(evtnum)+1;
            while VI(tempIdx) < thresh2 % loop one at time to check VI of each breath, until first breath >= 0.8
                EvEndIdxNew(evtnum) = tempIdx;
                tempIdx = tempIdx+1;
                if tempIdx - EvEndIdx(evtnum) > BrSrchWin % break the loop after 5 breaths
                    EvEndIdxNew(evtnum) = nan;
                    continue
                end
            end            
        end
        
        %%% Adjust event start location %%%
        if VI(EvStIdx(evtnum)) > thresh3 %IF event start breath has high VI, event started too early, move it forward
             diffIdx = find(VI(EvStIdx(evtnum):1:EvEndIdx(evtnum)) < thresh3, 1, 'first');
            if ~isempty(diffIdx)
                EvStIdxNew(evtnum) = EvStIdx(evtnum)-diffIdx+1;
            end    
        elseif VI(EvStIdx(evtnum)-1) < thresh4 %IF breath before event start breath has low VI, event started late, move it back
            tempIdx = EvStIdx(evtnum)-1;
            while VI(tempIdx) < thresh4 % loop one at time to check VI of each breath, until first breath >= 0.8
                EvStIdxNew(evtnum) = tempIdx;
                tempIdx = tempIdx-1;
                if  EvStIdx(evtnum)-tempIdx > BrSrchWin % break the loop after 5 breaths
                    EvStIdxNew(evtnum) = nan;
                    continue
                end
            end   
        end
    catch
        disp(['Problem adjusting events for Event number ', num2str(evtnum)])
    end
end
EvStIdx(~isnan(EvStIdxNew)) = EvStIdxNew(~isnan(EvStIdxNew));
EvEndIdx(~isnan(EvEndIdxNew)) = EvEndIdxNew(~isnan(EvEndIdxNew));

%%% Update Evts.RespT
EvtsOut = EvtsIn; 
EvtsOut.RespT.EventEnd = EvtsOut.RespT.EventStart + EvtsOut.RespT.EventDuration; % get original event end first
EvtsOut.RespTOrig = EvtsIn.RespT;
EvtsOut.RespT.EventStart(~rowsoutside) = BreathDataTable2.Time_start(EvStIdx);
EvtsOut.RespT.EventEnd(~rowsoutside) = BreathDataTable2.Time_end(EvEndIdx);
EvtsOut.RespT.EventDuration = EvtsOut.RespT.EventEnd - EvtsOut.RespT.EventStart;
EvtsOut.RespT.starttimesi = EvtsOut.RespT.EventStart*Fs;
EvtsOut.RespT.StartDiffT(~rowsoutside) = EventStartOrig - EvtsOut.RespT.EventStart(~rowsoutside);
EvtsOut.RespT.EndDiffT(~rowsoutside) = EventEndOrig - EvtsOut.RespT.EventEnd(~rowsoutside);
EvtsOut.RespT.RowsOutBrDataTbl = rowsoutside;

EvtsOut.Table1(EvtsOut.Table1.EventCodes > 1,:) = [];
EvtsOut.Table1 = vertcat(EvtsOut.Table1,EvtsOut.RespT(:,1:4));



% 
