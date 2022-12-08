function criteria = CriteriaSelect(subj,settings,hyp,pos,evt)

%EvtCriteria(subj,settings,EvtsData{1}.RespT.state,EvtsData{1}.RespT.PositionRaw)

%hyp = EvtsData{1}.RespT.state
%pos = EvtsData{1}.RespT.PositionRaw
%evt = EvtsData{1}.RespT.EventCodes

%expected settings:
%settings.selectstate = 8

%settings.selectposition = 'All'
%settings.poscodesdatabase (ARRAY) must be imported already
%settings.protocol (1 column) must be importad already
%settings.positioncodesout must be importad already %[1     2     3     4     5     6     1     1     4]

switch settings.selectstate
    case 4
        hypok = [0 1 2]; %NREM
    case 1
        hypok = [2]; %N1
    case 2
        hypok = [1]; %N2
    case 3
        hypok = [0]; %N3
    case 5
        hypok = [3]; %REM
    case 8
        hypok = [0 1 2 3]; %NREM and REM
    case 9
        hypok = [0 1 2 3 4]; %NREM and REM and Wake
    case 0
        hypok = [4]; %Wake 
end

criteriaState = sum(hyp==hypok,2)>0; % criteria for sleep stages in breath analysis.

% settings.poscodesdatabase is from Pos database 
% Changed to text variable: Supine, Left, Right, Prone, Unknown, Upright or All (all positions considered)
positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{subj});
posTranslated = PositionTranslator(positioncodes,settings.positioncodesout,pos);
criteriaPos = PositionSelector(posTranslated,settings.selectposition);

 
switch settings.selecteventtype
    case 0 %all time, whether in events or not
        evtok = [0:20];
    case 1 %all event types (not arousals)
        evtok = [2:9];
    case 2 %hypopnea
        evtok = [4 ];
    case 3 %just obstructive apnoeas
        evtok = [2];
end
criteriaEventType = sum(evt==evtok,2)>0;

criteria = criteriaState&criteriaPos&criteriaEventType;


