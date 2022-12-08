function BreathDataTable = Evts2BrTable(EvtsTable, BreathDataTable,ColumnLabel)

EvtsTable.EventEnd = EvtsTable.EventStart + EvtsTable.EventDuration;
NewVar = cell(length(BreathDataTable.Time0),1);

for evtnum = 1:size(EvtsTable,1)
%     BrIdx = BreathDataTable.Time_start > EvtsTable.EventStart(evtnum) & ...
%          BreathDataTable.Time_end < EvtsTable.EventEnd(evtnum);
    BrIdx = BreathDataTable.Time_mid > EvtsTable.EventStart(evtnum) & ...
     BreathDataTable.Time_mid < EvtsTable.EventEnd(evtnum);
    NewVar(BrIdx) = EvtsTable.EventName(evtnum);
end

BreathDataTable.(ColumnLabel) = NewVar;
