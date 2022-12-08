function BreathDataTableOut = LabelBreaths(BreathDataTable,Evts)

BreathDataTableOut = BreathDataTable;
for ii = 1:length(BreathDataTable)
    if istable(BreathDataTableOut{1,ii})
        BreathDataTableOut{1,ii} = Evts2BrTable(Evts.RespT,BreathDataTableOut{1,ii},'SiteOfCollapse');
        BreathDataTableOut{1,ii}  = Evts2BrTable(Evts.MouthOpenClosed,BreathDataTableOut{1,ii},'MouthOpenClosed');
        BreathDataTableOut{1,ii}  = Evts2BrTable(Evts.SideSleep,BreathDataTableOut{1,ii},'SideSleep');
    end
end
