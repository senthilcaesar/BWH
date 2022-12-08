function [AllPositionsEventsTable,AllPositionsEventsTableN,SupineEventsTable,EvtDurTable] = AHIdata2Tbls(AHIdata)

states = {'AllSleep','NREM','REM','W','N1','N2','N3'};

AllPositionsEventsTable = array2table(reshape(AHIdata(73:128),8,7));
AllPositionsEventsTable.Properties.VariableNames = states;
rows = {'Duration','ArI','OAI','CAI','OHI','MAI','CHI','AHI'};
AllPositionsEventsTable.Properties.RowNames = rows;

AllPositionsEventsTableN = AllPositionsEventsTable;
AllPositionsEventsTableN{2:end,:} = AllPositionsEventsTable{2:end,:}.*AllPositionsEventsTable{1,:}/60;

AllPositionsEventsTable_i = array2table(reshape(73:128,8,7));
AllPositionsEventsTable_i.Properties.VariableNames = states;
AllPositionsEventsTable_i.Properties.RowNames = rows;

SupineEventsTable = array2table(reshape(AHIdata(1:56),8,7));
SupineEventsTable.Properties.VariableNames = states;
rows = {'Duration','ArI','OAI','CAI','OHI','MAI','CHI','AHI'};
SupineEventsTable.Properties.RowNames = rows;

SupineEventsTable_i = array2table(reshape(1:56,8,7));
SupineEventsTable_i.Properties.VariableNames = states;
SupineEventsTable_i.Properties.RowNames = rows;

EvtDurTable = array2table(reshape(AHIdata(129:184),8,7));
EvtDurTable.Properties.VariableNames = states;
rows = {'Ar','OA','CA','OH','MA','CH','OAH','AH'};
EvtDurTable.Properties.RowNames = rows;

EvtDurTable_i = array2table(reshape(129:184,8,7));
EvtDurTable_i.Properties.VariableNames = states;
rows = {'Ar','OA','CA','OH','MA','CH','OAH','AH'};
EvtDurTable_i.Properties.RowNames = rows;







