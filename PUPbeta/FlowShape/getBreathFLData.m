function [FL_Data] = getBreathFLData(BreathFLDataTable, n, attrib)
FL_Data = [];
for w=1:length(BreathFLDataTable{1,n}) % for all the windows for this pt
    if isempty(BreathFLDataTable{1,n}{1,w})
        continue
    end
    FL_Data_Temp = table2array(BreathFLDataTable{1,n}{1,w}(:,attrib));
    FL_Data = cat(1,FL_Data, FL_Data_Temp);
end
end