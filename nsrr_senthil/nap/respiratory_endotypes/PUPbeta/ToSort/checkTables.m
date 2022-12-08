function [Table1, Table2] = ...
    checkTables(Table1, Table2)

% This function was added to fix issue with Subject 1830
% Basically, BreathFLDataTable had 3 nans where other tables had data
% This function replaces data with nans if no data exists in any of the
% three tables

for i=1:length(Table1{1})
    if sum(size(Table1{1}{i})==[1 1])==2 && ...
            sum(size(Table2{1}{i})==[1 1])==2
        
        continue
        
    elseif sum(size(Table1{1}{i})==[1 1])~=2 && ...
            sum(size(Table2{1}{i})==[1 1])~=2
        continue
    else
        Table1{1}{i} = nan; Table2{1}{i} = nan;
    end
end