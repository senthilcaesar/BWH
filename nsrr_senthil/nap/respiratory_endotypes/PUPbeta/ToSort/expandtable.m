function [Tnew, included] = expandtable(T)

if sum(size(T))==2 % fixes nested cell situation
    if iscell(T{1})
        T=T{1};
    end
end

if isa(T,'table')
    T_ = T;
    T = cell(1,1);
    T{1} = T_;
end

NewTableExists=0;
count = 0; countempty = 0;
for i=1:length(T)
    % this just counts the number of empty cells
    if sum(size(T{i})==[1 1])==2 || isempty(T{i})%size = [1,1] %second 
        countempty = countempty+1;
        continue
    end 
    
    % Occasionally, tables in T are a different size, this fixes that
    if i>1 && size(T{i-1},2) - size(T{i},2) > 1
       varsmissing = setdiff(T{i-1}.Properties.VariableNames,...
           T{i}.Properties.VariableNames);
       varcols = find(ismember(T{i-1}.Properties.VariableNames, varsmissing));
       for vars = 1:length(varsmissing) %find location of missing vars
           T{i}.(varsmissing{vars}) = nan(size(T{i},1),1);
       end
    end
    
    % Add a unique ID to the table
    % Consists of concatenated Windonw number*1000 + breath number
    BrNum = 1:size(T{i},1);
    UniqueID = i*1000 + BrNum;
    if ~isempty(UniqueID)
        T{i}.UniqueID = UniqueID';
    end
    
    if ~NewTableExists
        Tnew = T{i};
        NewTableExists=1;
        count = count+1;
        included(count,1) = i;
    else
        %add code to merge tables with different columns of data
        Tnew = [Tnew;T{i}];
        count = count+1;
        included(count,1) = i;
    end
    
    
end

% if all cells were empty or NaN, then return NaNs
if countempty==length(T)
    Tnew = NaN;
    included = NaN;
end