function T = LargeTableFromCellsOfTables(BreathDataTable)

%takes ~135 s for 135kx33 table

tic
%initializing Table T:
%find required N rows
temp2=[];
tempN=0;
for n=1:length(BreathDataTable)
temp = BreathDataTable{n};
for i=1:length(temp)
    if ~istable(temp{i})
        continue
    end
    if isempty(temp2)
    temp2 = temp{i};
    end
    tempN = tempN + size(temp{i},1);
end
end
VariableNames = [temp2.Properties.VariableNames,'subj','win']; 
A = zeros(tempN,length(VariableNames));
T = array2table(A,...
    'VariableNames',VariableNames);


tempN=1;
for n=1:length(BreathDataTable)
for i=1:length(BreathDataTable{n})
    if ~istable(BreathDataTable{n}{i})
        continue
    end
    xx = size(BreathDataTable{n}{i});
    B = table([n*ones(xx(1),1)],[i*ones(xx(1),1)],...
    'VariableNames',{'subj','win'});
    T(tempN:tempN+xx(1)-1,1:xx(2))=BreathDataTable{n}{i};
    T(tempN:tempN+xx(1)-1,end-1:end)=B;
    tempN = tempN + xx(1);
end
end
toc
