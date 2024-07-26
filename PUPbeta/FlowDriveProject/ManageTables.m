function FDuplicated = FindDuplicated(BreathDataTable)

T = BreathDataTable;
for i=1:length(T{1})
    if ~istable(T{1}{i})
        continue
    end 
    Win = ones(size(T{1}{i},1),1)*i;
    CenterTime = (T{1}{i}.Time_end(end) - T{1}{i}.Time_start(1))/2 + T{1}{i}.Time_start(1);
    TimeFromCenter =  abs(T{1}{i}.Time_start - CenterTime);
    T{1}{i}.Win = Win;
    T{1}{i}.TimeFromCenter = TimeFromCenter;
end
BreathDataTableNew = T; %has TimeFromCenter and Win (number) 

%Expand Tables
BreathDataTableNewE = expandtable(BreathDataTableNew);

clear TimeFromCenter Win T CenterTime

%% Find duplicated breaths
M=size(BreathDataTableNewE,1);
FDuplicated = zeros(size(BreathDataTableNewE,1),1);

for i=1:M
    %do we want to reject breath i as a duplicate? compare with others
    id = BreathDataTableNewE.TimeFromCenter(i);
    a = BreathDataTableNewE.Time_start(i);
    b = BreathDataTableNewE.Time_end(i);
    x = BreathDataTableNewE.Time_start;
    y = BreathDataTableNewE.Time_end;
    test = x>=a&x<=b|y>=a&y<=b|a>=x&a<=y|b>=x&b<=y;
    I = find(test & (id > BreathDataTableNewE.TimeFromCenter));
    I=I(:)';
    if ~isempty(I)
    for j=I
        if i==j
            continue
        end
        x = BreathDataTableNewE.Time_start(j);
        y = BreathDataTableNewE.Time_end(j);
        if x>=a&&x<=b||y>=a&&y<=b||a>=x&&a<=y||b>=x&&b<=y
            %no overlap
            if a>=x&&a<=y && b>=x&&b<=y
                FDuplicated(i) = FDuplicated(i) + 1;
            elseif x>=a&&x<=b
                FDuplicated(i) = FDuplicated(i) + (b-x)/(b-a);
            elseif y>=a&&y<=b
                FDuplicated(i) = FDuplicated(i) + (y-a)/(b-a);
            elseif x>=a&&x<=b && y>=a&&y<=b
                FDuplicated(i) = FDuplicated(i) + (y-x)/(b-a);
            end
        end
    end
    end
end

BreathDataTableNewE.FDuplicated = FDuplicated;

clear a b x y i I id j M test         
