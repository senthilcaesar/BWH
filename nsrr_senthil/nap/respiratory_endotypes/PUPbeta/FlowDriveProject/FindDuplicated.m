function FDuplicated = FindDuplicated(T)

if sum(size(T))==2
    if iscell(T{1})
        T=T{1};
    end
end

if isa(T,'table')
    T_ = T;
    T = cell(1,1);
    T{1} = T_;
end

for i=1:length(T)
    if ~istable(T{i})
        continue
    end 
    Win = ones(size(T{i},1),1)*i;
    CenterTime = (T{i}.Time_end(end) - T{i}.Time_start(1))/2 + T{i}.Time_start(1);
    TimeFromCenter =  abs(T{i}.Time_start - CenterTime);
    T{i}.Win = Win;
    T{i}.TimeFromCenter = TimeFromCenter;
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
    
    test = x>=a&x<=b|y>=a&y<=b|a>=x&a<=y|b>=x&b<=y; %finds breaths with which breath i overlaps
    I = find(test & (id > BreathDataTableNewE.TimeFromCenter)); % find overlapping breaths that are closer to the center of the window than breath i
    I=I(:)';
    
    if ~isempty(I)
    for j=I
        if i==j % ignore identical breaths
            continue
        end
        
        x = BreathDataTableNewE.Time_start(j);
        y = BreathDataTableNewE.Time_end(j);
        
        if x>=a&&x<=b||y>=a&&y<=b||a>=x&&a<=y||b>=x&&b<=y % this double checks that overlap exists
            
            % Identify the degree of overlap (this prevents us from calling a breath a duplicate if there is only a small overlap)
            if a>=x&&a<=y && b>=x&&b<=y % breath i completely within breath j or perfectly aligned
                FDuplicated(i) = FDuplicated(i) + 1;
            elseif x>a&&x<b&&b<=y % breath i left shifted overlap with j 
                FDuplicated(i) = FDuplicated(i) + (b-x)/(b-a);
            elseif y>a&&y<b&&a>=x % breath i right shifted overlap with j
                FDuplicated(i) = FDuplicated(i) + (y-a)/(b-a);
            elseif x>a&&x<b && y>a&&y<b % breath j within breath i
                FDuplicated(i) = FDuplicated(i) + (y-x)/(b-a);
            end
            
        end
    end
    end
end

BreathDataTableNewE.FDuplicated = FDuplicated; % nnz(FDuplicated)

clear a b x y i I id j M test         
