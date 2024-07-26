function BreathDataTableNewE_ = RemoveDuplicatesTidyTable(BreathDataTableNewE,...
                                    Time_start, Time_end, FDuplicated, cut)
%table needs Fduplicated, Time_start, Time_end
minoverlap=1; %seconds of acceptable error
BreathDataTableNewE_ = BreathDataTableNewE;

BreathDataTableNewE_(FDuplicated>cut,:) = [];
Time_start(FDuplicated>cut,:) = [];
Time_end(FDuplicated>cut,:) = [];

Tjump = [Time_start] - [NaN;Time_end(1:end-1)];

for i=1:size(BreathDataTableNewE_,1)
    if abs(Tjump(i))>minoverlap
        BreathDataTableNewE_ = ...
            [BreathDataTableNewE_(1:i,:); ...
            BreathDataTableNewE_(i:end,:)];
        Tjump = ...
            [Tjump(1:i); ...
            Tjump(i)*NaN; ...
            Tjump(i+1:end)];
        BreathDataTableNewE_(i,:) = array2table(NaN);
        Time_start(i) = Time_end(i-1);
        
        i=i+1;
        
        BreathDataTableNewE_ = ...
            [BreathDataTableNewE_(1:i,:); ...
            BreathDataTableNewE_(i:end,:)];
        Tjump = ...
            [Tjump(1:i); ...
            Tjump(i)*NaN; ...
            Tjump(i+1:end)];
        BreathDataTableNewE_(i,:) = array2table(NaN);
        i=i+1;
    end
end
