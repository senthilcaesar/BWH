function FDuplicated2 = FindDuplicatedPt2(BreathDataTableNewE_)
% Part 2 does a double check of the breaths between window jumps
% This was mainly to correct of a bug in FindDuplicated
% Ideally, FindDuplicated should just incorporate this

% At the jumps between windows, check again for overlaps missed by
% FDuplicated

diff = BreathDataTableNewE_.Time0 - circshift(BreathDataTableNewE_.Time0, 1);
winChange = find(diff > 0)';
FDuplicated2 = zeros(size(BreathDataTableNewE_,1),1);

for i = winChange
    a = BreathDataTableNewE_.Time_start(i-1);
    b = BreathDataTableNewE_.Time_end(i-1);
    x = BreathDataTableNewE_.Time_start(i);
    y = BreathDataTableNewE_.Time_end(i);
    
    if (x - b) >= 0
        continue
    elseif a>=x&&a<=y && b>=x&&b<=y % a and b within x and y or aligned
        FDuplicated2(i) = 1;
    elseif a<x && b>x && b<=y %left shifted overlap
        overlap = (b-x)/(b-a);
        if overlap > 0.25
            FDuplicated2(i) = 1;
        end
    elseif a>=x && a<y && b>y % right shifted overlap
        overlap = (y-a)/(b-a);
        if overlap > 0.25
            FDuplicated2(i) = 1;
        end
    elseif a<x && b>y %x and y completely within a and b
        FDuplicated2(i) = 1;
    end 
end
