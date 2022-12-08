function [BB, BB_] = getAllorSleepOnlyBB(SleepOnly, IncludePnasal, PtData, PtData_2)

if SleepOnly
    BB = (PtData.Hypnog<4)&(PtData.Ar==0); % sleep only, nnz(BB)
    if IncludePnasal
        BB_ = (PtData_2.Hypnog<4)&(PtData_2.Ar==0); % sleep only, nnz(BB_)
    else
        BB_ = [];
    end
else
    BB = true(height(PtData),1);
    if IncludePnasal
        BB_ = (~isnan(PtData_2{:,1}));
    else
        BB_ = [];
    end
end

end