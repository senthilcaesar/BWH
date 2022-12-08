function [FlowDriveOut,FlowDriveIn] = computeFlowDrive(FinalModelTable, BreathFLDataTable)

FlowDriveIn = nan(size(BreathFLDataTable,1), size(FinalModelTable,1));
FlowDriveIn(:,1) = FinalModelTable.Betas(1); % add intercept

for ii = 2:size(FinalModelTable)
    varname = FinalModelTable.Feature{ii};
    trExp = str2num([varname(end-2),'.',varname(end)]); %transform exponent
    
    % limit feature term to within upper and lower limits
    if 1
        temp = BreathFLDataTable.(varname(1:end-4));
        upper = FinalModelTable.upper(ii);
        lower = FinalModelTable.lower(ii);
        temp(temp<lower)=lower;
        temp(temp>upper)=upper;
        FlowDriveIn(:,ii) = sign(temp).*...
            abs(temp).^trExp.*...
            FinalModelTable.Betas(ii);
        
    else % no limit applied to feture terms
        FlowDriveIn(:,ii) = sign(BreathFLDataTable.(varname(1:end-4))).*...
            abs(BreathFLDataTable.(varname(1:end-4))).^trExp.*...
            FinalModelTable.Betas(ii);
        
    end
end

% add product of betas and corresponding feature terms
FlowDriveOut = sum(FlowDriveIn,2);

% Apply the limits for model output (0-1.5)
FlowDriveOut(FlowDriveOut<0)=0;
FlowDriveOut(FlowDriveOut>1.5)=1.5;