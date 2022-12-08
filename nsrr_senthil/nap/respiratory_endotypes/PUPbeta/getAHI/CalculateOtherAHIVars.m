function OtherVar = CalculateOtherAHIVars(AHIData,Var)

if strcmp(Var,'FHypopnea')

    NumApneas = (AHIData.OAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
        (AHIData.CAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
        (AHIData.MAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)); 

    NumHypopneas = (AHIData.OHIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
        (AHIData.CHIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)); 

    FHypopnea = NumHypopneas./(NumHypopneas + NumApneas);
    OtherVar = FHypopnea;
end

if strcmp(Var, 'ApneaDuration')
    AHIData = fillmissing(AHIData,'constant',0,'DataVariables',@isnumeric);
    ApneaDuration = AHIData.OA_DurationAllSleep + ...
        AHIData.MA_DurationAllSleep + AHIData.CA_DurationAllSleep;
    OtherVar = ApneaDuration;
end

if strcmp(Var, 'ApneaIndex')
    NumApneas = (AHIData.OAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
        (AHIData.CAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
        (AHIData.MAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)); 
    ApneaIndex = NumApneas./(AHIData.DurationAllPosAllSleep./60);
    OtherVar = ApneaIndex;
end