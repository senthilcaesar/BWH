function FHypopnea = CalculateFHypopnea(AHIData)

NumApneas = (AHIData.OAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
    (AHIData.CAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
    (AHIData.MAIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)); 

NumHypopneas = (AHIData.OHIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)) + ...
    (AHIData.CHIAllPosAllSleep.*(AHIData.DurationAllPosAllSleep./60)); 

FHypopnea = NumHypopneas./(NumHypopneas + NumApneas);