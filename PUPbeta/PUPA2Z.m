function PUPA2Z(Mrange)

MrangeOverride = Mrange;
ConvertToMat(MrangeOverride);
AnalysisN(MrangeOverride);

SummaryAnalysisN(MrangeOverride);
EventAnalysis(MrangeOverride);

getData(MrangeOverride);
getOvernightHR(MrangeOverride);
getFL(MrangeOverride);
