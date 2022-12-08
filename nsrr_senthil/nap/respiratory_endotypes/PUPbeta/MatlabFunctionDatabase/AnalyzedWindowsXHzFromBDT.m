function AnalyzedXHz = AnalyzedWindowsXHzFromBDT(BreathDataTable,Time)

%Function generates a signal that describes whether data were "analyzed" or not. 
%Used to find a modified "total sleep time" for event frequency (AHI) analysis, etc.
%e.g. it is inappropriate if only half the study has usable flow signals,
%but denominator time for AHI is the whole study. 
%AnalyzedXHz can be used to identify the "region of interest" period for analysis

for i=1:length(BreathDataTable)
    if istable(BreathDataTable{i})
        break
    end
end

StartTimeFirstWindow = BreathDataTable{i}.Time0(1);

WinStart=nan(length(BreathDataTable),1);
WinEnd=nan(length(BreathDataTable),1);

AnalyzedXHz = 0*Time;
for i=1:length(BreathDataTable)
    if istable(BreathDataTable{i})
        WinStart(i)=BreathDataTable{i}.Time_start(1);
        WinEnd(i)=BreathDataTable{i}.Time_end(end);
        AnalyzedXHz(Time>=WinStart(i)&Time<WinEnd(i))=1;
    end
end