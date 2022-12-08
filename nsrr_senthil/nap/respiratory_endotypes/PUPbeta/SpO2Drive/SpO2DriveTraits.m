function SpO2DriveTraitsT = SpO2DriveTraits(MrangeOverride)
global settings AMasterSpreadsheet

%First run StartHere for a project, then...
M=size(settings.patients,1);
Mrange=1:M;
if exist('MrangeOverride')
    Mrange=MrangeOverride;
end

dir=[settings.workdir 'Analyzed\'];
SpO2DriveTraitsT = table([]);

for i=1:length(Mrange)
    clear BreathDataTable BreathFLDataTable
    try
        success(i)=1; 
        filen = [dir settings.savename '_' num2str(Mrange(i))];
        load(filen,'BreathDataTable','BreathFLDataTable');
        
        if 1
            [T1,T2]=GetNonOvlappedVE(BreathDataTable); %T1 has overlapping data
            
            temp = T1.spo2(T1.hypnog_B==4);
            I = temp>prctile(temp,10)&temp<prctile(temp,90);
            spo2wakemean10to90 = nanmean(temp);
            
            run SpO2DriveOne
        end
    end
end
