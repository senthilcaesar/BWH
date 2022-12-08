function [AHI_perPT, AHI_perPT_table] = getAHI_postanalysis()

%load('..\FeatureSpaces\FlowDrive_only125Hz_All.mat', 'AHIdata');
load('C:\PSG_Data\FlowDrive\Analyzed\FlowDrive_only25Hz_Normalized_All.mat', 'AHIdata');

AHI_perPT = NaN(54,2);

for pt=1:54
    if size(AHIdata{1,pt},2)==184
        
        %displaytext=['Total AHI: ' num2str(AHIdata{1,pt}(58),4) ' per hr']; disp(displaytext);
        %displaytext=['NREM supine AHI: ' num2str(AHIdata{1,pt}(16),4) ' per hr']; disp(displaytext);
        
        AHI_perPT(pt,1) = AHIdata{1,pt}(58); %Total AHI
        AHI_perPT(pt,2) = AHIdata{1,pt}(16); %NREM supine AHI
    else
        continue
    end
end
AHI_perPT_table = array2table(AHI_perPT, 'VariableNames', {'TotalAHI','NREMSupineAHI'});

end

