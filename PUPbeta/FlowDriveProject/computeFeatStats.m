function FeatStatsTable = computeFeatStats(subjectInfo, FeatNum, BreathDataTableFinal,...
                    BreathFLDataTableFinal, AHICutIdx)
                
%Compute: Mean, median, 90th percentile, 80th percentile, 10th percentile,
%20th percentile data
%Overall and for each subject

%Subject level
%Initialize vars
subcount = 1;
numsubs = sum(subjectInfo.TOTAL_SUP_AIH >=20);
Stats = zeros(numsubs, 7);
LabelsSub = zeros(numsubs, 1);
dAHIper = zeros(numsubs, 1);
SubsIn = cell(numsubs, 1);

%Run loop
for subnum = 1:size(subjectInfo,1)
    subIdx = BreathDataTableFinal.Subject == ...
        str2num(subjectInfo.MATFilename{subnum}(1:end-8)) &...
        AHICutIdx == 1;
    
    if sum(subIdx) == 0
        continue
    end
    
    
    %Mean
    Stats(subcount,1) = nanmean(BreathFLDataTableFinal{subIdx, FeatNum});
    
    %SD
    Stats(subcount,2) = nanstd(BreathFLDataTableFinal{subIdx, FeatNum});
    
    %Median
    Stats(subcount,3) = median(BreathFLDataTableFinal{subIdx, FeatNum}, 'omitnan');
        
    %90
    Stats(subcount,4) = prctile(BreathFLDataTableFinal{subIdx, FeatNum}, 90);
    
    %80
    Stats(subcount,5) = prctile(BreathFLDataTableFinal{subIdx, FeatNum}, 80);
    
    %20
    Stats(subcount,6) = prctile(BreathFLDataTableFinal{subIdx, FeatNum}, 20);
    
    %10
    Stats(subcount,7) = prctile(BreathFLDataTableFinal{subIdx, FeatNum}, 10);
    
    %Labels
    LabelsSub(subcount,1) = subjectInfo.labels(subnum);
    dAHIper(subcount,1) = subjectInfo.PercentReduction(subnum);
    
    %Subject cell array
    SubsIn{subcount,1} = subjectInfo.MATFilename{subnum}(1:end-8); 
    
    %Count
    subcount = subcount + 1;
end

%Create table
varnames = {'Subject', 'Mean', 'SD', 'Median',...
    'prc90', 'prc80', 'prc20', 'prc10'};
FeatStatsTable = cell2table([SubsIn, num2cell(Stats)]);
FeatStatsTable.Properties.VariableNames = varnames;
