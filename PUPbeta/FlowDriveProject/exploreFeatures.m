function exploreFeatures(subjectInfo, BreathDataTableFinal,...
                    BreathFLDataTableFinal, LabelArray, AHICutIdx)

                     
                     
%Compute new features
 % Absolute skewness (no directionality to it)
BreathFLDataTableFinal.AbsSkewInsp_O = abs(BreathFLDataTableFinal.SkewDistInsp_O);
BreathFLDataTableFinal.AbsSkewInsp_T = abs(BreathFLDataTableFinal.SkewDistInsp_T);

 % Absolute asymmetry (no directionality to it)
BreathFLDataTableFinal.AbsAsymInsp_O = abs(BreathFLDataTableFinal.AsymmetryInsp_O);
BreathFLDataTableFinal.AbsAsymInsp_T = abs(BreathFLDataTableFinal.AsymmetryInsp_T);

%Subject level means
%Initialize vars
subcount = 1;
numsubs = sum(subjectInfo.TOTAL_SUP_AIH >=20);
BreathFLDataSub = zeros(numsubs, size(BreathFLDataTableFinal,2));
LabelsSub = zeros(numsubs, 1);
dAHIper = zeros(numsubs, 1);
%Run loop
for subnum = 1:size(subjectInfo,1)
    subIdx = BreathDataTableFinal.Subject == ...
        str2num(subjectInfo.MATFilename{subnum}(1:end-8)) &...
        AHICutIdx == 1;
    
    if sum(subIdx) == 0
        continue
    end
    
    %MIFL data
    BreathFLDataSub(subcount,:) = nanmean(BreathFLDataTableFinal{subIdx,:},1);
    
    %Traits
    
    %Labels
    LabelsSub(subcount,1) = subjectInfo.labels(subnum);
    dAHIper(subcount,1) = subjectInfo.PercentReduction(subnum);
    
    SubsIn(subcount,1) = str2num(subjectInfo.MATFilename{subnum}(1:end-8)); 
    %Count
    subcount = subcount + 1;
end

%Absolute Skewness
figure('pos', [176 582 464 285]);
boxplot(table2array(BreathFLDataTableFinal(:, 189)), LabelArray, 'Symbol', '')
% ylim([-0.005 max(meanp70to100)+0.005])
hold on
scatter(LabelsSub+1, BreathFLDataSub(:,189))
xticklabels({'NR', 'R'})
ylabel("Absolute Skewness")
set(gca,'Box','On','FontSize',12)


%Absolute Asym
figure('pos', [176 582 464 285]);
boxplot(table2array(BreathFLDataTableFinal(:, 191)), LabelArray, 'Symbol', '')
% ylim([-0.005 max(meanp70to100)+0.005])
hold on
scatter(LabelsSub+1, BreathFLDataSub(:,191))
xticklabels({'NR', 'R'})
ylabel("Absolute Asymmetry")
set(gca,'Box','On','FontSize',12)

for ii = 1:size(BreathFLDataTableFinal,2)
    figure('pos', [176 582 464 285]);
    boxplot(table2array(BreathFLDataTableFinal(:, ii)), LabelArray, 'Symbol', '')
%     ylim([-0.1 max(BreathFLDataSub(:,ii))+0.1])
    hold on
    scatter(LabelsSub+1, BreathFLDataSub(:,ii))
    xticklabels({'NR', 'R'})
    ylabel(['Variable', num2str(ii)])
    set(gca,'Box','On','FontSize',12)
    
    close all
end
