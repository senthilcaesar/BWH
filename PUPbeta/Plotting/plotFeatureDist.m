function plotFeatureDist(BreathFLDataTableFinal, featNum, labelArray, AHICutIdx,...
    BreathDataTableFinal)

figfolder = 'J:\PEOPLE\POST DOC\VENA\Flow shape\Results and Processed Data\Figures\PSG MM Scope';

%Responder breath
RBreathIdx = labelArray == 1 & AHICutIdx == 1;
RBreathFeat = BreathFLDataTableFinal(RBreathIdx,featNum);

%Non-responder breath indices
NRBreathIdx = labelArray == 0 & AHICutIdx == 1;
NRBreathFeat = BreathFLDataTableFinal(NRBreathIdx,featNum);

nbins = 25;

%Responder distribution for all subjects
figure
r = histogram(RBreathFeat{:,:}, nbins);
%format 
title('Responders')
xlabel('Flow')
ylabel('Count')
set(gca,'Box','On','FontSize',15)
% xlim([0 0.02])
% ylim([0 2000])

%Non-responder
figure
nr = histogram(NRBreathFeat{:,:}, nbins);
%format 
title('Non-responders')
xlabel('Flow')
ylabel('Count')
set(gca,'Box','On','FontSize',15)
% xlim([0 0.02])
% ylim([0 2000])

saveas(r, [figfolder,'\RHist'], 'png')
saveas(nr, [figfolder, '\NRHist'], 'png')



% Subject-level plots
f1 = figure('Name', 'Responders', 'pos', [-1644 41 1613 913]);
f2 = figure('Name', 'Non-Responders', 'pos', [26 49 1613 913]);
rcount = 0;
nrcount = 0;
subjects = unique(BreathDataTableFinal.Subject);

for subnum = 1:length(subjects)
    subIdx = BreathDataTableFinal.Subject == subjects(subnum);
    
    if sum(labelArray(subIdx)) > 0
        %Responder breaths
        RBreathIdx = BreathDataTableFinal.Subject == subjects(subnum) &...
        labelArray == 1 & AHICutIdx == 1;
        RBreathFeat = BreathFLDataTableFinal(RBreathIdx, featNum);
        
        if isempty(RBreathFeat)
            continue
        end
        
        %Plot histogran
        rcount = rcount + 1;
        figure(f1)
        subplot(3,4,rcount)
        r = histogram(RBreathFeat{:,:}, nbins);
        
        %format 
        title(subjects(subnum))
        xlabel('Flow')
        ylabel('Count')
        set(gca,'Box','On','FontSize',15)
%         xlim([0 0.02])
        
%         saveas(r, [figfolder, '\', num2str(subjects(subnum)),'Rhist'], 'png')
    else
        %Non-responder breaths
        NRBreathIdx = BreathDataTableFinal.Subject == subjects(subnum)  &...
            labelArray == 0 & AHICutIdx == 1;
        NRBreathFeat = BreathFLDataTableFinal(NRBreathIdx, featNum);
        
        if isempty(NRBreathFeat)
            continue
        end
        
        %Plot histogran
        nrcount = nrcount + 1;
        figure(f2)
        subplot(3,4,nrcount)
        nr = histogram(NRBreathFeat{:,:}, nbins);
        
        %format 
        title(subjects(subnum))
        xlabel('Flow')
        ylabel('Count')
        set(gca,'Box','On','FontSize',15)
%         xlim([0 0.02])
        
%         saveas(nr, [figfolder, '\', num2str(subjects(subnum)),'NRhist'], 'png')
    end
end

saveas(r, [figfolder, '\', 'Rhist_sub'], 'png')
saveas(nr, [figfolder, '\', 'NRhist_sub'], 'png')