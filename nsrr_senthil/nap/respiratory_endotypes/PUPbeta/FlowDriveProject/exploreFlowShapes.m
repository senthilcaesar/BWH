function exploreFlowShapes(BreathDataTableFinal, InspArrayFinal, labelArray, AHICutIdx, ylimit)

figfolder = 'J:\PEOPLE\POST DOC\VENA\Flow shape\Results and Processed Data\Figures\PSG MM Scope';

% Plot responders and non-responders
%Plot all subjects in R and NR categories
%Responder breaths
RbreathIdx = labelArray == 1 & AHICutIdx == 1;
Rbreath = InspArrayFinal(RbreathIdx, :);

%Responder breaths
figure('pos', [722 602 464 285]);
for breathNum = 1:size(Rbreath,1)
    r = plot(Rbreath(breathNum,:), 'k-');
    r.Color(4) = 0.01;
    hold on
end

%plot mean curve
RbreathMean = nanmean(Rbreath,1);
r = plot(RbreathMean, 'k-', 'LineWidth', 2);

%format 
xlabel('Index')
ylabel('Flow (unitless)')
set(gca,'Box','On','FontSize',15)
ylim(ylimit)
xlim([0 size(InspArrayFinal,2)])


%Non-responder breaths
NRbreathIdx = labelArray == 0 & AHICutIdx == 1;
NRbreath = InspArrayFinal(NRbreathIdx, :);

figure('pos', [722 602 464 285]);
for breathNum = 1:size(NRbreath,1)
    nr = plot(NRbreath(breathNum,:), 'k-');
    nr.Color(4) = 0.01;
    hold on
end

%plot mean curve
NRbreathMean = nanmean(NRbreath,1);
nr = plot(NRbreathMean, 'k-', 'LineWidth', 2);

%format 
xlabel('Index')
ylabel('Flow (unitless)')
set(gca,'Box','On','FontSize',15)
ylim(ylimit)
xlim([0 size(InspArrayFinal,2)])

saveas(r, [figfolder,'\Rbreaths'], 'png')
saveas(nr, [figfolder, '\NRbreaths'], 'png')

% Subject-level plots
subjects = unique(BreathDataTableFinal.Subject);
for subnum = 1:length(subjects)
    subIdx = BreathDataTableFinal.Subject == subjects(subnum);
    
    if sum(labelArray(subIdx)) > 0
        %Responder breaths
        RbreathIdx = BreathDataTableFinal.Subject == subjects(subnum) &...
        labelArray == 1 & AHICutIdx == 1;
        Rbreath = InspArrayFinal(RbreathIdx, :);
        
        if isempty(Rbreath)
            continue
        end
        
        %Responder breaths
        figure('pos', [722 673 368 230]);
        for breathNum = 1:size(Rbreath,1)
            r = plot(Rbreath(breathNum,:), 'k-');
            r.Color(4) = 0.05;
            hold on
        end

        %plot mean curve
        RbreathMean = nanmean(Rbreath,1);
        r = plot(RbreathMean, 'k-', 'LineWidth', 2);

        %format 
        title(subjects(subnum))
        xlabel('Index')
        ylabel('Flow (unitless)')
        set(gca,'Box','On','FontSize',15)
        ylim(ylimit)
        xlim([0 size(InspArrayFinal,2)])
        
        saveas(r, [figfolder, '\', num2str(subjects(subnum)),'R'], 'png')
    else
        %Non-responder breaths
        NRbreathIdx = BreathDataTableFinal.Subject == subjects(subnum)  &...
            labelArray == 0 & AHICutIdx == 1;
        NRbreath = InspArrayFinal(NRbreathIdx, :);
        
        if isempty(NRbreath)
            continue
        end
        
        figure('pos', [722 673 368 230]);
        for breathNum = 1:size(NRbreath,1)
            nr = plot(NRbreath(breathNum,:), 'k-');
            nr.Color(4) = 0.05;
            hold on
        end

        %plot mean curve
        NRbreathMean = nanmean(NRbreath,1);
        nr = plot(NRbreathMean, 'k-', 'LineWidth', 2);

        %format 
        title(subjects(subnum))
        xlabel('Index')
        ylabel('Flow (unitless)')
        set(gca,'Box','On','FontSize',15)
        ylim(ylimit)
        xlim([0 size(InspArrayFinal,2)])
        
        saveas(nr, [figfolder, '\', num2str(subjects(subnum)),'NR'], 'png')
    end
end