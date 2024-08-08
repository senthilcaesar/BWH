addpath('/Users/sq566/Desktop/pats/edfs/PUPbeta/External/');
id = "pats-813922-baseline";
eannot = "/Users/sq566/Desktop/pats/inspect/eannot/" + id + ".eannot";
fileID = fopen(eannot, 'r');
labels = textscan(fileID, '%s');

original = "/Users/sq566/Desktop/pats/inspect/" + id + ".edf";
[header, signalHeader, signalCell] = blockEdfLoad(char(original));

% Get SpO2 signal Index position
targetString = 'Pleth';
signalLabels = {signalHeader.signal_labels};
index = find(strcmp(signalLabels, targetString));

physioData = signalCell{1,index};
sr = signalHeader(index).samples_in_record;

epochDuration = 30;
timeSec = 1:length(physioData);
labelMap = containers.Map({'L', 'W', 'N1', 'N2', 'N3', 'R'}, [0, 1, 2, 3, 4, 5]);
numericStages = cellfun(@(x) labelMap(x), labels{1});
expandedLabels = repelem(numericStages, epochDuration * sr); % Expand labels to match the resolution of physiological data

% Now plot the data
fig = figure;
set(fig, 'Position', [100, 100, 8000, 600]); 
subplot(2, 1, 1);
stairs(timeSec, expandedLabels, 'LineWidth', 2);
yticks(0:5);
yticklabels({'L', 'W', 'N1', 'N2', 'N3', 'R'});
ylabel('Sleep Stage');
xlabel('Time (seconds)');
title('Hypnogram');
grid on;

subplot(2, 1, 2);
plot(timeSec, physioData, 'b-', 'LineWidth', 2);
ylabel('SpO2 Percent');
xlabel('Time (seconds)');
title('SpO2 Signal');
grid on

annotation('textbox', [0.5, 0.9, 0.1, 0.1], 'String', id, ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 24, 'FontWeight', 'bold');

print(fig, '/Users/sq566/Desktop/pats/inspect/' + id + '.png', '-dpng', '-r300');
