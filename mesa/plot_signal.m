addpath('/Users/sq566/Desktop/pats/edfs/PUPbeta/External/');
id = "mesa-sleep-0010";

original = "/Users/sq566/Desktop/mesa-5/edfs/" + id + ".edf";
[header, signalHeader, signalCell] = blockEdfLoad(char(original));

% Get signal Index position
targetString = 'DHR';
signalLabels = {signalHeader.signal_labels};
index = find(strcmp(signalLabels, targetString));

physioData = signalCell{1,index};
sr = signalHeader(index).samples_in_record;
timeSec = 1:length(physioData);

% Now plot the data
fig = figure;
set(fig, 'Position', [100, 100, 5000, 500]);
plot(timeSec, physioData, 'b-', 'LineWidth', 2);
xlabel('Time (seconds)');
title(targetString);
grid on

annotation('textbox', [0.5, 0.9, 0.1, 0.1], 'String', id, ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold');

print(fig, '/Users/sq566/Desktop/mesa-5/plots/' + id + targetString +'.png', '-dpng', '-r300');




