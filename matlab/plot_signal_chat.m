subjects = ListSubjects();
electrodes = ListElect();
DATA_DIR = '/Users/sq566/Downloads/nsrr_studies';
study = 'chat';

for idx=1:9
    display(['Processing interation ' num2str(idx)])

    subID = idx;
    chanID = randi([4 7]);
    
    subIS = subjects{subID};
    elecIS = electrodes{chanID};
    file = [DATA_DIR '/' study '/data/edf/' subIS '.edf'];
    
    [hdr, record] = edfread(file);
    fs = hdr.samples(chanID);
    y = size(record,2); %randi([1 size(record,2)]);
    x = y - 60001;
    
    % Randomly plot 5 mins of signal
    f = figure('visible','off');
    set(gcf,'Position',[300 300 1400 500])
    plot(record(chanID,x:y))
    xlim([0 (y-x)]);
    xticklabels((xticks+x)/fs)
    title(['Subject: ' subIS ', Channel: ' elecIS], 'FontSize',14);
    fname = [DATA_DIR '/' study '/data/plots/' subIS '_' elecIS '.png'];
    print(f,fname,'-r600','-dpng');
    %saveas(f,fname);  
    close(f);
    
    clear hdr record x y subID elecID subID elecIS;
end