%{

This program reads H5 data file and writes the signal data to a text file
The text files can be converted to EDF format using luna module.

Output text file format (ASCII-formatted, tab-delimited plain-text files)
------------------------------------------------------------------------
signal1    signal2    signal3  ................... signalN
  .           .
  .           .
  .           .
  .           .
  .           .
  .           .

%}

data_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data';
h5_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data/h5';
subs = readlines([data_dir '/h5/test.txt'], "EmptyLineRule","skip");
hypno_map = containers.Map({'9','0','1','2','3','4'}, {'NS', 'W', 'N1', 'N2', 'N3', 'R'});
eeg_chans = {'C3_M2', 'F3_F4', 'F3_M2', 'F3_O1', 'F4_M1', 'F4_O2', 'FP1_F3', 'FP1_M2', 'FP1_O1',... 
             'FP2_F4', 'FP2_M1', 'FP2_O2'};
emg_chans = {'ECG', 'EMG'};
eog_chans = {'EOG1', 'EOG2'};

%h5disp(filename, '/signals');

for idx=1:length(subs)
    disp(['Processing subject ' num2str(idx)])
    subname = subs{idx};
    filename = [h5_dir '/' subname];
    fs = h5readatt(filename,'/signals/eeg/', 'fs');
    unit = h5readatt(filename,'/signals/eeg/', 'unit');
    hypno_data = h5read(filename,'/hypnogram');
    hypno_data(hypno_data == -1) = 9;
    hypno_val = values(hypno_map, cellstr(num2str(hypno_data)));
    record_size = h5info(filename, '/signals/eeg/C3_M2').Dataspace.Size;
    electrodes_data = zeros(record_size,16); % 16 electrodes

    % EEG data
    for i = 1:length(eeg_chans)
        eeg_data = h5read(filename, ['/signals/eeg/' eeg_chans{i}]);
        electrodes_data(:,i) = eeg_data;
        clear eeg_data;
    end
    % EMG data
    for j = 1:length(emg_chans)
        emg_data = h5read(filename, ['/signals/emg/' emg_chans{j}]);
        electrodes_data(:,j+12) = emg_data;
        clear emg_data;
    end
    % EOG data
    for k = 1:length(eog_chans)
        eog_data = h5read(filename, ['/signals/eog/' eog_chans{k}]);
        electrodes_data(:,k+14) = eog_data;
        clear eog_data;
    end
    
    % Write data to text file, each column is a signal
    out_file = erase(subname,".h5");
    writematrix(electrodes_data,[data_dir '/txt/' out_file '.txt'],'Delimiter','tab')
    writecell(hypno_val,[data_dir '/annot/' out_file '._eannot.txt'])
    movefile([data_dir '/annot/' out_file '._eannot.txt'], [data_dir '/annot/' out_file '.eannot']);
    clear electrodes_data;
    
end


