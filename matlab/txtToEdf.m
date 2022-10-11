%{

This program read signals from ASCII-formatted, tab-delimited plain-text files
and writes compressed EDF file

%}

fs = 250;
data_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data';
subs = readlines([data_dir '/txt/sub.txt'], "EmptyLineRule","skip");
edf_dir = [data_dir '/edf'];
elect = ['C3_M2,F3_F4,F3_M2,F3_O1,F4_M1,F4_O2,FP1_F3,FP1_M2,FP1_O1,' ...
         'FP2_F4,FP2_M1,FP2_O2,ECG,EMG,EOG1,EOG2'];

for idx=1:length(subs)
    disp(['Processing subject ' num2str(idx)])
    subname = subs{idx};
    filename = [data_dir '/txt/' subname];

    % Luna command to convert a text file into an EDF
    luna_cmd = ['/usr/local/bin/luna ' filename ' --fs=' num2str(fs) ' --chs=' elect... 
                ' -s WRITE edf edf-dir=' edf_dir];

    % Run luna command in bash shell
    system(luna_cmd);

    out_file = erase(subname,".txt");
    movefile([data_dir '/edf/' subname '.edf'], [data_dir '/edf/' out_file '.edf']);

    % Set EDF header units paramater to Millivolts
    set_header = ['/usr/local/bin/luna ' data_dir '/edf/' out_file '.edf -s "SET-HEADERS unit=mV & WRITE edf-tag=edit"'];
    system(set_header)

end