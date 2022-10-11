% Luna SOAP analysis

annot_dir = '/Users/sq566/Desktop/scorer_5';
db_dir = '/Users/sq566/Desktop/scorer_5';
edf_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data/edf';
subs = readlines([annot_dir '/sub_name.txt'], "EmptyLineRule","skip");
pat = digitsPattern(1) + "-level classification: kappa = " + digitsPattern(1) + "." + digitsPattern(2);

kappaTable(25) = struct('Kappa5',[], 'Kappa3',[]);

for idx=1:length(subs)
    disp(['Processing subject ' num2str(idx)])
    subname = subs{idx};

    % Luna SOAP command
    cmd = ['/usr/local/bin/luna ' edf_dir '/' subname '-edit.edf -o ' db_dir '/' subname '.db annot-file=' annot_dir '/' subname '.eannot "alias=C4_M1|C3_M2" -s SOAP epoch'];

    % Run luna command in bash shell
    [status,cmdout] = system(cmd);

    kappa = extract(cmdout,pat);
    disp(subname)
    disp(kappa)

    % Extract the kappa float value from string
    kap5 = regexp(kappa{1},'\d+\.?\d*','match');
    kap3 = regexp(kappa{2},'\d+\.?\d*','match');

    kappaTable(idx).Kappa5 = str2double(kap5{2});
    kappaTable(idx).Kappa3 = str2double(kap3{2});
    clear kappa subname cmd;
end



