% Luna PLACE analysis

annot_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data/annot/scorer_1';
edf_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data/edf';
db_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data/annot/scorer_1';
subs = readlines([annot_dir '/place_sub.txt'], "EmptyLineRule","skip");

for idx=1:length(subs)
    disp(['Processing subject ' num2str(idx)])
    subname = subs{idx};

    % Luna SOAP command
    cmd = ['/usr/local/bin/luna ' edf_dir '/' subname '-edit.edf' ' alias="C4_M1|C3_M2" -s PLACE stages=' annot_dir '/' subname '.eannot'];

    % Run luna command in bash shell
    [status,cmdout] = system(cmd);

    disp(cmdout)
end

