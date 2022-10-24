%{
Also, could you also look at the rater-by-rater agreements in staging?   
(vs. each other, and vs. the consensus).  And summarize across EDFs?   
(Again, we may need to confirm that “rater 2” for one EDF is actually “rater 2” for the others, etc.
 
You can use a Luna convenience function to do this -  normal the “stage comparison” code 
(i.e. to calculate accuracy, kappa, etc) runs against the “internal” staging and whatever SOAP or POPS has predicted,.    
You can instead co-opt this by manually specifying two sets and the EVAL-STAGES command.    (i.e. this won’t even look at the signal data).
%}

base_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data';
luna_binDir = '/Users/sq566/Downloads/luna-base-0.27';
annot_dir1 = [ base_dir '/annot/scorer_4'];
annot_dir2 = [ base_dir '/annot/consensus'];
db_dir = [ base_dir '/annot/db'];
edf_dir = [ base_dir '/edf'];
subs = readlines([base_dir '/sub.txt'], "EmptyLineRule","skip");
kappaTable(25) = struct('Kappa3',[]);

for idx=1:length(subs)
    disp(['Processing subject ' num2str(idx)])
    subname = subs{idx};

    % Luna SOAP command
    cmd = [ luna_binDir '/luna ' edf_dir '/' subname '-edit.edf -o ' db_dir '/' subname '.db annot-file=' annot_dir1 '/' subname '.eannot -s EVAL-STAGES file=' annot_dir2 '/' subname '.eannot'];
    
    % Run luna command in bash shell
    [status,cmdout] = system(cmd);

    start_idx = strfind(cmdout, '3-class kappa');
    end_idx = strfind(cmdout, '(n =');
    cmdout = cmdout(start_idx:end_idx-2);

    % Extract the kappa float value from string
    kappa_loc = strfind(cmdout, '=');
    kappa_coeff = cmdout(kappa_loc+2:end);

    disp(subname)
    disp(kappa_coeff)

    kappaTable(idx).Kappa3 = str2double(kappa_coeff);
    clear kappa subname cmd kappa_coeff kappa_loc start_idx end_idx;
end




