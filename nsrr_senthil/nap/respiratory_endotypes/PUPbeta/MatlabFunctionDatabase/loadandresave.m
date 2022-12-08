function loadandresave(t_filename,t_backup,t_instruction)
disp('started loading');
load(t_filename);
disp('finished loading');
t_i=strfind(t_filename,'.');
t_tempfile = [t_filename(1:t_i-1) '_temp' t_filename(t_i:end)];
t_status = movefile(t_filename,t_tempfile);
eval(t_instruction);
disp('finished saving backup');
save(t_filename, '-regexp', '^(?!t_.*$).','-v7.3');
disp('finished resaving');
if ~t_backup
   delete(t_tempfile); 
   disp('deleted backup');
end


