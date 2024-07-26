

default = '\\Cifs2\dscd_j$\PEOPLE\FACULTY\SANDS\O2PSG\TraitsAnalysis\Source';
[file,folder]=uigetfile(default)
%%
t_filename = [folder file];
%instruction='clear Flow';
instruction='';
savebackup=1;
loadandresave(t_filename,savebackup,instruction);