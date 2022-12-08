%% 
%first, set current directory as source file directory

Dir = uigetdir();

%then run this code
x=dir(Dir); clear temp; for i=3:length(x), Filename{i-2,1} = x(i).name; end
% 
% if everysecondfilename
%     temp = temp(2:2:end);
% end


clear Name Extension
for i=1:length(Filename)
    I = find(Filename{i}=='.',1,'first');
    Name{i,1} = Filename{i}(1:I-1);
    Extension{i,1} = Filename{i}(I+1:end);
end

T = table(Filename,Name,Extension);

PDir = Dir(1:find(Dir=='\',1,'last'));
DirName = Dir(find(Dir=='\',1,'last')+1:end);

writetable(T,[PDir DirName 'Files.csv']);

