function T = fixNoxAudioLabel()
global settings

T=table();
M = length(settings.patients(:,1));
blankcol = zeros(M,1);
T.Index = [1:M]';
T.HaveSource=blankcol;
SourceDir=strcat(string(settings.Filenames(:,4)),string(settings.Filenames(:,1)));
SourceDirE=strcat(string(settings.Filenames(:,5)),string(settings.Filenames(:,2)));
SourceDirH=strcat(string(settings.Filenames(:,6)),string(settings.Filenames(:,3)));
ConvertedFileDir = string(strcat(settings.ConvertedDirectory,settings.patients(:,1),'.mat'));
AnalyzedFileDir = string(strcat(settings.AnalyzedDirectory,settings.savename,'_',string([1:M]'),'.mat'));


for i=1:M
    clear w
    try
        w = load(ConvertedFileDir(i));
        ind = find(w.SigT.Properties.VariableNames=="AudioVolume");
        if length(ind)==1
            disp('Found AudioVolume and replacing name with NoxAudio');
            w.SigT.Properties.VariableNames{ind}='NoxAudio';
            disp(strcat("saving ",ConvertedFileDir(i)));
            save(ConvertedFileDir(i),'-struct','w','-v7.3');
            disp('done');
        end
    end
end