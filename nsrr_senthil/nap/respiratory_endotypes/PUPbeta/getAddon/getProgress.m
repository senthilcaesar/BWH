function T = getProgress(WriteReady,WriteHaveSource,OldNeedsRerun,ConvertReady,AnalyzeReady)
global settings AMasterSpreadsheet
if ~exist('WriteReady','var')
    WriteReady=0;
end
if ~exist('WriteHaveSource','var')
    WriteHaveSource=0;
end
if ~exist('ConvertReady','var')
    ConvertReady=0;
end
if ~exist('AnalyzeReady','var')
    AnalyzeReady=0;
end
if ~exist('OldNeedsRerun','var')
    OldNeedsRerun=1;
end

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

T.HaveSource=blankcol;
T.HaveSourceE=blankcol;
T.HaveSourceH=blankcol;
T.HaveConverted=blankcol;
T.HaveAnalyzed=blankcol;
T.SourceDate=strings(M,1);
T.ConvertedDate=strings(M,1);
T.AnalyzedDate=strings(M,1);

for i=1:M
    T.HaveSource(i)=exist(SourceDir(i),'file')==2;
    T.HaveSourceE(i)=exist(SourceDirE(i),'file')==2;
    T.HaveSourceH(i)=exist(SourceDirH(i),'file')==2;
    T.HaveConverted(i)=exist(ConvertedFileDir(i),'file')==2;
    T.HaveAnalyzed(i)=exist(AnalyzedFileDir(i),'file')==2;
    if T.HaveConverted(i)==1
        temp = dir(SourceDir(i));
        T.SourceDate(i) = string(temp.date);
    end
    if T.HaveConverted(i)==1
        temp = dir(ConvertedFileDir(i));
        T.ConvertedDate(i) = string(temp.date);
    end
    if T.HaveAnalyzed(i)==1
        temp = dir(AnalyzedFileDir(i));
        T.AnalyzedDate(i) = string(temp.date);
    end
end

T.NeedSource=~(T.HaveSource & T.HaveSourceE & T.HaveSourceH);
T.NeedConverted=~T.HaveConverted;
T.NeedAnalyzed=~T.HaveAnalyzed;

INeedSource = find(T.NeedSource)';
INeedSource = find(T.NeedSource)';
MissingSourceFilesFor = settings.Filenames(INeedSource,1)

T.HaveConvertedOld = datetime(T.SourceDate)>datetime(T.ConvertedDate);
T.HaveConvertedOld = T.HaveConvertedOld | datetime(T.ConvertedDate)<datetime("06-Feb-2022 17:07:21") & (startsWith(settings.patients(:,1),"P")==0);
if OldNeedsRerun
    T.ReadyToConvert = (T.HaveSource & T.NeedConverted) | T.HaveConvertedOld;
else
    T.ReadyToConvert = T.HaveSource & T.NeedConverted;
end
IReadyToConvert = find(T.ReadyToConvert)'
ReadyToConvertF = settings.Filenames(T.ReadyToConvert,1)

T.HaveAnalyzedOld = datetime(T.ConvertedDate)>datetime(T.AnalyzedDate);
if OldNeedsRerun
    T.ReadyToAnalyze = (T.HaveConverted & T.NeedAnalyzed) | T.HaveAnalyzedOld;
else
    T.ReadyToAnalyze = (T.HaveConverted & T.NeedAnalyzed);
end
IReadyToAnalyze = find(T.ReadyToAnalyze)'

if WriteReady
    xlswrite(AMasterSpreadsheet,double(T.ReadyToConvert),1,'AB4');
    xlswrite(AMasterSpreadsheet,double(T.ReadyToAnalyze),1,'AH4');
    %run code to open file also
end
if WriteHaveSource
    xlswrite(AMasterSpreadsheet,double(T.HaveSource),1,'Q4'); %careful not to replace data in col Q
    %run code to open file also
end
if ConvertReady
   ConvertN(IReadyToConvert);
end
if AnalyzeReady
   AnalyzeN(IReadyToAnalyze);
end

