function [ArousalIntensity,ArScorecontinuous2]=ArousalIntensityRun(EEGsignals,ArousalSig2,SleepStages2,Time,TimeEEG,WakeSleepInfoAUC_M)
warning off;
%% make ArousalSig2 and SleepStages2 from ArousalSig,SleepStages,Time,TimeEEG
if length(Time)~=length(TimeEEG)
    ArousalSig = interp1(Time,ArousalSig2,TimeEEG,'nearest','extrap');
    SleepStages = interp1(Time,SleepStages2,TimeEEG,'nearest','extrap');
else
    ArousalSig=ArousalSig2;
    SleepStages=SleepStages2;
end
%%

Fs=1/(Time(2)-Time(1));

AnalyzeBothEEGs=0;
addpath([pwd '\ArousalIntensity']);
TrainMtx=xlsread('TrainingSet.xlsx');

ArousalSegs=FindSegs(ArousalSig);

%%%removing arousals starting in wakefullness and arousals with less then
%%%10 seconds apart
ArInWakeIdx=find(SleepStages(ArousalSegs(1,:))>3);
ArSegsInWake=ArousalSegs(:,ArInWakeIdx);
st0=1;
IdxToRemove=[];
for ii=1:length(ArInWakeIdx)
   if mean(SleepStages(st0:ArSegsInWake(1,ii)))>3 
       IdxToRemove=[IdxToRemove;ArInWakeIdx(ii)];
   end
   st0=ArSegsInWake(2,ii); 
end

ArousalSegs(:,IdxToRemove)=[];
TooCloseArousals=find(ArousalSegs(1,2:end)-ArousalSegs(2,1:end-1)<=10*Fs);
ArousalSegs(:,TooCloseArousals+1)=[];

%%%removing arousals longer than 15 seconds
ArousalSegs(2,ArousalSegs(2,:)-ArousalSegs(1,:)>15*Fs)=ArousalSegs(1,ArousalSegs(2,:)-ArousalSegs(1,:)>15*Fs)+15*Fs;

ArStrts=round(ArousalSegs(1,:));
ArEnds=round(ArousalSegs(2,:));

ArousalIntensity=FindVecs(ArousalSegs,length(ArousalSig));

Nsignals=length(EEGsignals);
n=1;

for i=1:Nsignals
    try
        signal{n}=evalin('caller',EEGsignals{i}); %isempty(eval(signallist{i}))
        eval([EEGsignals{i} '=signal{n};']);
        n=n+1;
    catch me
        disp(me.message);
        EEGsignals{i}=[];
    end
end

ArStruct.training_set =TrainMtx;

[~,BestEEGs]=sort(WakeSleepInfoAUC_M,'descend');
C3_A2=signal{BestEEGs(1)};
C4_A1=signal{BestEEGs(2)};

for ii=1:length(ArStrts)
    ArDur=ArEnds(ii)-ArStrts(ii)+1;
    if isempty(find(SleepStages(ArStrts(ii)-ArDur:ArStrts(ii)-1)==0,1))
        ArStruct.sig=C3_A2(ArStrts(ii):ArEnds(ii));
        ArStruct.presig=C3_A2(ArStrts(ii)-ArDur:ArStrts(ii)-1);
        ArScore1=ScoreArousal(ArStruct,Fs);
        
        ArStruct.sig=C4_A1(ArStrts(ii):ArEnds(ii));
        ArStruct.presig=C4_A1(ArStrts(ii)-ArDur:ArStrts(ii)-1);
        ArScore2=ScoreArousal(ArStruct,Fs);
        
        ArScore(ii)=max(ArScore1,ArScore2);
        ArousalIntensity(ArStrts(ii):ArEnds(ii))=ArScore(ii);
    else
        ArScore(ii)=NaN;
    end
end

decomLevel=5; %%% Wavelet decomposition level
wavefun='db4'; %%% Wavelet function

WinSize=round(6*Fs);
Ovlap=round(5*Fs);
C3_A2_Mat = buffer(C3_A2,WinSize,Ovlap);

Time_Mat = buffer(TimeEEG+1,WinSize,Ovlap);
Time_Mat(Time_Mat==0)=NaN;
Time_Mat=Time_Mat-1;


if AnalyzeBothEEGs
    C4_A1_Mat = buffer(C4_A1,WinSize,Ovlap);
end

SleepStages_Mat = buffer(SleepStages,WinSize,Ovlap);
Sleep_Mat=SleepStages_Mat<4;

ArousalIntensity_Mat = buffer(ArousalIntensity,WinSize,Ovlap);
NoArousal_Mat=ArousalIntensity_Mat<1;

SleepNoArousalMat=sum(NoArousal_Mat&Sleep_Mat)==WinSize;
SleepNoArousalCSum=cumsum(SleepNoArousalMat);

for ii=1:size(C3_A2_Mat,2)
    [C3_A2_ftrs(:,ii)]=FindWaveFtrs(C3_A2_Mat(:,ii),decomLevel,wavefun,Fs);
    
    if AnalyzeBothEEGs
        [C4_A1_ftrs(:,ii)]=FindWaveFtrs(C4_A1_Mat(:,ii),decomLevel,wavefun,Fs);
    end
end


[ClassificationModels]=TrainArousalIntensity(TrainMtx);

winlen=45;
StartPoint=find(SleepNoArousalCSum==winlen);
StartPoint=StartPoint(1);
EndPoint=find(SleepNoArousalCSum==SleepNoArousalCSum(end)-winlen);
EndPoint=EndPoint(end);
ArScorecontinuous=nan(size(SleepNoArousalCSum));

f = waitbar(0,'','Name','Arousal Intensity...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for ii=StartPoint+1:EndPoint-1
    StIdx=find(SleepNoArousalCSum==SleepNoArousalCSum(ii)-winlen);
    EndIdx=find(SleepNoArousalCSum==SleepNoArousalCSum(ii)+winlen);
    if isempty(StIdx)
        StIdx=1;
    end
    if isempty(EndIdx)
        EndIdx=size(C3_A2_ftrs,2);
    end
    LocalWin=StIdx:EndIdx;
        
    C3_A2_baseftrs=prctile(C3_A2_ftrs(:,LocalWin(SleepNoArousalMat(LocalWin)==1)),25,2)';
    C3_A2_ar_ftrs=C3_A2_ftrs(:,ii)';
    C3_A2_Score=TestArousalIntensity(C3_A2_ar_ftrs,C3_A2_baseftrs,ClassificationModels);
    
    if AnalyzeBothEEGs
        C4_A1_baseftrs=prctile(C4_A1_ftrs(:,LocalWin(SleepNoArousalMat(LocalWin)==1)),25,2)';
        C4_A1_ar_ftrs=C4_A1_ftrs(:,ii)';
        C4_A1_Score=TestArousalIntensity(C4_A1_ar_ftrs,C4_A1_baseftrs,ClassificationModels);
        ArScorecontinuous(ii)=max(C3_A2_Score,C4_A1_Score);
    else
        ArScorecontinuous(ii)=C3_A2_Score;
    end
    if getappdata(f,'canceling')
        break
    end
    waitbar(ii/(EndPoint-1),f)
end

Time_ArCnt=nanmean(Time_Mat);
ArScorecontinuous2 = interp1(Time_ArCnt,ArScorecontinuous,Time,'linear');

if length(Time)~=length(TimeEEG)
    ArousalIntensity = interp1(TimeEEG,ArousalIntensity,Time,'nearest');
end

delete(f)
