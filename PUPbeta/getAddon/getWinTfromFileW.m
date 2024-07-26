function WinT = getWinTfromFileW(W)

%legacy handling, does this: DataOut = DataOut{1} ... etc
datatoloadlist = fieldnames(W);
for i=1:length(datatoloadlist)
    if isfield(W,datatoloadlist{i}) && iscell(W.(datatoloadlist{i})) && size(W.(datatoloadlist{i}),1)==1 ...
            && size(W.(datatoloadlist{i}),2)==1  %exist, be a 1x1 cell, contains a cell
        W.(datatoloadlist{i}) = W.(datatoloadlist{i}){1}; 
    end
end      

WinT=table;

AnalysisIndexT = array2table(W.AnalysisIndex);
AnalysisIndexT.Properties.VariableNames = {'AnalyzeIndL','AnalyzeIndR'};
WinT = [WinT AnalysisIndexT];


SleepDataT = array2table(W.SleepData);
SleepDataT.Properties.VariableNames = {'FWake','FNREM','FNREM1','FNREM2','FNREM3','FREM','LongestWake'};
WinT = [WinT SleepDataT];

LGplusinfoT = array2table(W.LGplusinfo);
LGplusinfoT.Properties.VariableNames = {'TimeB1','LG0','tau','tau2','LGn','Tn','LG1','LG2','delay','VRA','VRA2','ArThres','MSE','TtotMean','TtotStd','TtotMedian','TtotIQR'};
WinT = [WinT LGplusinfoT];

LGQualityInfoT = array2table(W.LG_QualityInfo(:,1:3)); %only use first 3 cols
LGQualityInfoT.Properties.VariableNames = {'Narousals','Nevents','meanNotE'};
WinT = [WinT LGQualityInfoT];

clear OneMinusRsq
for jj=1:size(W.fitQual,2)
    temp3pt1=W.fitQual{jj};
    if length(temp3pt1)>1
        OneMinusRsq(jj,:)=1-temp3pt1(2); %Element #2 of FitQual
    else
        OneMinusRsq(jj,:)=NaN;
    end
end
WinT.OneMinusRsq = OneMinusRsq;

StoNDataT = array2table(W.StoNData.Fnoise);
StoNDataT.Properties.VariableNames = {'FNoise1','FNoise2','FNoise3'};
WinT = [WinT StoNDataT];


CPAPDataT = array2table(W.CPAPData);
CPAPDataT.Properties.VariableNames = {'CPAPoff','CPAPmedian','CPAPstd','CPAPabs95'};
WinT = [WinT CPAPDataT]; 

