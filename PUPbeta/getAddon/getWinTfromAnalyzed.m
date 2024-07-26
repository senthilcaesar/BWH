function getWinTfromAnalyzed(MrangeOverride)
% RUN StartHere.m first

% this function creates WinT from analyzed files.

global settings AMasterSpreadsheet

t_startGetData = clock;


settings = ImportSettings(settings,AMasterSpreadsheet);

%% LOAD AMASTERSPREADSHEET
[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,1),'UniformOutput',false)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
AnalyzeList=cell2mat(MasterWorksheet(:,5));
M = size(AnalyzeList,1); % using all analyzed files.

ID=1:M;

% override if mrange present
if exist('MrangeOverride')
    ID = MrangeOverride(:);
end

success=zeros(length(ID),1);

for i=2%1:length(ID)
    try
        AnalyzedDir = [settings.AnalyzedDirectory settings.savename '_' num2str(i) '.mat'];
        
        if exist(AnalyzedDir)==2
            disp(['Processing: ' num2str(i) ': ' settings.savename '_' num2str(i) '.mat']);
            W = load(AnalyzedDir);
            
            %% simplifies cell array for each variable
            %e.g. CPAPData = CPAPData{1} if needed.
            datatoloadlist = fieldnames(W);
            for jj=1:length(datatoloadlist)
                if isfield(W,datatoloadlist{jj}) && iscell(W.(datatoloadlist{jj})) && size(W.(datatoloadlist{jj}),1)==1 ...
                        && size(W.(datatoloadlist{jj}),2)==1  %exist, be a 1x1 cell, contains a cell
                    W.(datatoloadlist{jj}) = W.(datatoloadlist{jj}){1};
                end
            end
            
            %% convert to tables
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
            
                    
            % decided to remove EventsInfo.so not adding that in WinT
            
            W.WinT = WinT;
            save(AnalyzedDir,'-struct','W','-v7.3');
            success(i)=1;
            disp(['success: WinT created for: ' settings.savename '_' num2str(i) '.mat']);
            
            % to do ---remove the fields in original data once we decide
            % this looks okay.
            
        end
    end
end



%%





%%









