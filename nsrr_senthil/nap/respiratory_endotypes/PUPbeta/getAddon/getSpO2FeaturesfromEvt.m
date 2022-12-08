function getSpO2FeaturesfromEvt(MrangeOverride)

global AMasterSpreadsheet settings

t_startGetData = clock;

%% read excel spreadsheet

[~,txt] = xlsread(AMasterSpreadsheet,1,'AD4:AE6000');
M = size(txt,1);

[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');

Study = settings.patients(:,1);


ID=1:M;

% override if mrange present
if exist('MrangeOverride')
    ID = MrangeOverride(:);
end


blank = nan(length(ID),15);
SpO2features =array2table(blank);

for i=1:length(ID)
    try
        clear W Evts
        
        AnalyzedDir = [path{:} settings.savename '_' num2str(i) '.mat'];
        %         W = load([txt{ID(i),2} txt{ID(i),1}]);
        W=load(AnalyzedDir,'Evts');
        disp(['processing: ' txt{ID(i),1}])
        Evts=W.Evts;
        if isfield(Evts,'SpO2')
             
            SpO2features(i,:)=Evts.SpO2;
            if i==1
                SpO2features.Properties.VariableNames=Evts.SpO2.Properties.VariableNames;
            end
            
        else
            disp('run getEvtsAddOn before proceeding');
            disp ('failed to extract spo2')
        end
        
    catch
        disp ('failed to extract spo2')
    end
end
StudyId=table(Study);
OvernightSpO2featuresT=[StudyId SpO2features];
save([settings.workdir,'Summary' filesep 'OvernightSpO2features.mat'],'OvernightSpO2featuresT','-v7.3');

