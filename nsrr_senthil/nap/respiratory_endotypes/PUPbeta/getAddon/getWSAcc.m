function WSInfoT = getWSAcc(MrangeOverride)
global AMasterSpreadsheet settings

Study = settings.patients(:,1);
M = size(Study,1);
Mrange = 1:558; %az added code

if ~isfield(settings,'Mrange')
    Mrange=1:M;
end
if exist('MrangeOverride')
    Mrange = MrangeOverride;
end

blank = nan(M,12);
Acc = blank;
AccPred = blank;

for i=Mrange
    i
    try
        clear temp SigT
        %temp=load([txt{ID(i),2} txt{ID(i),1}]); %Errors, edited by SS 8/4/21
        filename = [settings.patients{i,2} settings.patients{i,1} '.mat'];
        w=load(filename,'WakeSleepInfo','Info');
        clear WSinfo
        try 
            WSinfo = w.WakeSleepInfo.WSinfo;
        catch
            try
                WSinfo = w.Info.WSinfo;
            end
        end
        Acc(i,1:length(WSinfo.Acc))=WSinfo.Acc;
        AccPred(i,1:length(WSinfo.AccPred))=WSinfo.AccPred;
        
    catch
        i
        'failed'
    end
end


AccT = array2table(Acc);
AccPredT = array2table(AccPred);

WSInfoT = [AccT AccPredT];
WSInfoT.MaxAcc = max(Acc')';
WSInfoT.MaxAccPred = max(AccPred')';

