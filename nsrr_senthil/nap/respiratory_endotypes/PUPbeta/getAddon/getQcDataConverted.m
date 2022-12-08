function T = getQcDataConverted(MrangeOverride)
global settings
ConvertedFileDir = string(strcat(settings.ConvertedDirectory,settings.patients(:,1),'.mat'));
M = length(ConvertedFileDir);

T = table();
blankcol = nan(M,1);

T.HaveConverted = blankcol;
T.HaveFlowQc = blankcol;
T.PrUpright = blankcol;
T.SNRwindow = blankcol;
T.HaveWSinfo = blankcol;
T.WSMaxAcc = blankcol;
T.WSMaxAccPred = blankcol;

if exist('MrangeOverride')
    Mrange=MrangeOverride;
else
    Mrange=1:M;
end

for i=Mrange
    disp(i);
    try 
        T.HaveConverted(i) = 1*(exist(ConvertedFileDir(i))==2);
        if T.HaveConverted(i)==1
            W = load(ConvertedFileDir(i),'Info');

            T.HaveFlowQc(i) = isfield(W.Info,'FlowQ');
            if T.HaveFlowQc(i)
                try
                    T.PrUpright(i) = W.Info.FlowQ.PrUpright;
                    T.SNRwindow(i) = W.Info.FlowQ.SNRwindow;
                end
            end
            T.HaveWSinfo(i) = isfield(W.Info,'WSinfo');
            if T.HaveWSinfo(i)
                try
                    T.WSMaxAcc(i) = max(W.Info.WSinfo.Acc);
                    T.WSMaxAccPred(i) = max(W.Info.WSinfo.AccPred);
                end
            end
        end
    end
end




