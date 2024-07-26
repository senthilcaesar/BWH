function HRmeansleepT = getOvernightHR(MrangeOverride)
global AMasterSpreadsheet settings

%example: HRmeansleepT = getOvernightHR([1,2,3,6,7:15])

t_startGetData = clock;

%%

%[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26'); % only used incase of analyzed file

% [~,txt] = xlsread(AMasterSpreadsheet,1,'AD4:AE5000');
Study = settings.patients(:,1);
M = size(Study,1);

if ~isfield(settings,'Mrange')
    Mrange=1:M;
end
if exist('MrangeOverride')
    Mrange = MrangeOverride;
end

blank = nan(M,1);
FHRavailable = blank;
HRmeansleep = blank;
PRmeansleep = blank;
FPRavailable = blank;

for i=Mrange
    try
        clear temp SigT
        %temp=load([txt{ID(i),2} txt{ID(i),1}]); %Errors, edited by SS 8/4/21
        temp=load([Study{i}]);
        
        if isfield(temp, 'DataEventHypnog_Mat') % convert DataEventHypnog_Mat to SigT
            SigT=array2table(temp.DataEventHypnog_Mat);
            SigT.Properties.VariableNames = temp.ChannelsList;
        else
            SigT=temp.SigT;
        end
        clear temp 
        Epochs =SigT.Epochs;
        criteria = Epochs<4 & Epochs>=0; %sleep epochs only
        
        if any(SigT.Properties.VariableNames=="HR")
        HR = SigT.HR;
        
        HRmeansleep(i,1) = nanmean(HR(criteria));
        FHRavailable(i,1) = mean(~isnan(HR(criteria)));
        else
            HRmeansleep(i,1) = NaN;
            FHRavailable(i,1) = NaN;
        end
        
        if any(SigT.Properties.VariableNames=="Pulse")
        Pulse = SigT.Pulse;
        PRmeansleep(i,1) = nanmean(Pulse(criteria));
        FPRavailable(i,1) = mean(~isnan(Pulse(criteria)));
        else
            PRmeansleep(i,1) = NaN;
            FPRavailable(i,1) = NaN;
        end
    catch
        i
        'failed'
        HRmeansleep(i,1) = NaN;
            FHRavailable(i,1) = NaN;
            PRmeansleep(i,1) = NaN;
            FPRavailable(i,1) = NaN;
    end
end

%HRmeansleep = HRmeansleep(:);
HRmeansleepT = table(Study,HRmeansleep,FHRavailable,PRmeansleep,FPRavailable);

%QC
minFHR=0.05;
HRmeansleepT.HRmeansleep(HRmeansleepT.FHRavailable<minFHR)=NaN;
HRmeansleepT.PRmeansleep(HRmeansleepT.FPRavailable<minFHR)=NaN;
HRmeansleepT.HRmeansleep2 = nanmean([HRmeansleepT.HRmeansleep HRmeansleepT.PRmeansleep],2); %mean of both signals

