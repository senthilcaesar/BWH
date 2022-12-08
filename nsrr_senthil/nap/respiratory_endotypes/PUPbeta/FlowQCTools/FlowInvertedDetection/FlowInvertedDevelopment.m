clear all; close all; clc;

if 0
%% Training set
files=dir('D:\MrOS\Visit1\Analyzed\*.mat');

for jj=1:200
    jj
    clear BreathDataTable BreathFLDataTable
    load([files(jj).folder '\' files(jj).name])

clear diffcurr diffprev COVVTi
for i = 1:length(BreathDataTable{1})
    if istable(BreathDataTable{1}{i})
        Ttot = BreathDataTable{1}{i}.Time_end - BreathDataTable{1}{i}.Time_start;
        VI = BreathDataTable{1}{i}.VI;
        VT = VI./Ttot;
        %VT = ((VTi*Te)+(VTe*Ti))/(Ti + Te);
        Ti = BreathDataTable{1}{i}.Time_mid - BreathDataTable{1}{i}.Time_start;
        Te = Ttot - Ti;
        
        VTi_VT = BreathFLDataTable{1}{i}.VTi_VT_O;
        VTe_VT = BreathFLDataTable{1}{i}.VTe_VT_O;
        
        VTi = VTi_VT.*VT;
        VTe = VTe_VT.*VT;
        VTi(VT==0)=0;
        VTe(VT==0)=0;
        
        VTtest = ((VTi.*Te)+(VTe.*Ti))./(Ti + Te);
        
%         temp = [VTi VTe]
        
        Apnea_B = (VT==0)*1;
        
        %% Check if inverted trace
        [diffcurr(i,1),diffprev(i,1),COVVTi(i,1)] = FlowInvertedParameters(VTi,VTe,VT,Apnea_B);

        if diffcurr(i)>diffprev(i)
            if settings.verbose
                disp(['Warning: Flow signal appears upside down, F=' num2str(-100*(diffprev(i)-diffcurr(i)),2)]);
            end
        end
        
    else
        diffcurr(i,1)=NaN;
        diffprev(i,1)=NaN;
        COVVTi(i,1)=NaN;
    end
end


appearsupsidedown=1*(diffcurr>diffprev);
appearsupsidedown(isnan(COVVTi))=NaN;
upright = 1 + 0*appearsupsidedown;

diffcurroverprev = diffcurr./diffprev;
SubId=repmat(string(files(jj).name),size(diffcurroverprev,1),1);
SubNum=jj*ones(size(diffcurroverprev,1),1);
T1temp = table(SubNum,SubId,diffcurr,diffprev,diffcurroverprev,COVVTi,appearsupsidedown,upright);
if exist('T1')
    T1=[T1;T1temp];
else
    T1=T1temp;
end

end

T2 = T1;
T2.diffprev = T1.diffcurr;
T2.diffcurr = T1.diffprev;
T2.diffcurroverprev = 1./T1.diffcurroverprev;
T2.upright = 1-T1.upright;

T3 = [T1;T2];

%build model with large table all patients.

mdlUpright = fitglm(T3,'upright ~ diffcurroverprev ','Distribution','Binomial')
%test model on subjects separately;
T1.PappearsInverted = 1-predict(mdlUpright,T1);
T2.PappearsInverted = 1-predict(mdlUpright,T2);

T3.AbsError = abs(T3.upright - predict(mdlUpright,T3));
mdlError = fitglm(T3,'AbsError ~ COVVTi','Distribution','Binomial')

T1.Perror = predict(mdlError,T1);
T2.Perror = predict(mdlError,T2);
%figure(); histogram(T1.Perror);

[temp,indx] = unique(T1.SubNum);
clear Tall;
for i=1:length(temp)
    I = T1.SubNum==temp(i);
    Tall.MeanP(2*(i-1)+1,1)=nanmean(T1.PappearsInverted(I));
    Tall.MeanP(2*(i-1)+2,1)=nanmean(T2.PappearsInverted(I));
    Tall.MedianP(2*(i-1)+1,1)=nanmedian(T1.PappearsInverted(I));
    Tall.MedianP(2*(i-1)+2,1)=nanmedian(T2.PappearsInverted(I));
    Tall.MeanPw(2*(i-1)+1,1)=nanmean(T1.PappearsInverted(I).*(1-T1.Perror(I)))/nanmean(1-T1.Perror(I));
    Tall.MeanPw(2*(i-1)+2,1)=nanmean(T2.PappearsInverted(I).*(1-T2.Perror(I)))/nanmean(1-T2.Perror(I));
%     err1 = T1.Perror(I); err1 = err1*2; err1(err1>1)=1; 
%     err2 = T2.Perror(I); err2 = err2*2; err2(err2>1)=1;
%     Tall.MeanPw2(2*(i-1)+1,1)=nanmean(T1.PappearsInverted(I).*(1-err1))/nanmean(1-err1);
%     Tall.MeanPw2(2*(i-1)+2,1)=nanmean(T2.PappearsInverted(I).*(1-err2))/nanmean(1-err2);
    Tall.Upright(2*(i-1)+1,1)=1;
    Tall.Upright(2*(i-1)+2,1)=0;
    Tall.Subj(2*(i-1)+1,1)=T1.SubId(indx(i));
    Tall.Subj(2*(i-1)+2,1)=T1.SubId(indx(i));
end
Tall = struct2table(Tall);
%%

mdlSubjUpright = fitglm(Tall,'Upright ~ MeanPw','Distribution','Binomial')
Tall.PappearsUpright = predict(mdlSubjUpright,Tall);
Tall.Error = abs(Tall.Upright-Tall.PappearsUpright)>=0.5;
Pred = (Tall.PappearsUpright>0.5)*1;
    Pred(isnan(Tall.PappearsUpright))=NaN;
ErrIDs = unique(Tall.Subj(find(Tall.Error==1)));
[PerfT,raw,~] = PredictiveValue(Tall.Upright,Pred,Tall.Upright)
%% Note these were visually inspected and are truly upright already:
%     "MrOS_Apr2020_1056.mat"
%     "MrOS_Apr2020_1120.mat"


%% Debug
if 0
%load converted file here
figure(3)
plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,2))
end

%% Keep these things
FlowInvertedDetector.mdlUpright = mdlUpright;
FlowInvertedDetector.mdlError = mdlError;
FlowInvertedDetector.mdlSubjUpright = mdlSubjUpright;
save FlowInvertedDetectorFull FlowInvertedDetector

FlowInvertedDetector.mdlUpright = compact(mdlUpright);
FlowInvertedDetector.mdlError = compact(mdlError);
FlowInvertedDetector.mdlSubjUpright = compact(mdlSubjUpright);
save FlowInvertedDetector FlowInvertedDetector
end

%% Apply to one subject, output is PrUpright (0 means inverted, flow insp is down)

