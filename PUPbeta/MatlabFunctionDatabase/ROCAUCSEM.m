function [thresX,AUC,SEM,p,posclass,SensSpec]=ROCAUCSEM(labels,y,orig)

%labels = x<thres;

%test to find posclass
[~,~,~,AUC0,~] = perfcurve(labels,y,0);
[~,~,~,AUC1,~] = perfcurve(labels,y,1);
posclass = AUC1>AUC0;   

%rerun AUC
clear data
[data.X,data.Y,data.T,AUC,OPTTHRES] = perfcurve(labels,y,posclass);
data = struct2table(data);

data.Spec = 1-data.X;
data.SensPlusSpec = data.Y + (1-data.X);
data.SensPlusSpec = data.Y + 0.90*(1-data.X);

[~,I]=max(data.SensPlusSpec);

    if exist('orig') & orig==1
        %find best Xthreshold and Sens/Spec
        thresX = data.T((data.X==OPTTHRES(1))&(data.Y==OPTTHRES(2)));
        SensSpec = [OPTTHRES(2) 1-OPTTHRES(1)]; %same
    else
        thresX=mean(data.T(I:(I+1)));
        SensSpec = [data.Y(I) 1-data.X(I)];
    end

%Calculate SEM and p
N1 = sum(labels==1); %must be the positive, i.e. abnormal, group being detected
N2 = sum(labels==0);
Q1 = AUC/(2-AUC);
Q2 = 2*AUC^2/(1+AUC);
SEM = ((AUC*(1-AUC)+(N1-1)*(Q1-AUC^2)+(N2-1)*(Q2-AUC^2))/(N1*N2))^0.5;
%from Hanley McNeill 1982 
%accessed free at http://www.med.mcgill.ca/epidemiology/hanley/software/Hanley_McNeil_Radiology_82.pdf

p = 2*[1-normcdf(abs(AUC-0.5)./SEM,0,1)];