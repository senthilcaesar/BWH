function [thresX,AUC,SEM,p,posclass,SensSpec]=ROCAUCSEM(labels,y,ShowFigures)

%labels = x<thres;

%test to find posclass
[~,~,~,AUC0,~] = perfcurve(labels,y,0);
[~,~,~,AUC1,~] = perfcurve(labels,y,1);
posclass = single(AUC1>AUC0);   

%rerun AUC
[X,Y,T,AUC,OPTTHRES] = perfcurve(labels,y,AUC1>AUC0);
if ShowFigures
    %figure(500); clf(figure(500)); 
    plot(X,Y); hold on; box('off');
    axis square; r(1) = refline(1,0); r(1).Color=[0.5 0.5 0.5];
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title('ROC');
end

%find best Xthreshold and Sens/Spec
thresX = T((X==OPTTHRES(1))&(Y==OPTTHRES(2)));
SensSpec = [OPTTHRES(2) 1-OPTTHRES(1)]; 

%Calculate SEM and p
N1 = sum(labels==1); %must be the positive, i.e. abnormal, group being detected
N2 = sum(labels==0);
Q1 = AUC/(2-AUC);
Q2 = 2*AUC^2/(1+AUC);
SEM = ((AUC*(1-AUC)+(N1-1)*(Q1-AUC^2)+(N2-1)*(Q2-AUC^2))/(N1*N2))^0.5;
%from Hanley McNeill 1982 
%accessed free at http://www.med.mcgill.ca/epidemiology/hanley/software/Hanley_McNeil_Radiology_82.pdf

p = 2*[1-normcdf(abs(AUC-0.5)./SEM,0,1)];