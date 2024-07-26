function [PerfT,Raw,Cont] = PredictiveValue(criteriaR,PredT,Yvariable)

% disp('Note: Output Structure Variable Modified 7/7/2020')
%%
I = isnan(PredT)|isnan(criteriaR);

criteriaR=1*criteriaR; %actual

%criteriaR = criteriaR(randperm(length(criteriaR)));

TP = nansum(criteriaR(:).*PredT(:));
FP = nansum((1-criteriaR(:)).*PredT(:));
TN = nansum((1-criteriaR(:)).*(1-PredT(:)));
FN = nansum(criteriaR(:).*(1-PredT(:)));
N = TP+FP+TN+FN;

Raw = [TP FP TN FN]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUT
chancePPV = (TP+FN)/N; %likelihood of being positive
chanceNPV = (TN+FP)/N;
chanceSens = (TP+FP)/N; %likelihood of guessing positive
chanceSpec = (TN+FN)/N;
chanceAcc = chancePPV*chanceSens+chanceNPV*chanceSpec; 

Sens = TP/(TP+FN);
    sem = (Sens*(1-Sens)/(TP+FN)).^0.5;
    p = 2*[1-normcdf(abs(chanceSens-Sens)./sem,0,1)];
    performance.T.Sens = [Sens sem chanceSens p TP TP+FN]';
Spec = TN/(TN+FP);
    sem = (Spec*(1-Spec)/(TN+FP)).^0.5;
    p = 2*[1-normcdf(abs(chanceSpec-Spec)./sem,0,1)];
    performance.T.Spec = [Spec sem chanceSpec p TN TN+FP]';
Acc = (TP+TN)/N;
    sem = (Acc*(1-Acc)/(N)).^0.5;
    p = 2*[1-normcdf(abs(chanceAcc-Acc)./sem,0,1)];
    performance.T.Acc = [Acc sem chanceAcc p TP+TN N]';
PPV = TP/(TP+FP);
    sem = (PPV*(1-PPV)/(TP+FP)).^0.5;
    p = 2*[1-normcdf(abs(chancePPV-PPV)./sem,0,1)]; %recently added the 2* because p=0.05 at normcdf=0.975, i.e. 1-alpha/2 for 2-tailed test.
    performance.T.PPV = [PPV sem chancePPV p TP TP+FP]';
NPV = TN/(TN+FN);
    sem = (NPV*(1-NPV)/(TN+FN)).^0.5;
    p = 2*[1-normcdf(abs(chanceNPV-NPV)./sem,0,1)];
    performance.T.NPV = [NPV sem chanceNPV p TN TN+FN]';

PerfT = struct2table(performance.T);
PerfT.Properties.RowNames = {'Value','SEM','ChanceValue','P','Ncorrect','Ntotal'};

predY_mean=nanmean(Yvariable(PredT==1));
predY_sem=nanstd(Yvariable(PredT==1))/nansum(PredT==1)^0.5;
predN_mean=nanmean(Yvariable(PredT==0));
predN_sem=nanstd(Yvariable(PredT==0))/nansum(PredT==0)^0.5;
[~,phighvslow_ttest]=ttest2(Yvariable(PredT==1),Yvariable(PredT==0));
predY_median=nanmedian(Yvariable(PredT==1));
predN_median=nanmedian(Yvariable(PredT==0));
try
    phighvslow_ranksum=ranksum(Yvariable(PredT==1),Yvariable(PredT==0));
catch me
    phighvslow_ranksum=NaN;      
end

PredY = [predY_mean predY_sem predY_median TP TP+FP]';
PredN = [predN_mean predN_sem predN_median TN TN+FN]';
P = [phighvslow_ttest NaN phighvslow_ranksum NaN NaN]';

Cont = table(PredY,PredN,P);
Cont.Properties.RowNames = {'Mean','SEM','Median','Ncorrect','Ntotal'};




