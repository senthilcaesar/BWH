function performance = PredictiveValue(criteriaR,PredT,Yvariable)
%%
criteriaR=1*criteriaR; %actual

%criteriaR = criteriaR(randperm(length(criteriaR)));

TP = nansum(criteriaR(:).*PredT(:));
FP = nansum((1-criteriaR(:)).*PredT(:));
TN = nansum((1-criteriaR(:)).*(1-PredT(:)));
FN = nansum(criteriaR(:).*(1-PredT(:)));
N = TP+FP+TN+FN;

performance.TP_FP_TN_FN = [TP FP TN FN];
chancePPV = (TP+FN)/N; %likelihood of being positive
chanceNPV = (TN+FP)/N;
chanceSens = (TP+FP)/N; %likelihood of guessing positive
chanceSpec = (TN+FN)/N;
chanceAcc = chancePPV*chanceSens+chanceNPV*chanceSpec; 

PPV = TP/(TP+FP);
    sem = (PPV*(1-PPV)/(TP+FP)).^0.5;
    p = 2*[1-normcdf(abs(chancePPV-PPV)./sem,0,1)]; %recently added the 2* because p=0.05 at normcdf=0.975, i.e. 1-alpha/2 for 2-tailed test.
    performance.PPV_sem_chance_p = [PPV sem chancePPV p];
NPV = TN/(TN+FN);
    sem = (NPV*(1-NPV)/(TN+FN)).^0.5;
    p = 2*[1-normcdf(abs(chanceNPV-NPV)./sem,0,1)];
    performance.NPV_sem_chance_p = [NPV sem chanceNPV p];
Sens = TP/(TP+FN);
    sem = (Sens*(1-Sens)/(TP+FN)).^0.5;
    p = 2*[1-normcdf(abs(chanceSens-Sens)./sem,0,1)];
    performance.Sens_sem_chance_p = [Sens sem chanceSens p];
Spec = TN/(TN+FP);
    sem = (Spec*(1-Spec)/(TN+FP)).^0.5;
    p = 2*[1-normcdf(abs(chanceSpec-Spec)./sem,0,1)];
    performance.Spec_sem_chance_p = [Spec sem chanceSpec p];
Acc = (TP+TN)/N;
    sem = (Acc*(1-Acc)/(N)).^0.5;
    p = 2*[1-normcdf(abs(chanceAcc-Acc)./sem,0,1)];
    performance.Acc_sem_chance_p = [Acc sem chanceAcc p];

performance.predY_mean=nanmean(Yvariable(PredT==1));
performance.predY_sem=nanstd(Yvariable(PredT==1))/nansum(PredT==1)^0.5;
performance.predN_mean=nanmean(Yvariable(PredT==0));
performance.predN_sem=nanstd(Yvariable(PredT==0))/nansum(PredT==0)^0.5;
[~,performance.phighvslow_ttest]=ttest2(Yvariable(PredT==1),Yvariable(PredT==0));
performance.predY_median=nanmedian(Yvariable(PredT==1));
performance.predN_median=nanmedian(Yvariable(PredT==0));
    
try
    performance.phighvslow_ranksum=ranksum(Yvariable(PredT==1),Yvariable(PredT==0));
catch me
    disp('me.message');
    performance.phighvslow_ranksum=NaN;
end


