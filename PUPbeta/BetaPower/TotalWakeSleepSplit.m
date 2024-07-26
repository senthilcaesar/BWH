
mdl4Alt = compact(fitglm(WStbl,['WakeNoAR ~ Pbeta_ref + Palpha_ref + Ptheta_ref + Pdelta_ref + Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeAR),'weights',weights))
mdl4Alt.Rsquared.Ordinary

mdlA = mdl4Alt;
%mdlA = mdlARrestart;


%%

nanmean(WStbl.Ptheta_ref(logical(WStbl.ExcludeAR)==0 & WStbl.WakeNoAR==1)) %constant = 1.2263
nanmean(WStbl.Pdelta_ref(logical(WStbl.ExcludeAR)==0 & WStbl.WakeNoAR==1)) %constant = 1.4154

T = WStbl(:,{'Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'});
T{:,{'Ptheta_ref'}} = T{:,{'Ptheta_ref'}}*0 + 1.2263;
T{:,{'Pdelta_ref'}} = T{:,{'Pdelta_ref'}}*0 + 1.4154;

WStbl.WakeIntensity = logit(predict(mdlA,T));
WStbl.WSBalance = logit(predict(mdlA,WStbl));
WStbl.SleepIntensity =  WStbl.WakeIntensity - WStbl.WSBalance;
% 
% 
% T = WStbl(:,{'Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'});
% T{:,{'Ptheta_ref'}} = T{:,{'Ptheta_ref'}}*0 + nanmean(WStbl.Ptheta_ref(logical(WStbl.ExcludeAR)==0 & WStbl.WakeNoAR==1));
% T{:,{'Pdelta_ref'}} = T{:,{'Pdelta_ref'}}*0 + nanmean(WStbl.Pdelta_ref(logical(WStbl.ExcludeAR)==0 & WStbl.WakeNoAR==1));
% T{:,{'Pbeta_ref'}} = T{:,{'Pbeta_ref'}}*0 ;
% T{:,{'Palpha_ref'}} = T{:,{'Palpha_ref'}}*0 ;
% temp = logit(predict(mdlA,T));

plusN=100;
figure(11+plusN); clf(11+plusN); set(gcf,'color',[1 1 1]);
XvalH = WStbl.WakeIntensity;
dStep=0.1;
Centers=-10:dStep:10;
run SleepHistograms

figure(12+plusN); clf(12+plusN); set(gcf,'color',[1 1 1]);
XvalH = WStbl.SleepIntensity;
dStep=0.1;
Centers=-10:dStep:10;
run SleepHistograms

figure(13+plusN); clf(13+plusN); set(gcf,'color',[1 1 1]);
XvalH = WStbl.WSBalance;
dStep=0.1;
Centers=-10:dStep:10;
run SleepHistograms


%%
temp = unique(WStbl.Subjn);

%Exclude_ = WStbl.NoiseBinary==1 | WStbl.SpO2off==1 | WStbl.Pbeta==-Inf | isnan(WStbl.Epochs) | WStbl.Epochs<0 | WStbl.Epochs>4;
%WStbl.ExcludeAR = Exclude | WStbl.Epochs==4 & WStbl.EventsAr==0 | WStbl.Epochs~=4 & WStbl.EventsAr==1;
%WStbl.ExcludeAR = logical(WStbl.ExcludeAR);

WStbl.ExcludeW = WStbl.Epochs==4 | WStbl.Exclude;% & string(WStbl.Subj)~=string(UniqueSubjList(1:10:100));

adjustedbalance=1.685; %1.5 before adding MESA; higher number increases FP and lowers FN
balanceAR = nanmean(WStbl.EventsAr(~WStbl.ExcludeW))*adjustedbalance %adjusted to balance FP and FN using *adjustedbalance

weightsAR = 0*WStbl.EventsAr;
weightsAR(WStbl.EventsAr>0.5)=1-balanceAR;
weightsAR(WStbl.EventsAr<=0.5)=balanceAR;
b(1)=sum(weightsAR(WStbl.EventsAr>0.5 & ~WStbl.ExcludeW))
b(2)=sum(weightsAR(WStbl.EventsAr<0.5 & ~WStbl.ExcludeW))
b(1)/(sum(b))

WStbl.L1 = [NaN;WStbl.WSBalance(1:end-1)];
WStbl.L3 = [NaN;NaN;NaN;WStbl.WSBalance(1:end-3)];
WStbl.L5 = [NaN;NaN;NaN;NaN;NaN;WStbl.WSBalance(1:end-5)];
WStbl.L7 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.WSBalance(1:end-7)];
WStbl.L9 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.WSBalance(1:end-9)];
WStbl.N1 = [WStbl.WSBalance(2:end);NaN];
WStbl.N3 = [WStbl.WSBalance(4:end);NaN;NaN;NaN];
WStbl.N5 = [WStbl.WSBalance(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7 = [WStbl.WSBalance(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9 = [WStbl.WSBalance(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];

WStbl.L1w = [NaN;WStbl.WakeIntensity(1:end-1)];
WStbl.L3w = [NaN;NaN;NaN;WStbl.WakeIntensity(1:end-3)];
WStbl.L5w = [NaN;NaN;NaN;NaN;NaN;WStbl.WakeIntensity(1:end-5)];
WStbl.L7w = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.WakeIntensity(1:end-7)];
WStbl.L9w = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.WakeIntensity(1:end-9)];
WStbl.N1w = [WStbl.WakeIntensity(2:end);NaN];
WStbl.N3w = [WStbl.WakeIntensity(4:end);NaN;NaN;NaN];
WStbl.N5w = [WStbl.WakeIntensity(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7w = [WStbl.WakeIntensity(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9w = [WStbl.WakeIntensity(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];


WStbl.L1s = [NaN;WStbl.SleepIntensity(1:end-1)];
WStbl.L3s = [NaN;NaN;NaN;WStbl.SleepIntensity(1:end-3)];
WStbl.L5s = [NaN;NaN;NaN;NaN;NaN;WStbl.SleepIntensity(1:end-5)];
WStbl.L7s = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.SleepIntensity(1:end-7)];
WStbl.L9s = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.SleepIntensity(1:end-9)];
WStbl.N1s = [WStbl.SleepIntensity(2:end);NaN];
WStbl.N3s = [WStbl.SleepIntensity(4:end);NaN;NaN;NaN];
WStbl.N5s = [WStbl.SleepIntensity(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7s = [WStbl.SleepIntensity(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9s = [WStbl.SleepIntensity(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];



sum(WStbl.ExcludeW)
UniqueSubjListi = [1;(find(diff(WStbl.Subjn)>0)+1);height(WStbl)];
for i=1:length(UniqueSubjListi)-1
    WStbl.ExcludeW(UniqueSubjListi(i):UniqueSubjListi(i)+15)=1;
    WStbl.ExcludeW(UniqueSubjListi(i+1):-1:UniqueSubjListi(i+1)-15)=1;
end
sum(WStbl.ExcludeW)

%%
T = WStblTest(:,{'Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'});
T{:,{'Ptheta_ref'}} = T{:,{'Ptheta_ref'}}*0 + 1.2263;
T{:,{'Pdelta_ref'}} = T{:,{'Pdelta_ref'}}*0 + 1.4154;

WStblTest.WakeIntensity = logit(predict(mdlA,T));
WStblTest.WSBalance = logit(predict(mdlA,WStblTest));
WStblTest.SleepIntensity =  WStblTest.WakeIntensity - WStblTest.WSBalance;

WStblTest.L1w = [NaN;WStblTest.WakeIntensity(1:end-1)];
WStblTest.L3w = [NaN;NaN;NaN;WStblTest.WakeIntensity(1:end-3)];
WStblTest.L5w = [NaN;NaN;NaN;NaN;NaN;WStblTest.WakeIntensity(1:end-5)];
WStblTest.L7w = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.WakeIntensity(1:end-7)];
WStblTest.L9w = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.WakeIntensity(1:end-9)];
WStblTest.N1w = [WStblTest.WakeIntensity(2:end);NaN];
WStblTest.N3w = [WStblTest.WakeIntensity(4:end);NaN;NaN;NaN];
WStblTest.N5w = [WStblTest.WakeIntensity(6:end);NaN;NaN;NaN;NaN;NaN];
WStblTest.N7w = [WStblTest.WakeIntensity(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N9w = [WStblTest.WakeIntensity(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];


WStblTest.L1s = [NaN;WStblTest.SleepIntensity(1:end-1)];
WStblTest.L3s = [NaN;NaN;NaN;WStblTest.SleepIntensity(1:end-3)];
WStblTest.L5s = [NaN;NaN;NaN;NaN;NaN;WStblTest.SleepIntensity(1:end-5)];
WStblTest.L7s = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.SleepIntensity(1:end-7)];
WStblTest.L9s = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.SleepIntensity(1:end-9)];
WStblTest.N1s = [WStblTest.SleepIntensity(2:end);NaN];
WStblTest.N3s = [WStblTest.SleepIntensity(4:end);NaN;NaN;NaN];
WStblTest.N5s = [WStblTest.SleepIntensity(6:end);NaN;NaN;NaN;NaN;NaN];
WStblTest.N7s = [WStblTest.SleepIntensity(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N9s = [WStblTest.SleepIntensity(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];



sum(WStblTest.ExcludeW)
UniqueSubjListi = [1;(find(diff(WStblTest.Subjn)>0)+1);height(WStblTest)];
for i=1:length(UniqueSubjListi)-1
    WStblTest.ExcludeW(UniqueSubjListi(i):UniqueSubjListi(i)+15)=1;
    WStblTest.ExcludeW(UniqueSubjListi(i+1):-1:UniqueSubjListi(i+1)-15)=1;
end
sum(WStblTest.ExcludeW)

%% Playing

WStbl.WIminusL1w = WStbl.WakeIntensity - WStbl.L1w;
WStbl.WIminusL3w = WStbl.WakeIntensity - WStbl.L3w;
WStbl.WIminusL5w = WStbl.WakeIntensity - WStbl.L5w;
WStbl.WIminusL7w = WStbl.WakeIntensity - WStbl.L7w;
WStbl.WIminusL9w = WStbl.WakeIntensity - WStbl.L9w;
WStbl.WIminusN1w = WStbl.WakeIntensity - WStbl.N1w;
WStbl.WIminusN3w = WStbl.WakeIntensity - WStbl.N3w;
WStbl.WIminusN5w = WStbl.WakeIntensity - WStbl.N5w;
WStbl.WIminusN7w = WStbl.WakeIntensity - WStbl.N7w;
WStbl.WIminusN9w = WStbl.WakeIntensity - WStbl.N9w;
WStbl.SIminusL1s = WStbl.SleepIntensity - WStbl.L1s;
WStbl.SIminusL3s = WStbl.SleepIntensity - WStbl.L3s;
WStbl.SIminusL5s = WStbl.SleepIntensity - WStbl.L5s;
WStbl.SIminusL7s = WStbl.SleepIntensity - WStbl.L7s;
WStbl.SIminusL9s = WStbl.SleepIntensity - WStbl.L9s;
WStbl.SIminusN1s = WStbl.SleepIntensity - WStbl.N1s;
WStbl.SIminusN3s = WStbl.SleepIntensity - WStbl.N3s;
WStbl.SIminusN5s = WStbl.SleepIntensity - WStbl.N5s;
WStbl.SIminusN7s = WStbl.SleepIntensity - WStbl.N7s;
WStbl.SIminusN9s = WStbl.SleepIntensity - WStbl.N9s;


WStbl.L1b = [NaN;WStbl.Pbeta_ref(1:end-1)];
WStbl.L3b = [NaN;NaN;NaN;WStbl.Pbeta_ref(1:end-3)];
WStbl.L5b = [NaN;NaN;NaN;NaN;NaN;WStbl.Pbeta_ref(1:end-5)];
WStbl.L7b = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Pbeta_ref(1:end-7)];
WStbl.L9b = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Pbeta_ref(1:end-9)];
WStbl.N1b = [WStbl.Pbeta_ref(2:end);NaN];
WStbl.N3b = [WStbl.Pbeta_ref(4:end);NaN;NaN;NaN];
WStbl.N5b = [WStbl.Pbeta_ref(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7b = [WStbl.Pbeta_ref(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9b = [WStbl.Pbeta_ref(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.L1a = [NaN;WStbl.Palpha_ref(1:end-1)];
WStbl.L3a = [NaN;NaN;NaN;WStbl.Palpha_ref(1:end-3)];
WStbl.L5a = [NaN;NaN;NaN;NaN;NaN;WStbl.Palpha_ref(1:end-5)];
WStbl.L7a = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Palpha_ref(1:end-7)];
WStbl.L9a = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Palpha_ref(1:end-9)];
WStbl.N1a = [WStbl.Palpha_ref(2:end);NaN];
WStbl.N3a = [WStbl.Palpha_ref(4:end);NaN;NaN;NaN];
WStbl.N5a = [WStbl.Palpha_ref(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7a = [WStbl.Palpha_ref(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9a = [WStbl.Palpha_ref(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.L1t = [NaN;WStbl.Ptheta_ref(1:end-1)];
WStbl.L3t = [NaN;NaN;NaN;WStbl.Ptheta_ref(1:end-3)];
WStbl.L5t = [NaN;NaN;NaN;NaN;NaN;WStbl.Ptheta_ref(1:end-5)];
WStbl.L7t = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Ptheta_ref(1:end-7)];
WStbl.L9t = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Ptheta_ref(1:end-9)];
WStbl.N1t = [WStbl.Ptheta_ref(2:end);NaN];
WStbl.N3t = [WStbl.Ptheta_ref(4:end);NaN;NaN;NaN];
WStbl.N5t = [WStbl.Ptheta_ref(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7t = [WStbl.Ptheta_ref(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9t = [WStbl.Ptheta_ref(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.L1d = [NaN;WStbl.Pdelta_ref(1:end-1)];
WStbl.L3d = [NaN;NaN;NaN;WStbl.Pdelta_ref(1:end-3)];
WStbl.L5d = [NaN;NaN;NaN;NaN;NaN;WStbl.Pdelta_ref(1:end-5)];
WStbl.L7d = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Pdelta_ref(1:end-7)];
WStbl.L9d = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.Pdelta_ref(1:end-9)];
WStbl.N1d = [WStbl.Pdelta_ref(2:end);NaN];
WStbl.N3d = [WStbl.Pdelta_ref(4:end);NaN;NaN;NaN];
WStbl.N5d = [WStbl.Pdelta_ref(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N7d = [WStbl.Pdelta_ref(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9d = [WStbl.Pdelta_ref(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.L1td = WStbl.L1t.*WStbl.L1d;
WStbl.L3td = WStbl.L3t.*WStbl.L3d;
WStbl.L5td = WStbl.L5t.*WStbl.L5d;
WStbl.L7td = WStbl.L7t.*WStbl.L7d;
WStbl.L9td = WStbl.L9t.*WStbl.L9d;
WStbl.N1td = WStbl.N1t.*WStbl.N1d;
WStbl.N3td = WStbl.N3t.*WStbl.N3d;
WStbl.N5td = WStbl.N5t.*WStbl.N5d;
WStbl.N7td = WStbl.N7t.*WStbl.N7d;
WStbl.N9td = WStbl.N9t.*WStbl.N9d;

%%

mdlARorig = compact(fitglm(WStbl,...
    ['EventsAr ~ WSBalance + L9 + L7 + L5 + L3 + L1 + N1 + N3 + N5 + N7 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARorig.Rsquared.Ordinary

mdlARwsbalance = compact(fitglm(WStbl,...
    ['EventsAr ~ WSBalance'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARwsbalance.Rsquared.Ordinary


%% Playing
if 0
%mdlAR = compact(fitglm(WStbl,...
%    ['EventsAr ~ WSBalance + L9 + L8 + L7 + L6 + L5 + L4 + L3 + L2 + L1 + N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))

mdlARw = compact(fitglm(WStbl,...
    ['EventsAr ~ WakeIntensity + L9w + L7w + L5w + L3w + L1w + N1w + N3w + N5w + N7w + N9w'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARw.Rsquared.Ordinary

mdlARb = compact(fitglm(WStbl,...
    ['EventsAr ~ Pbeta_ref + L9b + L7b + L5b + L3b + L1b + N1b + N3b + N5b + N7b + N9b'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARb.Rsquared.Ordinary


mdlARba = compact(fitglm(WStbl,...
    ['EventsAr ~ Palpha_ref + L9a + L7a + L5a + L3a + L1a + N1a + N3a + N5a + N7a + N9a + Pbeta_ref + L9b + L7b + L5b + L3b + L1b + N1b + N3b + N5b + N7b + N9b'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARba.Rsquared.Ordinary

mdlARbat = compact(fitglm(WStbl,...
    ['EventsAr ~ Ptheta_ref + L9t + L7t + L5t + L3t + L1t + N1t + N3t + N5t + N7t + N9t + Palpha_ref + L9a + L7a + L5a + L3a + L1a + N1a + N3a + N5a + N7a + N9a + Pbeta_ref + L9b + L7b + L5b + L3b + L1b + N1b + N3b + N5b + N7b + N9b'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARbat.Rsquared.Ordinary

mdlARbatd = compact(fitglm(WStbl,...
    ['EventsAr ~ Pdelta_ref + L9d + L7d + L5d + L3d + L1d + N1d + N3d + N5d + N7d + N9d + Ptheta_ref + L9t + L7t + L5t + L3t + L1t + N1t + N3t + N5t + N7t + N9t + Palpha_ref + L9a + L7a + L5a + L3a + L1a + N1a + N3a + N5a + N7a + N9a + Pbeta_ref + L9b + L7b + L5b + L3b + L1b + N1b + N3b + N5b + N7b + N9b'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARbatd.Rsquared.Ordinary


mdlARrel = compact(fitglm(WStbl,...
    ['EventsAr ~ WIminusL1w + WIminusL3w + WIminusL5w + WIminusL7w + WIminusL9w + WIminusN1w + WIminusN3w + WIminusN5w + WIminusN7w + WIminusN9w + SIminusL1s + SIminusL3s + SIminusL5s + SIminusL7s + SIminusL9s + SIminusN1s + SIminusN3s + SIminusN5s + SIminusN7s + SIminusN9s'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARrel.Rsquared.Ordinary

mdlARbatd = compact(fitglm(WStbl,...
    ['EventsAr ~ Pdelta_ref*Ptheta_ref + L9td + L7td + L5td + L3td + L1td + N1td + N3td + N5td + N7td + N9td + Pdelta_ref + L9d + L7d + L5d + L3d + L1d + N1d + N3d + N5d + N7d + N9d + Ptheta_ref + L9t + L7t + L5t + L3t + L1t + N1t + N3t + N5t + N7t + N9t + Palpha_ref + L9a + L7a + L5a + L3a + L1a + N1a + N3a + N5a + N7a + N9a + Pbeta_ref + L9b + L7b + L5b + L3b + L1b + N1b + N3b + N5b + N7b + N9b'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARbatd.Rsquared.Ordinary


mdlARs = compact(fitglm(WStbl,...
    ['EventsAr ~ SleepIntensity + L9s + L7s + L5s + L3s + L1s + N1s + N3s + N5s + N7s + N9s'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARs.Rsquared.Ordinary

mdlAR = compact(fitglm(WStbl,...
    ['EventsAr ~ WSBalance + L9 + L7 + L5 + L3 + L1 + N1 + N3 + N5 + N7 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))

mdlARrestart = compact(fitglm(WStbl,['EventsAr ~ Pbeta_ref + Palpha_ref + Ptheta_ref + Pdelta_ref + Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeW),'weights',weightsAR))
mdlARrestart.Rsquared.Ordinary

[PrArousal,~] = predict(mdlARrestart,WStbl);
[x,y,t,AUC_WS,~] = perfcurve(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>0.5),1*(ARieF_pred(~WStbl.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_WS


WStbl.ARpred = predict(mdlAR,WStbl);
WStbl.ARpredw = predict(mdlARw,WStbl);
WStbl.ARpreds = predict(mdlARs,WStbl);

mdlARcombo = compact(fitglm(WStbl,...
    ['EventsAr ~ ARpredw*ARpreds'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARcombo.Rsquared.Ordinary

mdlARcombo = compact(fitglm(WStbl,...
    ['EventsAr ~ ARpredw*ARpreds + ARpredw*ARpred + ARpred*ARpreds'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARcombo.Rsquared.Ordinary

end
%% This one, ~1 min for fit
mdlAR = compact(fitglm(WStbl,...
    ['EventsAr ~ WakeIntensity + L9w + L7w + L5w + L3w + L1w + N1w + N3w + N5w + N7w + N9w + SleepIntensity + L9s + L7s + L5s + L3s + L1s + N1s + N3s + N5s + N7s + N9s'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlAR.Rsquared.Ordinary

%% Performance tests; on training data at present

mdlARtest = mdlARwsbalance;
%mdlARtest = mdlARorig;
%mdlARtest = mdlAR;
%+*mdlARtest = mdlARw;
mdlARtest.Rsquared.Ordinary
newARpred_ = predict(mdlARtest,WStbl);
WStbl.newARpredF = newARpred_;
WStbl.newARpred_logit = logit(newARpred_);
WStbl.newARpred = 1*(newARpred_>0.5);
WStbl.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceAR = PredictiveValue(1*(WStbl.EventsAr(~WStbl.ExcludeW)>0.5),1*(newARpred_(~WStbl.ExcludeW)>thres),WStbl.EventsAr(~WStbl.ExcludeW))
%PredictiveValue(criteriaR,PredT,Yvariable)
[thresX,AUC,SEM,p,posclass,SensSpec]=ROCAUCSEM(1*(WStbl.EventsAr(~WStbl.ExcludeW)>0.5),1*(newARpred_(~WStbl.ExcludeW)>thres));

%% Test
newARpred_ = predict(mdlARtest,WStblTest);
WStblTest.newARpredF = newARpred_;
WStblTest.newARpred_logit = logit(newARpred_);
WStblTest.newARpred = 1*(newARpred_>0.5);
WStblTest.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceARtest = PredictiveValue(1*(WStblTest.EventsAr(~WStblTest.ExcludeW)>0.5),1*(newARpred_(~WStblTest.ExcludeW)>thres),WStblTest.EventsAr(~WStblTest.ExcludeW))
%PredictiveValue(criteriaR,PredT,Yvariable)
[thresX,AUC,SEM,p,posclass,SensSpec]=ROCAUCSEM(1*(WStblTest.EventsAr(~WStblTest.ExcludeW)>0.5),1*(newARpred_(~WStblTest.ExcludeW)>thres));

%% Plot AR constants themselves
figure(210); clf(210);
set(gcf,'color',[1 1 1]);

mdlARtable = mdlAR.Coefficients;
%mdlARtable = mdl9b.Coefficients;
mdlARtable = mdlARtable([1 8:-1:4 2 9:13 18:-1:14 3 19:23],:)
categorical(mdlARtable.Properties.RowNames)

x = categorical(mdlARtable.Properties.RowNames);
x = reordercats(x,mdlARtable.Properties.RowNames);

h0=bar(x,mdlARtable.Estimate);

hold('on');

set(gca,'box','off','tickdir','out')

h=errorbar(x,mdlARtable.Estimate,mdlARtable.SE,'CapSize',18)
h.LineStyle='none';
h.Color=[0 0 0];
%h.YNegativeDelta(mdlARtable.Estimate>0)=0;
%h.YPositiveDelta(mdlARtable.Estimate<0)=0;
h.CapSize=12;

hold('on');

%redraw the bars 
h0=bar(x,mdlARtable.Estimate);
h0.FaceColor = [ 0.8510    0.3294    0.1020];
h0.EdgeColor=h0.FaceColor;

%%
if 0
   save mdlARwsbalance mdlAR performanceAR
end