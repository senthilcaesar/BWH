addpath(genpath('G:\Dropbox (Personal)\PUPbeta_git\PUPbeta20190629\'));
logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p));
normalci = @(p,n) 1.96*((p.*(1-p)./n).^0.5);

%% Choose source table
if 0
    WStbl = BreathDataFullTableChosen;
else
    WStbl = BreathDataFullTablebest;    
end

%% Remove extra RICCADSA data
if 1
    I = find(MultiScorer>1)+40000;
    I2 = find(sum(WStbl.Subjn==I',2)>0); %Exp 4 is RICCADSA   
    WStbl(I2,:)=[];    
end

%% Exclude Phenotype
if 1
    I = find(floor(WStbl.Subjn/10000)==1); %Remove Phenotype
    WStbl(I,:)=[];
end
%% Exclude MESA
if 1
    I = find(floor(WStbl.Subjn/10000)==5); %Remove MESA
    WStbl(I,:)=[];
end


%% Choose exclusions, general, removal of AR, removal of W
WStbl.Exclude = WStbl.NoiseBinary==1 | WStbl.SpO2off==1 | WStbl.Pbeta==-Inf | isnan(WStbl.Epochs) | WStbl.Epochs<0 | WStbl.Epochs>4;
WStbl.ExcludeAR = WStbl.Exclude | WStbl.Epochs==4 & WStbl.EventsAr==0 | WStbl.Epochs~=4 & WStbl.EventsAr==1;
WStbl.ExcludeAR = logical(WStbl.ExcludeAR);
%WStblTest.Exclude = WStblTest.NoiseBinary==1 | WStblTest.SpO2off==1 | WStblTest.Pbeta==-Inf | isnan(WStblTest.Epochs) | WStblTest.Epochs<0 | WStblTest.Epochs>4;

%% Convert mV to uV in RICCADSA
if 1 
    I = find(floor(WStbl.Subjn/10000)==4); %Exp 4 is RICCADSA        
    if nanmedian(WStbl.Pbeta(I))<-2 %test of whether we have already corrected mV to uV
        WStbl.Pbeta(I)=WStbl.Pbeta(I)+6;
        WStbl.Palpha(I)=WStbl.Palpha(I)+6;
        WStbl.Ptheta(I)=WStbl.Ptheta(I)+6;
        WStbl.Pdelta(I)=WStbl.Pdelta(I)+6;
    end
end
%% Convert mV to uV in MESA
if 1 
    I = find(floor(WStbl.Subjn/10000)==5); %Exp 4 is RICCADSA        
    if nanmedian(WStbl.Pbeta(I))<-2 %test of whether we have already corrected mV to uV
        WStbl.Pbeta(I)=WStbl.Pbeta(I)+6;
        WStbl.Palpha(I)=WStbl.Palpha(I)+6;
        WStbl.Ptheta(I)=WStbl.Ptheta(I)+6;
        WStbl.Pdelta(I)=WStbl.Pdelta(I)+6;
    end
end

%%
NN=size(WStbl,1);
newcollist = {'Pbeta','Palpha','Ptheta','Pdelta'}
centilexlist = [5 25 50 75 95];

%% UniqueSubjList
UniqueSubjListi = [1;(find(diff(WStbl.Subjn)>0)+1)];
RefSubj = WStbl.Subjn(UniqueSubjListi);
sum(UniqueSubjListi)

%% START Generate Centiles data (few min) -- speed up?
clear temp
for i=1:length(newcollist)
    for j=1:length(centilexlist)
        newcol = newcollist{i};
        %newcol = 'Pdelta'
        centilex=centilexlist(j);
        newcol_ = [newcol num2str(centilex) 'p'];
        %disp([newcol_]);
        for n=1:length(RefSubj)
            disp([newcol_ ' ' num2str(n)]);
            rowrange = WStbl.Subjn==RefSubj(n);
            %tempdata = eval(['WStbl.' newcol '(rowrangeincl,:)']);           
            centiletemp = prctile(eval(['WStbl.' newcol '(rowrange & ~WStbl.Exclude,:)']),centilex);
            if any(strcmp(newcol_,WStbl.Properties.VariableNames))==0
                nrow = size(WStbl,1);
                eval(['WStbl.' newcol_ '=zeros(nrow,1)+NaN;']);
            end
            eval(['WStbl.' newcol_ '(rowrange)=centiletemp;']);
        end
    end
end

%% STEP 2: Make cols of prctile constants referenced to Pbeta5p
count = 1;
clear constantslist
constantsliststring=[];
ref = {'Pbeta5p'}%,'Palpha1p','Palpha5p','Palpha25p','Palpha75p','Palpha10p','Palpha75p','Palpha90p','Ptheta','Pdelta','Psigma'};
for i=1:length(newcollist)
    for j=1:length(centilexlist)
        temp=[newcollist{i} num2str(centilexlist(j)) 'p'];
        if strcmp(temp,ref)==1
            continue
        end
        constantslist{count}=temp;
        newcol = constantslist(count);
        count=count+1;
        eval(['WStbl.' ref{1} '_' newcol{1} '= WStbl.' ref{1} '- WStbl.' newcol{1} ';']);
        constantsliststring = [constantsliststring '+ ' ref{1} '_' newcol{1} ' '];
    end
end
clear constantslist2
for i=1:length(constantslist)
    constantslist2{1,i} = [ref{1} '_' constantslist{i}];
end


%% STEP 3: Reference all powers to Pbeta5p
%newcol = {'Pbeta','Palpha','Pdelta','Ptheta'};
refs = {'Pbeta5p'};
for j=1:length(refs)
    for i=1:length(newcollist)
        eval(['WStbl.' newcollist{i} '_' refs{j} '= WStbl.' newcollist{i} '- WStbl.' refs{j} ';']);
    end
end

%% Odds and Evens for devel + validation
if 1
    WStblBackup=WStbl;
    WStblTest=WStbl;
    removeevenforholdout=1;
    if removeevenforholdout
        I = find(rem(WStbl.Subjn,2)==0);
        WStbl(I,:)=[];
        I = find(rem(WStbl.Subjn,2)==1);
        WStblTest(I,:)=[];
    end
end

%% UniqueSubjList
UniqueSubjListi = [1;(find(diff(WStbl.Subjn)>0)+1)];
RefSubj = WStbl.Subjn(UniqueSubjListi);

%% STEP 4: run model with constants

WStbl.WakeNoAR = (WStbl.Epochs==4)*1;

% Weights
balancefactor=1%0.86
balance = nanmean(WStbl.WakeNoAR(~WStbl.ExcludeAR));
weights = zeros(length(WStbl.WakeNoAR),1);
weights(WStbl.WakeNoAR>0.5)=1-balance*balancefactor;
weights(WStbl.WakeNoAR<=0.5)=balance*balancefactor;
sum(weights(WStbl.WakeNoAR>0.5&~WStbl.ExcludeAR))
sum(weights(WStbl.WakeNoAR<=0.5&~WStbl.ExcludeAR))

meanbalancedY = nansum(weights(~WStbl.ExcludeAR).*WStbl.WakeNoAR(~WStbl.ExcludeAR))/nansum(weights(~WStbl.ExcludeAR))


%%
%Pbeta_Pbeta5p*Palpha_Pbeta5p +
equationmodel = ['WakeNoAR ~ Ptheta_Pbeta5p*Pdelta_Pbeta5p + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p + Palpha_Pbeta5p + Pdelta_Pbeta5p + Pdelta_Pbeta5p^2',constantsliststring];
mdl3b = compact(fitglm(WStbl,equationmodel,'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeAR),'weights',weights))
mdl3b.Rsquared.Ordinary

%% STEP 5: convert coeffs for constants to generate single reference value
%mdl3b = fitglm(WStbl,['EventsAr ~ Pbeta_Pbeta5p^2 + Pdelta_Pbeta5p^2 + Palpha_Pbeta5p^2 + Ptheta_Pbeta5p^2 + Psigma_Psigma5p^2 + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p*Pdelta_Pbeta5p' constantsliststring],'Distribution','binomial','Link','logit','weights',weights)
%mdl3b=compact(mdl3b);
%mdlX = mdl3bAlt;

mdlX = mdl3b;
constant = 1*(sum((string(mdlX.Coefficients.Properties.RowNames) == string(constantslist2)),2)==1);
constant(1)=NaN;
I1 = find(constant==1);
constantnames = mdlX.Coefficients.Properties.RowNames(I1);
coeffsconstant = mdlX.Coefficients.Estimate(I1);
I2 = find(constant==0);
coeffs1 = mdlX.Coefficients.Estimate(I2);
terms = 0;
temptemp=1/sum(coeffs1)*coeffsconstant;
for i=1:length(I1)
    terms = terms - 1/sum(coeffs1)*coeffsconstant(i)*eval(['WStbl.' constantnames{i}]);
end
WStbl.ref = WStbl.Pbeta5p + terms;

Coeffs = [1-sum(temptemp);temptemp];

clear Constants
for i=1:length(constantnames)
    Constants{i,1} = constantnames{i}(9:end);
end
Constants = [{'Pbeta5p'};Constants];
RefTable = table(Constants,Coeffs);
refmean = nanmean(WStbl.ref(UniqueSubjListi));
RefTable = [RefTable(1,:);RefTable];
RefTable(1,1) = {'Constant'};
RefTable{1,2} = -refmean;
RefTable

WStbl.ref = WStbl.ref - refmean;
%How to apply:
ref_=0;
for i=2:length(RefTable.Constants)
    ref_ = ref_ + RefTable.Coeffs(i)*eval(['WStbl.' RefTable.Constants{i}]);
end
ref_ = ref_ + RefTable.Coeffs(1);

listlist=[];
%newcol = {'Pbeta','Palpha','Pdelta','Ptheta','Psigma'};
refs = {'ref'};
for j=1:length(refs)
    for i=1:length(newcollist)
        listlist{length(listlist)+1} = [newcollist{i} '_' refs{j}];
        eval(['WStbl.' newcollist{i} '_' refs{j} '= WStbl.' newcollist{i} '- WStbl.' refs{j} ';']);
    end
end

%% Repeat for WStblTest
%How to apply:
ref_=0;
for i=2:length(RefTable.Constants)
    ref_ = ref_ + RefTable.Coeffs(i)*eval(['WStblTest.' RefTable.Constants{i}]);
end
ref_ = ref_ + RefTable.Coeffs(1);
WStblTest.ref = ref_;

listlist=[];
%newcol = {'Pbeta','Palpha','Pdelta','Ptheta','Psigma'};
refs = {'ref'};
for j=1:length(refs)
    for i=1:length(newcollist)
        listlist{length(listlist)+1} = [newcollist{i} '_' refs{j}];
        eval(['WStblTest.' newcollist{i} '_' refs{j} '= WStblTest.' newcollist{i} '- WStblTest.' refs{j} ';']);
    end
end


%% Plot variability in Ref

dStep=0.1;
Centers=-1:dStep:1.55;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
temp = (WStbl.ref(UniqueSubjListi));
temp(isnan(temp))=[];
[hall,edges] = histcounts(temp,Edges);
figure(20); clf(20); 
set(gcf,'color',[1 1 1]);
bar(Centers,hall,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
set(gca,'box','off','tickdir','out')

% plot constants themselves
figure(21);
set(gcf,'color',[1 1 1]);
x = categorical(["bananas" "apples" "cherries"]);
x = reordercats(x,{'bananas' 'apples' 'cherries'});
y = [14,12,7];
bar(x,y);

categorical(RefTable.Constants)

x = categorical(RefTable.Constants);
x = reordercats(x,RefTable.Constants);

bar(x,RefTable.Coeffs')
set(gca,'box','off','tickdir','out')

%% Unknown purpose:
% consstrings = RefTable.Constants(2:end);
% for i=1:length(consstrings)
% ind = find(strcmp(WStbl.Properties.VariableNames,consstrings(i)));
% tempval(i)=mean(WStbl{UniqueSubjListi,ind});
% end

figure(56); clf(56);
nrange = [4 42]
tempref=[WStbl.ref(UniqueSubjListi) [1:length(UniqueSubjListi)]']
Centers=-1.5:dStep:3;
dStep=0.01;

powerslist = {'beta','alpha','theta','delta'};
for n=nrange
    rowsi = [(UniqueSubjListi(n):UniqueSubjListi(n+1)-1)]';
    WStbl.ref(UniqueSubjListi(n));

for j=1:4
    
temp=eval(['WStbl.P' powerslist{j} '(rowsi);']);
temp(WStbl.Exclude(rowsi))=[];


Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
Edges(end)=Inf;
[h1,edges] = histcounts(temp,Edges);

area = sum([h1])*dStep;
h1=h1/area;


subplot(2,2,j);
    bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
    box('off');
    hold('on')
end
end

figure(57); clf(57);

powerslist = {'beta_ref','alpha_ref','theta_ref','delta_ref'};
for n=nrange
    rowsi = [(UniqueSubjListi(n):UniqueSubjListi(n+1)-1)]';

for j=1:4
    
temp=eval(['WStbl.P' powerslist{j} '(rowsi);']);
temp(WStbl.Exclude(rowsi))=[];

% dStep=0.1;
% Centers=-3:dStep:7;
% Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
% Edges(end)=Inf;
[h1,edges] = histcounts(temp,Edges);

area = sum([h1])*dStep;
h1=h1/area;

subplot(2,2,j);
    bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
    box('off');
    hold('on')
end
end
%% STEP 5b: rerun using single ref (for presentation, communication), excl arousals

% mdl3bXrefTop4_full_noweights_nobetasq = compact(fitglm(WStbl,['EventsAr ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit'))
% mdl3bXrefTop4_full_noweights_nobetasq.Rsquared.Ordinary

% WStbl.WakeNoAR = (WStbl.Epochs==4)*1;
% ExcludeAR = 0*WStbl.WakeNoAR;
%     ExcludeAR(WStbl.Epochs==4&WStbl.EventsAr==0)=1;
%     ExcludeAR(WStbl.Epochs~=4&WStbl.EventsAr==1)=1;
%     ExcludeAR = logical(ExcludeAR);

%no weights
% mdl3bXrefTop4_full_noweights_nobetasq_W_noweights = compact(fitglm(WStbl,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',ExcludeAR))
% mdl3bXrefTop4_full_noweights_nobetasq_W_noweights.Rsquared.Ordinary

%weights
mdl4Alt = compact(fitglm(WStbl,['WakeNoAR ~ Pbeta_ref + Palpha_ref + Ptheta_ref + Pdelta_ref + Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeAR),'weights',weights))

mdl4Alt.Rsquared.Ordinary

mdlA = mdl4Alt;%mdlARrestart;
mdlA.Rsquared.Ordinary;
clear temp
[ARieF_pred,ARieF_predCI] = predict(mdlA,WStbl);
[x,y,t,AUC_WS,~] = perfcurve(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>0.5),1*(ARieF_pred(~WStbl.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_WS

% Total = Wake minus Sleep
% Wake - Total = Sleep;



if 0
    %lose ~6% Rsq without ref:
    temp = WStbl(:,{'WakeNoAR','Pbeta','Palpha','Ptheta','Pdelta'});
    temp.Properties.VariableNames = {'WakeNoAR','Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'};
    mdl4Alt_noref = compact(fitglm(temp,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeAR),'weights',weights))
    mdl4Alt_noref.Rsquared.Ordinary
end
%
% mdl4Alt = compact(fitglm(WStbl,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',ExcludeAR,'weights',weights))
% mdl4Alt.Rsquared.Ordinary

%% Try multinomial regression

WStbl.NREMstate = NaN*WStbl.Epochs;
    WStbl.NREMstate(WStbl.Epochs==4) = 0;
    WStbl.NREMstate(WStbl.Epochs==2) = 1;
    WStbl.NREMstate(WStbl.Epochs==1) = 2;
    WStbl.NREMstate(WStbl.Epochs==0) = 3;

% Weights
balancefactor=1%0.86
balance = nanmean(WStbl.WakeNoAR(~WStbl.ExcludeAR));
weights = zeros(length(WStbl.WakeNoAR),1);
weights(WStbl.WakeNoAR>0.5)=1-balance*balancefactor;
weights(WStbl.WakeNoAR<=0.5)=balance*balancefactor;
sum(weights(WStbl.WakeNoAR>0.5&~WStbl.ExcludeAR));
sum(weights(WStbl.WakeNoAR<=0.5&~WStbl.ExcludeAR));

mdl4Alt = compact(fitglm(WStbl,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeAR),'weights',weights))
mdl4Alt.Rsquared.Ordinary

Amatrix = [WStbl.Pbeta_ref WStbl.Palpha_ref WStbl.Ptheta_ref WStbl.Pdelta_ref WStbl.Ptheta_ref.*WStbl.Pdelta_ref  WStbl.Ptheta_ref.*WStbl.Ptheta_ref  WStbl.Pdelta_ref.*WStbl.Pdelta_ref];
Ymatrix = WStbl.WakeNoAR;
Amatrix(logical(WStbl.ExcludeAR),:)=[];
Ymatrix(logical(WStbl.ExcludeAR))=[];
W = weights;
W(logical(WStbl.ExcludeAR))=[];

if 0
    [~,~,c]=glmfit(Amatrix,Ymatrix,'binomial','weights',W); %,'weights',weights
    PredY = glmval(c.beta,Amatrix,'logit');
    SSE = nansum(W.*(PredY - Ymatrix).^2)/nansum(W);  
SStot = nansum(W.*(PredY - nanmean(PredY)).^2)/nansum(W);  
1 - SSE/SStot;
end

%not quite:

YmatrixOrdinal = ordinal(1-Ymatrix,{'1','2','3','4'},[],[-Inf,0.5,1.5,2.5,Inf]);

[~,~,c1] = mnrfit(Amatrix(1:50:end,:),YmatrixOrdinal(1:50:end),'model','ordinal');
c1.beta


mdl4Alt


YmatrixState = WStbl.NREMstate;
YmatrixState(logical(WStbl.ExcludeAR),:)=[];
YmatrixStateOrdinal = ordinal(YmatrixState,{'1','2','3','4'},[],[-Inf,0.5,1.5,2.5,Inf]);

[~,~,c_] = mnrfit(Amatrix(1:50:end,:),YmatrixStateOrdinal(1:50:end),'model','ordinal');
c_.beta

[~,~,c1_] = mnrfit2(Amatrix(1:50:end,:),YmatrixStateOrdinal(1:50:end),W(1:50:end));
c1_.beta

[~,~,c2_] = mnrfit2(Amatrix(1:50:end,:),YmatrixOrdinal(1:50:end),W(1:50:end));
c2_.beta

[~,~,c3_] = glmfit(Amatrix(1:50:end,:),YmatrixOrdinal(1:50:end),W(1:50:end));
c2_.beta


c.beta

c.beta ./ c2_.beta

y = grp2idx(y); 

[pihat,dlow,hi] = mnrval(c.beta,Amatrix(1:5:end,:),c,'model','ordinal');

Amatrix = [WStbl.Pbeta_ref WStbl.Palpha_ref WStbl.Ptheta_ref WStbl.Pdelta_ref WStbl.Ptheta_ref.*WStbl.Pdelta_ref  WStbl.Ptheta_ref.*WStbl.Ptheta_ref  WStbl.Pdelta_ref.*WStbl.Pdelta_ref];
PredYlogodds = sum(Amatrix.*c.beta(4:end)',2) + c.beta(1)*ones(length(Amatrix),1);

c.beta = [-0.6066;1.2451;5.2324;7.152;2.016;-2.260;-4.599;-5.479;4.036;2.295];
Amatrix = [WStbl.Pbeta_ref WStbl.Palpha_ref WStbl.Ptheta_ref WStbl.Pdelta_ref WStbl.Ptheta_ref.*WStbl.Pdelta_ref  WStbl.Ptheta_ref.*WStbl.Ptheta_ref  WStbl.Pdelta_ref.*WStbl.Pdelta_ref];
PredYlogodds = sum(Amatrix.*c.beta(4:end)',2) + c.beta(1)*ones(length(Amatrix),1);

performanceOrd = PredictiveValue(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>0.5),1*(PredYlogodds(~WStbl.ExcludeAR)>0),WStbl.WakeNoAR(~WStbl.ExcludeAR))
sensspecave=0.5*(performance.Spec_sem_chance_p(1)+performance.Sens_sem_chance_p(1))
performance

% Coefficients: full sample, weighted per wake/sleep
%       Value Std. Error t value
% Pb   -7.152    0.01683 -424.84
% Pa   -2.016    0.01455 -138.57
% Pt    2.260    0.05661   39.92
% Pd    4.599    0.05202   88.40
% PtPd  5.479    0.06362   86.12
% Pt2  -4.036    0.04156  -97.09
% Pd2  -2.295    0.03060  -75.00
% 
% Intercepts:
%     Value     Std. Error t value  
% 0|1   -0.6066    0.0312   -19.4660
% 1|2    1.2451    0.0313    39.7952
% 2|3    5.2324    0.0319   163.7959
% 
% Residual Deviance: 875124.37 
% AIC: 875144.37 


% limited smple, ds50
% Coefficients:
%       Value Std. Error t value
% Pb   -6.484    0.06763 -95.872
% Pa   -1.696    0.05943 -28.536
% Pt    2.011    0.23389   8.596
% Pd    3.764    0.21347  17.631
% PtPd  4.713    0.25363  18.583
% Pt2  -3.569    0.16715 -21.351
% Pd2  -1.753    0.12331 -14.213
% 
% Intercepts:
%     Value    Std. Error t value 
% 0|1  -1.3843   0.1267   -10.9240
% 1|2   0.8851   0.1271     6.9662
% 2|3   4.9045   0.1294    37.8979
% 
% Residual Deviance: 52510.29 
% AIC: 52530.29 
% >    test <- polr(factor(Ystate) ~Pb+Pa+Pt+Pd+PtPd+Pt2+Pd2 , weights = W,data = SWdata)
% non-integer #successes in a binomial glm!glm.fit: fitted probabilities numerically 0 or 1 occurred
% >    summary(test)
% 
% Re-fitting to get Hessian
% 
% Call:
% polr(formula = factor(Ystate) ~ Pb + Pa + Pt + Pd + PtPd + Pt2 + 
%     Pd2, data = SWdata, weights = W)
% 
% Coefficients:
%       Value Std. Error t value
% Pb   -7.176     0.1194 -60.111
% Pa   -1.916     0.1025 -18.703
% Pt    2.061     0.4013   5.135
% Pd    4.606     0.3635  12.671
% PtPd  5.655     0.4418  12.800
% Pt2  -4.082     0.2905 -14.052
% Pd2  -2.376     0.2129 -11.163
% 
% Intercepts:
%     Value    Std. Error t value 
% 0|1  -0.6842   0.2181    -3.1367
% 1|2   1.1636   0.2189     5.3152
% 2|3   5.1414   0.2235    23.0063
% 
% Residual Deviance: 17567.05 
% AIC: 17587.05 

% all data
% Call:
% polr(formula = factor(Ystate) ~ Pb + Pa + Pt + Pd + PtPd + Pt2 + 
%     Pd2, data = SWdata, weights = W)
% 
% Coefficients:
%       Value Std. Error t value
% Pb   -7.152    0.01683 -424.84
% Pa   -2.016    0.01455 -138.57
% Pt    2.260    0.05661   39.92
% Pd    4.599    0.05202   88.40
% PtPd  5.479    0.06362   86.12
% Pt2  -4.036    0.04156  -97.09
% Pd2  -2.295    0.03060  -75.00
% 
% Intercepts:
%     Value     Std. Error t value  
% 0|1   -0.6066    0.0312   -19.4660
% 1|2    1.2451    0.0313    39.7952
% 2|3    5.2324    0.0319   163.7959
% 
% Residual Deviance: 875124.37 
% AIC: 875144.37 

%% STEP 6: Performance W only
% Incl AR
% clear temp
% [ARieF_pred,ARieF_predCI] = predict(mdlA,WStbl);
% temp = [WStbl.EventsAr ARieF_pred];
% thres = 0.5;
% performance = PredictiveValue(1*(WStbl.EventsAr>thres),1*(ARieF_pred>thres),WStbl.EventsAr)
% WStbl.ARieF_predModel = ARieF_pred;

% Excl AR
mdlA = mdl4Alt;
%mdlA = mdlX;
%mdlA = mdlWSC;
mdlA.Rsquared.Ordinary;
clear temp
[ARieF_pred,ARieF_predCI] = predict(mdlA,WStbl);

[x,y,t,AUC_WS,~] = perfcurve(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>0.5),1*(ARieF_pred(~WStbl.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_WS

WStbl.SWSnotN1W = 1*(WStbl.Epochs==0);
WStbl.SWSnotN1W(isnan(WStbl.Epochs)) = NaN;
WStbl.SWSnotN1W(WStbl.ExcludeAR==1) = NaN;
WStbl.SWSnotN1W(WStbl.Epochs==1) = NaN;

[x,y,t,AUC_SWS,~] = perfcurve(1*(WStbl.SWSnotN1W(~isnan(WStbl.SWSnotN1W))>0.5),ARieF_pred(~isnan(WStbl.SWSnotN1W)),0); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFindSWS=mean(t(I:(I+1)))
AUC_SWS

if 1
    thres = 0.5;
else
    thres=thresFind;
end
performance = PredictiveValue(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>0.5),1*(ARieF_pred(~WStbl.ExcludeAR)>thres),WStbl.WakeNoAR(~WStbl.ExcludeAR))
sensspecave=0.5*(performance.Spec_sem_chance_p(1)+performance.Sens_sem_chance_p(1))

WStbl.ARieF_predModel = ARieF_pred;
WStbl.ARieF_pred_logit = logit(WStbl.ARieF_predModel);
%WStbl.ARieF_pred_logit = (WStbl.ARieF_predModel);

%adjust score so threshold is 0.5
if 0
    WStbl.ARieF_pred_logit = logit(WStbl.ARieF_predModel) - logit(thres);
    WStbl.ARieF_predModel = logitinverse(WStbl.ARieF_pred_logit);
end

%% STEP 6: Performance W only, test set
% Incl AR
% clear temp
% [ARieF_pred,ARieF_predCI] = predict(mdlA,WStbl);
% temp = [WStbl.EventsAr ARieF_pred];
% thres = 0.5;
% performance = PredictiveValue(1*(WStbl.EventsAr>thres),1*(ARieF_pred>thres),WStbl.EventsAr)
% WStbl.ARieF_predModel = ARieF_pred;

% Excl AR
[ARieF_pred,ARieF_predCI] = predict(mdlA,WStblTest); 

[x,y,t,AUC_WS,~] = perfcurve(1*(WStblTest.WakeNoAR(~WStblTest.ExcludeAR)>0.5),1*(ARieF_pred(~WStblTest.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_WS

WStblTest.SWSnotN1W = 1*(WStblTest.Epochs==0);
WStblTest.SWSnotN1W(isnan(WStblTest.Epochs)) = NaN;
WStblTest.SWSnotN1W(WStblTest.ExcludeAR==1) = NaN;
WStblTest.SWSnotN1W(WStblTest.Epochs==1) = NaN;

[x,y,t,AUC_SWS,~] = perfcurve(1*(WStblTest.SWSnotN1W(~isnan(WStblTest.SWSnotN1W))>0.5),ARieF_pred(~isnan(WStblTest.SWSnotN1W)),0); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFindSWS=mean(t(I:(I+1)))
AUC_SWS

if 1
    thres = 0.5;
else
    thres=thresFind;
end
performance_ = PredictiveValue(1*(WStblTest.WakeNoAR(~WStblTest.ExcludeAR)>0.5),1*(ARieF_pred(~WStblTest.ExcludeAR)>thres),WStblTest.WakeNoAR(~WStblTest.ExcludeAR))
sensspecave_=0.5*(performance.Spec_sem_chance_p(1)+performance.Sens_sem_chance_p(1))

WStblTest.ARieF_predModel = ARieF_pred;
WStblTest.ARieF_pred_logit = logit(WStblTest.ARieF_predModel);
%WStbl.ARieF_pred_logit = (WStbl.ARieF_predModel);

%adjust score so threshold is 0.5
if 0
    WStbl.ARieF_pred_logit = logit(WStbl.ARieF_predModel) - logit(thres);
    WStbl.ARieF_predModel = logitinverse(WStbl.ARieF_pred_logit);
end

%% SKIP: Younes Lookup Table Method, acc=83.94
if 1
    PbetaCut = prctile(WStbl.Pbeta(~WStbl.ExcludeAR),10:10:90);
    PalphaCut = prctile(WStbl.Palpha(~WStbl.ExcludeAR),10:10:90);
    PthetaCut = prctile(WStbl.Ptheta(~WStbl.ExcludeAR),10:10:90);
    PdeltaCut = prctile(WStbl.Pdelta(~WStbl.ExcludeAR),10:10:90);
    
    IDa = [sum(WStbl.Pdelta>PdeltaCut,2), ...
        sum(WStbl.Ptheta>PthetaCut,2), ...
        sum(WStbl.Palpha>PalphaCut,2), ...
        sum(WStbl.Pbeta>PbetaCut,2)];
    
    ID = IDa(:,1)*1000 + IDa(:,2)*100 + IDa(:,3)*10 + IDa(:,4)*1;
    IDs=(0:9999)';
    num=NaN*ones(10000,1);
    den=NaN*ones(10000,1);
    Fnumden=NaN*ones(10000,1);
%     for i=1:10000
%         I=find(ID==IDs(i)&~WStbl.ExcludeAR);
%         den(i)=length(I);
%         num(i)=sum(WStbl.WakeNoAR(I));
%         %num2(i)=sum(WStbl.ARieF_predModel(I));
%     end
    for i=1:10000
    I=find(ID==IDs(i)&~WStbl.ExcludeAR);
    den(i)=length(I);
    if 0
        num(i)=sum(WStbl.WakeNoAR(I));
        
    else
        Fnumden(i) = nansum(WStbl.WakeNoAR(I).*weights(I)) ./ nansum(weights(I));
        num(i)=den(i).*Fnumden(i);
    end
    num2(i)=sum(WStbl.ARieF_predModel(I)); %for plots only, comparison with model
%     betaID(i) = nanmedian(WStbl.Pbeta_ref(I)); %for plots only
%     alphaID(i) = nanmedian(WStbl.Palpha_ref(I)); %for plots only
%     thetaID(i) = nanmedian(WStbl.Ptheta_ref(I)); %for plots only
%     deltaID(i) = nanmedian(WStbl.Pdelta_ref(I)); %for plots only
    end
    
PrWID_unbalanced = num./den; %no correction for missing lookups
PrWID_unbalanced(den<20)=NaN;

%balanced probabilities (i.e. modified based on prevalence so that cutoff for wake sleep is at Pr=0.5)
% WonS_balance =  1/(1./balance - 1);
% SonW = 1./PrWID_unbalanced - 1;
% SonW_balanced = SonW*WonS_balance;
% PrWID_balanced = 1./(1 + SonW_balanced);
% PrWID = PrWID_balanced;

%select balanced for further analysis
PrWID = PrWID_unbalanced; %%new

    %raw (unbalanced) probability
%     PrWID_unbalanced = num./den; %no correction for missing lookups
%     PrWID_unbalanced(den<20)=NaN;
    
    %balanced probabilities (i.e. modified based on prevalence so that cutoff for wake sleep is at Pr=0.5)
%     WonS_balance =  1/(1./balance - 1);
%     SonW = 1./PrWID_unbalanced - 1;
%     SonW_balanced = SonW*WonS_balance;
%     PrWID_balanced = 1./(1 + SonW_balanced);
    
    %select balanced for further analysis
%     PrWID = PrWID_balanced;
    
    Tx = table(IDs,PrWID,den);
    
    %PrWIDmodel = num2./den; %no correction for missing lookups
    %PrWIDmodel(den<10)=NaN;
    
    IDsa=[floor(rem(IDs,10000)/1000) floor(rem(IDs,1000)/100) floor(rem(IDs,100)/10) floor(rem(IDs,10)/1)];
    IDdelta=IDsa(:,1);
    IDtheta=IDsa(:,2);
    IDalpha=IDsa(:,3);
    IDbeta=IDsa(:,4);
    
    % Performance
    WStbl.ARieF_predID_MY = PrWID(ID+1);
    ARieF_pred = WStbl.ARieF_predID_MY;
    clear temp
    thres = 0.5;
    performanceIDMY = PredictiveValue(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>thres),1*(ARieF_pred(~WStbl.ExcludeAR)>thres),WStbl.WakeNoAR(~WStbl.ExcludeAR))
    
    %Model error
    PrWIDmodelSSE = nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum(num2)
    PrWIDmodelRsq = 1-nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum((num2.*(PrWID-nanmean(PrWID)).^2))

    [x,y,t,AUC_WSMY,~] = perfcurve(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>thres),1*(ARieF_pred(~WStbl.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
    

end

%% STEP 7: RERUN_THETA Lookup Table Method, four, normalized by ref (acc=92.5p)
PdeltaCut = prctile(WStbl.Pdelta_ref(~WStbl.ExcludeAR),10:10:90);
PthetaCut = prctile(WStbl.Ptheta_ref(~WStbl.ExcludeAR),10:10:90);
PalphaCut = prctile(WStbl.Palpha_ref(~WStbl.ExcludeAR),10:10:90);
PbetaCut = prctile(WStbl.Pbeta_ref(~WStbl.ExcludeAR),10:10:90);

IDa = [         sum(WStbl.Pdelta_ref>PdeltaCut,2), ...
    sum(WStbl.Ptheta_ref>PthetaCut,2), ...
    sum(WStbl.Palpha_ref>PalphaCut,2), ...
    sum(WStbl.Pbeta_ref>PbetaCut,2)];

ID = IDa(:,1)*1000 + IDa(:,2)*100 + IDa(:,3)*10 + IDa(:,4)*1;
IDs=(0:9999)';

IDsa=[floor(rem(IDs,10000)/1000) floor(rem(IDs,1000)/100) floor(rem(IDs,100)/10) floor(rem(IDs,10)/1)];

num=NaN*ones(10000,1);
num2=NaN*ones(10000,1);
Fnumden=NaN*ones(10000,1);
den=NaN*ones(10000,1);
betaID = NaN*ones(10000,1);
alphaID = NaN*ones(10000,1);
thetaID = NaN*ones(10000,1);
deltaID = NaN*ones(10000,1);
for i=1:10000
    I=find(ID==IDs(i)&~WStbl.ExcludeAR);
    den(i)=length(I);
    if 0
        num(i)=sum(WStbl.WakeNoAR(I));
        
    else
        Fnumden(i) = nansum(WStbl.WakeNoAR(I).*weights(I)) ./ nansum(weights(I));
        num(i)=den(i).*Fnumden(i);
    end
    num2(i)=sum(WStbl.ARieF_predModel(I)); %for plots only, comparison with model
    betaID(i) = nanmedian(WStbl.Pbeta_ref(I)); %for plots only
    alphaID(i) = nanmedian(WStbl.Palpha_ref(I)); %for plots only
    thetaID(i) = nanmedian(WStbl.Ptheta_ref(I)); %for plots only
    deltaID(i) = nanmedian(WStbl.Pdelta_ref(I)); %for plots only
end

% Keep
DecilesLookupRef.PdeltaCut = PdeltaCut;
DecilesLookupRef.PthetaCut = PthetaCut;
DecilesLookupRef.PalphaCut = PalphaCut;
DecilesLookupRef.PbetaCut = PbetaCut;
DecilesLookupRef.Pr = num./den;
DecilesLookupRef.N = den;
DecilesLookupRef.ID = IDsa;
DecilesLookupRef.IDdelta=IDsa(:,1);
DecilesLookupRef.IDtheta=IDsa(:,2);
DecilesLookupRef.IDalpha=IDsa(:,3);
DecilesLookupRef.IDbeta=IDsa(:,4);

if 0
    save('DecilesLookupRef','DecilesLookupRef','-v7.3')
end

%raw (unbalanced) probability
PrWID_unbalanced = num./den; %no correction for missing lookups
PrWID_unbalanced(den<20)=NaN;

%balanced probabilities (i.e. modified based on prevalence so that cutoff for wake sleep is at Pr=0.5)
% WonS_balance =  1/(1./balance - 1);
% SonW = 1./PrWID_unbalanced - 1;
% SonW_balanced = SonW*WonS_balance;
% PrWID_balanced = 1./(1 + SonW_balanced);
% PrWID = PrWID_balanced;

%select balanced for further analysis
PrWID = PrWID_unbalanced; %%new

NemptyID = sum(isnan(PrWID))
Tx = table(IDs,PrWID,den);

PrWIDmodel = num2./den; %no correction for missing lookups
%PrWIDmodel(den<10)=NaN;

IDsa=[floor(rem(IDs,10000)/1000) floor(rem(IDs,1000)/100) floor(rem(IDs,100)/10) floor(rem(IDs,10)/1)];
IDdelta=IDsa(:,1);
IDtheta=IDsa(:,2);
IDalpha=IDsa(:,3);
IDbeta=IDsa(:,4);

% Performance
WStbl.ARieF_predID = PrWID(ID+1);
ARieF_pred = WStbl.ARieF_predID;
clear temp
thres = 0.5;
performanceID = PredictiveValue(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>thres),1*(ARieF_pred(~WStbl.ExcludeAR)>thres),WStbl.WakeNoAR(~WStbl.ExcludeAR))

[x,y,t,AUC_WSID,~] = perfcurve(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>thres),1*(ARieF_pred(~WStbl.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
    

%ERRORS
PrWIDmodelSSE = nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum(num2)
PrWIDmodelRsq = 1-nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum((num2.*(PrWID-nanmean(PrWID)).^2))

%% STEP 8: RERUN_THETA Plots
%wake vs beta, 9 levels of delta; two levels of theta; three levels of theta
y=[2 7]; %alpha
x=[5]; %theta
figure(1); clf(1);
pauselength = 0.05;
datacolor=[0.75 0.6 0.6];
for i=1:length(y)
    for w=[0:9]
        subplot(1,length(y),i);
        I = find(IDtheta==x&IDdelta==w&IDalpha==y(i));
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I); PrWIDplot2(isnan(PrWIDplot))=NaN;
        figure(1);
        %color_ = [0.5 0.5 z/10];
        %xvals = 0:9;
        xvals = betaID(I);
        plot(xvals,100*[PrWIDplot],'-','color',datacolor,'linewidth',2); hold('on');
        plot(xvals,100*[PrWIDplot2],'k-','linewidth',2); hold('on');
        pause(pauselength);
        box('off');
        set(gca,'tickdir','out');
    end
    %xlabel('beta rank');
    xlabel('beta power');
    title(['alpha rank = ' num2str(y(i))]);
    chH = get(gca,'Children');
    neworder = [1:2:length(chH) 2:2:length(chH)];
    set(gca,'Children',chH(neworder));
    xlim([0.4 1.7]);
end
set(gcf,'position',[81   669   533   237])


%wake vs beta, 9 levels of delta; mid range level of alpha; three levels of theta
y=[5]; %alpha
x=[2 7]; %theta
figure(2); clf(2);
for i=1:length(x)
    for w=[0:9]
        subplot(1,length(x),i);
        I = find(IDtheta==x(i)&IDdelta==w&IDalpha==y);
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I); PrWIDplot2(isnan(PrWIDplot))=NaN;
        figure(2);
        %xvals = 0:9;
        xvals = betaID(I);
        plot(xvals,100*[PrWIDplot],'-','color',datacolor,'linewidth',2); hold('on');
        plot(xvals,100*[PrWIDplot2],'k-','linewidth',2); hold('on');
        pause(pauselength);
        box('off');
        set(gca,'tickdir','out');
    end
    %xlabel('beta rank');
    xlabel('beta power');
    title(['theta rank = ' num2str(x(i))]);
    
    chH = get(gca,'Children');
    neworder = [1:2:length(chH) 2:2:length(chH)];
    set(gca,'Children',chH(neworder));
    xlim([0.4 1.7]);
end
set(gcf,'position',[81   669   533   237])
%wake vs theta, 9 levels of beta; mid alpha, three levels of delta;
y=[5]; %alpha
w=[2 7]; %delta
figure(3); clf(3);
for i=1:length(w)
    for z=[0:9]
        subplot(1,length(w),i);
        I = find(IDbeta==z&IDdelta==w(i)&IDalpha==y);
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I); PrWIDplot2(isnan(PrWIDplot))=NaN;
        %color_ = [0.5 0.5 z/10];
        %xvals = 0:9;
        xvals = thetaID(I);
        plot(xvals,100*[PrWIDplot],'-','color',datacolor,'linewidth',2); hold('on');
        plot(xvals,100*[PrWIDplot2],'k-','linewidth',2); hold('on');
        pause(pauselength);
        box('off');
        set(gca,'tickdir','out');
    end
    xlabel('theta power');
    title(['delta rank = ' num2str(w(i))]);
    
    chH = get(gca,'Children');
    neworder = [1:2:length(chH) 2:2:length(chH)];
    set(gca,'Children',chH(neworder));
    if i==1
    xlim([0.6 1.6]);
    else
    xlim([0.8 1.8]);    
    end
end
set(gcf,'position',[81   669   533   237])
%wake vs delta, 9 levels of beta; mid alpha, three levels of theta;
y=[5]; %alpha
x=[2 7]; %theta
figure(4); clf(4);
for i=1:length(x)
    for z=[0:9]
        subplot(1,length(w),i);
        I = find(IDbeta==z&IDtheta==x(i)&IDalpha==y);
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I);
        PrWIDplot2(isnan(PrWIDplot))=NaN;
        figure(4);
        %color_ = [0.5 0.5 z/10];
        %xvals = 0:9;
        xvals = deltaID(I);
        plot(xvals,100*[PrWIDplot],'-','color',datacolor,'linewidth',2); hold('on');
        plot(xvals,100*[PrWIDplot2],'k-','linewidth',2); hold('on');
        pause(pauselength);
        box('off');
        set(gca,'tickdir','out');
    end
    xlabel('delta power');
    title(['theta rank = ' num2str(x(i))]);
    
    chH = get(gca,'Children');
    neworder = [1:2:length(chH) 2:2:length(chH)];
    set(gca,'Children',chH(neworder));
    xlim([0.8 2.4]);
end
set(gcf,'position',[81   669   533   237])

%% plot3d, individual data against beta, alpha, delta

figure(1); clf(1);
ds=50;
maxballs = 2000;
plotsize=10;
PrctileTheta=50;
Pthetarefrange=5
PthetarefL=prctile(WStbl.Ptheta_ref,PrctileTheta-5);
PthetarefU=prctile(WStbl.Ptheta_ref,PrctileTheta+5);
I = find(WStbl.WakeNoAR==1 & ~WStbl.ExcludeAR & WStbl.Ptheta_ref>PthetarefL & WStbl.Ptheta_ref<=PthetarefU);
I=I(1:ds:end);
h1=scatter3(WStbl.Pdelta_ref(I),WStbl.Palpha_ref(I),WStbl.Pbeta_ref(I),plotsize,'filled')
hold('on');
axis([0.5 2.2 0 3 0 3]);
ballsize=0.4;
expandview = [1 1 1.25];
temp = [diff(get(gca,'xlim')) diff(get(gca,'ylim')) diff(get(gca,'zlim'))]./expandview;
daspect(temp)
if 0
    temp = temp/mean(temp);
end
compressxyz = 1./temp;
if length(I)>maxballs
    I=I(1:floor(length(I)/maxballs):end);
    I(maxballs+1:end)=[];
end
C = repmat([3 169 244]/255,length(I),1);
h=scatter3sph(WStbl.Pdelta_ref(I),WStbl.Palpha_ref(I),WStbl.Pbeta_ref(I),'size',ballsize,'color',C,'transp',1,'compress',compressxyz);
  
%I = find(WStbl.EventsAr<0.1,Npoints);
Pthetarefin = prctile(WStbl.Ptheta_ref,PrctileTheta);
I = find(WStbl.WakeNoAR==0 & ~WStbl.ExcludeAR & abs(WStbl.Ptheta_ref-Pthetarefin)<Pthetarefrange);
I=I(1:ds:end);


%scatter3(WStbl.Pdelta_ref(I),WStbl.Palpha_ref(I),WStbl.Pbeta_ref(I),plotsize,'filled')
if length(I)>maxballs
    I=I(1:floor(length(I)/maxballs):end);
    I(maxballs+1:end)=[];
end
C = repmat([255 81 30]/255,length(I),1);
h2=scatter3sph(WStbl.Pdelta_ref(I),WStbl.Palpha_ref(I),WStbl.Pbeta_ref(I),'size',ballsize,'color',C,'transp',1,'compress',compressxyz);
set(h,'facelighting','phong','ambientstrength',0.5);
set(h2,'facelighting','phong','ambientstrength',0.5);
h3=light('position',[prctile(WStbl.Pdelta_ref(I),99.9) prctile(WStbl.Palpha_ref(I),20) prctile(WStbl.Pbeta_ref(I),50)],'style','local');
%set(gca,'CameraPosition',camerapos,'projection','perspective','gridlinestyle','-','TickLength',[0 0],'fontsize',20);
set(gca,'fontname','arial narrow');          

%scatter3(WStbl.Pdelta_ref(I),WStbl.Palpha_ref(I),WStbl.Pbeta_ref(I),plotsize,'filled')
xlabel('delta');
ylabel('alpha');
zlabel('beta');
hold('off')


%make and plot surface
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
[X,Y] = meshgrid(xlims(1):0.1:xlims(2),ylims(1):0.1:ylims(2));

N_ = length(X(:));
T_ = [X(:) Y(:) Pthetarefin*ones(N_,1) zeros(N_,1)];
T = array2table(T_);
T.Properties.VariableNames = {'Pdelta_ref','Palpha_ref','Ptheta_ref','Pbeta_ref'};
betacoef = mdlA.Coefficients.Estimate(2);
temp = -logit(predict(mdlA,T))/betacoef;
Z = reshape(temp,size(X));

hold('on')
h=surface(X,Y,Z,'linestyle','none','facecolor',[0.5 0.5 0.5],'facealpha',0.6);

delete(h1);

if 0
campos = get(gca,'camerapos');
else
campos = [14.7828  -10.6307    9.1936];
set(gca,'camerapos',campos);
end

set(gcf,'position',[     81          69        1116         837]);
set(gca,'fontsize',20)

%% repeat plot3d; test data

figure(14); clf(14);
ds=50;
maxballs = 2000;
plotsize=10;
PrctileTheta=50;
PthetarefL=prctile(WStblTest.Ptheta_ref,PrctileTheta-5);
PthetarefU=prctile(WStblTest.Ptheta_ref,PrctileTheta+5);
I = find(WStblTest.WakeNoAR==1 & ~WStblTest.ExcludeAR & WStblTest.Ptheta_ref>PthetarefL & WStblTest.Ptheta_ref<=PthetarefU);
I=I(1:ds:end);
h1=scatter3(WStblTest.Pdelta_ref(I),WStblTest.Palpha_ref(I),WStblTest.Pbeta_ref(I),plotsize,'filled')
hold('on');
axis([0.5 2.2 0 3 0 3]);
ballsize=0.4;
expandview = [1 1 1.25];
temp = [diff(get(gca,'xlim')) diff(get(gca,'ylim')) diff(get(gca,'zlim'))]./expandview;
daspect(temp)
if 0
    temp = temp/mean(temp);
end
compressxyz = 1./temp;
if length(I)>maxballs
    I=I(1:floor(length(I)/maxballs):end);
    I(maxballs+1:end)=[];
end
C = repmat([3 169 244]/255,length(I),1);
h=scatter3sph(WStblTest.Pdelta_ref(I),WStblTest.Palpha_ref(I),WStblTest.Pbeta_ref(I),'size',ballsize,'color',C,'transp',1,'compress',compressxyz);
  
%I = find(WStblTest.EventsAr<0.1,Npoints);
I = find(WStblTest.WakeNoAR==0 & ~WStblTest.ExcludeAR & abs(WStblTest.Ptheta_ref-Pthetarefin)<Pthetarefrange);
I=I(1:ds:end);
Pthetarefin = prctile(WStblTest.Ptheta_ref,PrctileTheta);

%scatter3(WStblTest.Pdelta_ref(I),WStblTest.Palpha_ref(I),WStblTest.Pbeta_ref(I),plotsize,'filled')
if length(I)>maxballs
    I=I(1:floor(length(I)/maxballs):end);
    I(maxballs+1:end)=[];
end
C = repmat([255 81 30]/255,length(I),1);
h2=scatter3sph(WStblTest.Pdelta_ref(I),WStblTest.Palpha_ref(I),WStblTest.Pbeta_ref(I),'size',ballsize,'color',C,'transp',1,'compress',compressxyz);
set(h,'facelighting','phong','ambientstrength',0.5);
set(h2,'facelighting','phong','ambientstrength',0.5);
h3=light('position',[prctile(WStblTest.Pdelta_ref(I),99.9) prctile(WStblTest.Palpha_ref(I),20) prctile(WStblTest.Pbeta_ref(I),50)],'style','local');
%set(gca,'CameraPosition',camerapos,'projection','perspective','gridlinestyle','-','TickLength',[0 0],'fontsize',20);
set(gca,'fontname','arial narrow');          

%scatter3(WStblTest.Pdelta_ref(I),WStblTest.Palpha_ref(I),WStblTest.Pbeta_ref(I),plotsize,'filled')
xlabel('delta');
ylabel('alpha');
zlabel('beta');
hold('off')


%make and plot surface
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
[X,Y] = meshgrid(xlims(1):0.1:xlims(2),ylims(1):0.1:ylims(2));

N_ = length(X(:));
T_ = [X(:) Y(:) Pthetarefin*ones(N_,1) zeros(N_,1)];
T = array2table(T_);
T.Properties.VariableNames = {'Pdelta_ref','Palpha_ref','Ptheta_ref','Pbeta_ref'};
betacoef = mdlA.Coefficients.Estimate(2);
temp = -logit(predict(mdlA,T))/betacoef;
Z = reshape(temp,size(X));

hold('on')
h=surface(X,Y,Z,'linestyle','none','facecolor',[0.5 0.5 0.5],'facealpha',0.6);

delete(h1);

set(gca,'camerapos',campos);
set(gcf,'position',[     81          69        1116         837]);
set(gca,'fontsize',20)

%% Surface4

Nlines=10;
figure(9); clf(9);

Pdelta_refX=prctile(WStbl.Pdelta_ref,1:(99-1)/9:99);
Pbeta_refX=prctile(WStbl.Pbeta_ref,1:(99-1)/9:99);
Palpha_refX=prctile(WStbl.Palpha_ref,1:(99-1)/9:99);
Ptheta_refX=prctile(WStbl.Ptheta_ref,1:(99-1)/9:99);

[Pdelta_ref,Pbeta_ref] = meshgrid(Pdelta_refX,Pbeta_refX);

Palpha_ref = Pbeta_ref*0 + prctile(WStbl.Palpha_ref,5);
Ptheta_ref = Pbeta_ref*0 + prctile(WStbl.Ptheta_ref,50);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);
figure(9)
subplot(2,2,1)
surf(Pbeta_ref,Pdelta_ref,temp)
xlabel('beta')
ylabel('delta')

Palpha_ref = Pbeta_ref*0 + prctile(WStbl.Palpha_ref,95);
Ptheta_ref = Pbeta_ref*0 + prctile(WStbl.Ptheta_ref,50);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);
figure(9)
subplot(2,2,2)
surf(Pbeta_ref,Pdelta_ref,temp)
xlabel('beta')
ylabel('delta')

Palpha_ref = Pbeta_ref*0 + prctile(WStbl.Palpha_ref,50);
Ptheta_ref = Pbeta_ref*0 + prctile(WStbl.Ptheta_ref,5);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);
figure(9)
subplot(2,2,3)
surf(Pbeta_ref,Pdelta_ref,temp)
xlabel('beta')
ylabel('delta')
%
figure(10); clf(10);
[Ptheta_ref,Pbeta_ref] = meshgrid(Ptheta_refX,Pbeta_refX);
Palpha_ref = Pbeta_ref*0 + prctile(WStbl.Palpha_ref,50);
Pdelta_ref = Pbeta_ref*0 + prctile(WStbl.Pdelta_ref,5);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);

subplot(1,2,1)
surf(Pbeta_ref,Ptheta_ref,temp)
xlabel('beta')
ylabel('theta')

[Ptheta_ref,Pbeta_ref] = meshgrid(Ptheta_refX,Pbeta_refX);
Palpha_ref = Pbeta_ref*0 + prctile(WStbl.Palpha_ref,50);
Pdelta_ref = Pbeta_ref*0 + prctile(WStbl.Pdelta_ref,95);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);

subplot(1,2,2)
surf(Pbeta_ref,Ptheta_ref,temp)
xlabel('beta')
ylabel('theta')



%% AR detect


%Exclude_ = WStbl.NoiseBinary==1 | WStbl.SpO2off==1 | WStbl.Pbeta==-Inf | isnan(WStbl.Epochs) | WStbl.Epochs<0 | WStbl.Epochs>4;
%WStbl.ExcludeAR = Exclude | WStbl.Epochs==4 & WStbl.EventsAr==0 | WStbl.Epochs~=4 & WStbl.EventsAr==1;
%WStbl.ExcludeAR = logical(WStbl.ExcludeAR);

WStbl.ExcludeW = WStbl.Epochs==4 | WStbl.Exclude;% & string(WStbl.Subj)~=string(UniqueSubjList(1:10:100));

adjustedbalance=1.685; %1.5 before adding MESA; higher number increases FP and lowers FN
balanceAR = nanmean(WStbl.EventsAr(~WStbl.ExcludeW))*adjustedbalance %adjusted to balance FP and FN using *1.5

weightsAR = 0*WStbl.EventsAr;
weightsAR(WStbl.EventsAr>0.5)=1-balanceAR;
weightsAR(WStbl.EventsAr<=0.5)=balanceAR;
b(1)=sum(weightsAR(WStbl.EventsAr>0.5 & ~WStbl.ExcludeW))
b(2)=sum(weightsAR(WStbl.EventsAr<0.5 & ~WStbl.ExcludeW))
b(1)/(sum(b))

WStbl.L1 = [NaN;WStbl.ARieF_pred_logit(1:end-1)];
WStbl.L2 = [NaN;NaN;WStbl.ARieF_pred_logit(1:end-2)];
WStbl.L3 = [NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-3)];
WStbl.L4 = [NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-4)];
WStbl.L5 = [NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-5)];
WStbl.L6 = [NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-6)];
WStbl.L7 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-7)];
WStbl.L8 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-8)];
WStbl.L9 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-9)];
WStbl.L10 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-10)];
WStbl.L11 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-11)];
WStbl.L12 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-12)];
WStbl.L13 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-13)];
WStbl.L14 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-14)];
WStbl.L15 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStbl.ARieF_pred_logit(1:end-15)];

WStbl.N1 = [WStbl.ARieF_pred_logit(2:end);NaN];
WStbl.N2 = [WStbl.ARieF_pred_logit(3:end);NaN;NaN];
WStbl.N3 = [WStbl.ARieF_pred_logit(4:end);NaN;NaN;NaN];
WStbl.N4 = [WStbl.ARieF_pred_logit(5:end);NaN;NaN;NaN;NaN];
WStbl.N5 = [WStbl.ARieF_pred_logit(6:end);NaN;NaN;NaN;NaN;NaN];
WStbl.N6 = [WStbl.ARieF_pred_logit(7:end);NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N7 = [WStbl.ARieF_pred_logit(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N8 = [WStbl.ARieF_pred_logit(9:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N9 = [WStbl.ARieF_pred_logit(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N10 = [WStbl.ARieF_pred_logit(11:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N11 = [WStbl.ARieF_pred_logit(12:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N12 = [WStbl.ARieF_pred_logit(13:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N13 = [WStbl.ARieF_pred_logit(14:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N14 = [WStbl.ARieF_pred_logit(15:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStbl.N15 = [WStbl.ARieF_pred_logit(16:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];

sum(WStbl.ExcludeW)
UniqueSubjListi = [1;(find(diff(WStbl.Subjn)>0)+1);height(WStbl)];
for i=1:length(UniqueSubjListi)-1
    WStbl.ExcludeW(UniqueSubjListi(i):UniqueSubjListi(i)+15)=1;
    WStbl.ExcludeW(UniqueSubjListi(i+1):-1:UniqueSubjListi(i+1)-15)=1;
end
sum(WStbl.ExcludeW)

%% AR detect, repeat for Test

% Exclude_ = WStblTest.NoiseBinary==1 | WStblTest.SpO2off==1 | WStblTest.Pbeta==-Inf | isnan(WStblTest.Epochs) | WStblTest.Epochs<0 | WStblTest.Epochs>4;
% %WStblTest.ExcludeAR = Exclude | WStblTest.Epochs==4 & WStblTest.EventsAr==0 | WStblTest.Epochs~=4 & WStblTest.EventsAr==1;
% %WStblTest.ExcludeAR = logical(WStblTest.ExcludeAR);

WStblTest.ExcludeW = WStblTest.Epochs==4 | WStblTest.Exclude;% & string(WStblTest.Subj)~=string(UniqueSubjList(1:10:100));
% 
% adjustedbalance=1.5; %1.5 before adding MESA
% balanceAR = nanmean(WStblTest.EventsAr(~WStblTest.ExcludeW))*adjustedbalance %adjusted to balance FP and FN using *1.5
% weightsAR = 0*WStblTest.EventsAr;
% weightsAR(WStblTest.EventsAr>0.5)=1-balanceAR;
% weightsAR(WStblTest.EventsAr<=0.5)=balanceAR;
% sum(WStblTest.EventsAr>0.5 & ~WStblTest.ExcludeW)
% sum(WStblTest.EventsAr<0.5 & ~WStblTest.ExcludeW)

WStblTest.L1 = [NaN;WStblTest.ARieF_pred_logit(1:end-1)];
WStblTest.L2 = [NaN;NaN;WStblTest.ARieF_pred_logit(1:end-2)];
WStblTest.L3 = [NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-3)];
WStblTest.L4 = [NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-4)];
WStblTest.L5 = [NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-5)];
WStblTest.L6 = [NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-6)];
WStblTest.L7 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-7)];
WStblTest.L8 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-8)];
WStblTest.L9 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-9)];
WStblTest.L10 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-10)];
WStblTest.L11 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-11)];
WStblTest.L12 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-12)];
WStblTest.L13 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-13)];
WStblTest.L14 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-14)];
WStblTest.L15 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblTest.ARieF_pred_logit(1:end-15)];

WStblTest.N1 = [WStblTest.ARieF_pred_logit(2:end);NaN];
WStblTest.N2 = [WStblTest.ARieF_pred_logit(3:end);NaN;NaN];
WStblTest.N3 = [WStblTest.ARieF_pred_logit(4:end);NaN;NaN;NaN];
WStblTest.N4 = [WStblTest.ARieF_pred_logit(5:end);NaN;NaN;NaN;NaN];
WStblTest.N5 = [WStblTest.ARieF_pred_logit(6:end);NaN;NaN;NaN;NaN;NaN];
WStblTest.N6 = [WStblTest.ARieF_pred_logit(7:end);NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N7 = [WStblTest.ARieF_pred_logit(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N8 = [WStblTest.ARieF_pred_logit(9:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N9 = [WStblTest.ARieF_pred_logit(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N10 = [WStblTest.ARieF_pred_logit(11:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N11 = [WStblTest.ARieF_pred_logit(12:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N12 = [WStblTest.ARieF_pred_logit(13:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N13 = [WStblTest.ARieF_pred_logit(14:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N14 = [WStblTest.ARieF_pred_logit(15:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblTest.N15 = [WStblTest.ARieF_pred_logit(16:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];


sum(WStblTest.ExcludeW)
UniqueSubjListi = [1;(find(diff(WStblTest.Subjn)>0)+1);height(WStblTest)];
for i=1:length(UniqueSubjListi)-1
    WStblTest.ExcludeW(UniqueSubjListi(i):UniqueSubjListi(i)+15)=1;
    WStblTest.ExcludeW(UniqueSubjListi(i+1):-1:UniqueSubjListi(i+1)-15)=1;
end
sum(WStblTest.ExcludeW)
%% mta
% 
% WStbl.ARieF_pred_logit_MTA = (mean([WStbl.ARieF_pred_logit WStbl.L1 WStbl.L2 WStbl.L3 WStbl.L4 WStbl.L5 WStbl.L6 WStbl.L7 ...
%                                     WStbl.L8 WStbl.L9 WStbl.L10 WStbl.L11 WStbl.L12 WStbl.L13 WStbl.L14 WStbl.L15 ...
%                                    WStbl.N1 WStbl.N2 WStbl.N3 WStbl.N4 WStbl.N5 WStbl.N6 WStbl.N7 ...
%                                     WStbl.N8 WStbl.N9 WStbl.N10 WStbl.N11 WStbl.N12 WStbl.N13 WStbl.N14 WStbl.N15]'))';
% 
% WStblTest.ARieF_pred_logit_MTA = (mean([WStblTest.ARieF_pred_logit WStblTest.L1 WStblTest.L2 WStblTest.L3 WStblTest.L4 WStblTest.L5 WStblTest.L6 WStblTest.L7 ...
%                                     WStblTest.L8 WStblTest.L9 WStblTest.L10 WStblTest.L11 WStblTest.L12 WStblTest.L13 WStblTest.L14 WStblTest.L15 ...
%                                    WStblTest.N1 WStblTest.N2 WStblTest.N3 WStblTest.N4 WStblTest.N5 WStblTest.N6 WStblTest.N7 ...
%                                     WStblTest.N8 WStblTest.N9 WStblTest.N10 WStblTest.N11 WStblTest.N12 WStblTest.N13 WStblTest.N14 WStblTest.N15]'))';



%%

%mdlAR = compact(fitglm(WStbl,...
%    ['EventsAr ~ ARieF_pred_logit + L9 + L8 + L7 + L6 + L5 + L4 + L3 + L2 + L1 + N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlAR = compact(fitglm(WStbl,...
    ['EventsAr ~ ARieF_pred_logit + L9 + L7 + L5 + L3 + L1 + N1 + N3 + N5 + N7 + N9 '],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
% mdlAR = compact(fitglm(WStbl,...
%     ['EventsAr ~ ARieF_pred_logit + ARieF_pred_logit_MTA'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
% mdlAR = compact(fitglm(WStbl,...
%     ['EventsAr ~ ARieF_pred_logit + ARieF_pred_logit_MTA + L9 + L7 + L5 + L3 + L1 + N1 + N3 + N5 + N7 + N9 '],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlAR.Rsquared.Ordinary

if 0 %alternative arousal model that avoids use of positive coefficients 
    mdlARx = compact(fitglm(WStbl,...
        ['EventsAr ~ ARieF_pred_logit + L9 + L7 + L5 + L3 + N5 + N7 + N9 '],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
    mdlARx = compact(fitglm(WStbl,...
        ['EventsAr ~ ARieF_pred_logit + L9 + L7 + L5 + N5 + N7 + N9 '],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
    
    mdlARx.Rsquared.Ordinary
    WStbl
    %save mdlARx mdlARx -v7.3
end

if 1
mdlARa = compact(fitglm(WStbl,...
     ['EventsAr ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARa.Rsquared.Ordinary
thresOrig = logitinverse(-mdlARa.Coefficients.Estimate(1)/mdlARa.Coefficients.Estimate(2));
end
% mdlAR_ = compact(fitglm(WStbl,...
%      ['EventsAr ~ ARieF_pred_logit + L9 + L8 + L7 + L6 + L5 + L4 + L3 + L2 + L1 + N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 '],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))
% mdlAR_.Rsquared.Ordinary
 % mdlAR_rebuild = compact(fitglm(WStbl,...
%     ['EventsAr ~ ARieF_pred_logit + L10 + L8 + L6 + L4 + L2 + N2 + N4 + N6 + N8 + N10 '],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))
%
% mdlAR_rebuild = compact(fitglm(WStbl,...
%     ['EventsAr ~ ARieF_pred_logit + L6 + L5 + L4 + L3 + L2 + L1 + N1 + N2 + N3 + N4 + N5 + N6'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))

% mdlAR_rebuild = compact(fitglm(WStbl,...
%     ['EventsAr ~ ARieF_pred_logit + SWSpred_logit + ARieF_pred_logit_last1 + ARieF_pred_logit_last3 + ARieF_pred_logit_last5 + ARieF_pred_logit_last7 + ARieF_pred_logit_last9 + ARieF_pred_logit_next1 + ARieF_pred_logit_next3  + ARieF_pred_logit_next5  + ARieF_pred_logit_next7  + ARieF_pred_logit_next9'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))

% mdlAR_rebuild = compact(fitglm(WStbl,...
%     ['ARieF2 ~ Palpha_ref*Pbeta_ref + Pdelta_ref*Ptheta_ref + Pdelta_ref^2 + Ptheta_ref^2'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR2)) %
% mdlAR_rebuild.Rsquared.Ordinary

%% remake of Ar model without theta (dropout of theta seen as cause for false arousals)


sum(WStbl.ExcludeW)
UniqueSubjListi = [1;(find(diff(WStbl.Subjn)>0)+1);height(WStbl)];
for i=1:length(UniqueSubjListi)-1
    WStbl.ExcludeW(UniqueSubjListi(i):UniqueSubjListi(i)+15)=1;
    WStbl.ExcludeW(UniqueSubjListi(i+1):-1:UniqueSubjListi(i+1)-15)=1;
end
sum(WStbl.ExcludeW)

mdlexcltheta = compact(fitglm(WStbl,['WakeNoAR ~ Pbeta_ref + Palpha_ref + Pdelta_ref'],'Distribution','binomial','Link','logit','Exclude',logical(WStbl.ExcludeAR),'weights',weights))
mdlexcltheta.Rsquared.Ordinary

clear WStblx
WStblx = WStbl(:,{'EventsAr','WakeNoAR','Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'});
WStblx.ARieF_pred_logit = logitinverse(predict(mdlexcltheta,WStbl));

mdlexclthetaxB = compact(fitglm(WStblx,['EventsAr ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlexclthetaxB.Rsquared.Ordinary

WStblx.L1 = [NaN;WStblx.ARieF_pred_logit(1:end-1)];
WStblx.L3 = [NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-3)];
WStblx.L5 = [NaN;NaN;NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-5)];
WStblx.L7 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-7)];
WStblx.L9 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-9)];
WStblx.N1 = [WStblx.ARieF_pred_logit(2:end);NaN];
WStblx.N3 = [WStblx.ARieF_pred_logit(4:end);NaN;NaN;NaN];
WStblx.N5 = [WStblx.ARieF_pred_logit(6:end);NaN;NaN;NaN;NaN;NaN];
WStblx.N7 = [WStblx.ARieF_pred_logit(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblx.N9 = [WStblx.ARieF_pred_logit(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];


mdlARexcltheta = compact(fitglm(WStblx,...
        ['EventsAr ~ ARieF_pred_logit + L9 + L7 + L5 + L3 + L1 + N1 + N3 + N5 + N7 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARexcltheta.Rsquared.Ordinary

%fit directly to Ar
mdlexclthetax = compact(fitglm(WStblx,['EventsAr ~ Pbeta_ref + Palpha_ref + Pdelta_ref'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlexclthetax.Rsquared.Ordinary
WStblx.ARieF_pred_logit = logitinverse(predict(mdlexclthetax,WStblx));

WStblx.L1 = [NaN;WStblx.ARieF_pred_logit(1:end-1)];
WStblx.L3 = [NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-3)];
WStblx.L5 = [NaN;NaN;NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-5)];
WStblx.L7 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-7)];
WStblx.L9 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;WStblx.ARieF_pred_logit(1:end-9)];
WStblx.N1 = [WStblx.ARieF_pred_logit(2:end);NaN];
WStblx.N3 = [WStblx.ARieF_pred_logit(4:end);NaN;NaN;NaN];
WStblx.N5 = [WStblx.ARieF_pred_logit(6:end);NaN;NaN;NaN;NaN;NaN];
WStblx.N7 = [WStblx.ARieF_pred_logit(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
WStblx.N9 = [WStblx.ARieF_pred_logit(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];

%refine WSinfo
%mdlARbackup = mdlAR
mdlARexclthetaX = compact(fitglm(WStblx,...
        ['EventsAr ~ ARieF_pred_logit + L9 + L7 + L5 + L3 + N3 + N5 + N7 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeW,'weights',weightsAR))
mdlARexclthetaX.Rsquared.Ordinary

mdlWSexclthetaX = compact(fitglm(WStblx,...
        ['WakeNoAR ~ ARieF_pred_logit + L9 + L7 + L5 + L3 + L1 + N1 + N3 + N5 + N7 + N9'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeAR,'weights',weights))
mdlWSexclthetaX.Rsquared.Ordinary



%%
newARpred_ = predict(mdlAR,WStbl);
WStbl.newARpredF = newARpred_;
WStbl.newARpred_logit = logit(newARpred_);
WStbl.newARpred = 1*(newARpred_>0.5);
WStbl.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceAR = PredictiveValue(1*(WStbl.EventsAr(~WStbl.ExcludeW)>0.5),1*(newARpred_(~WStbl.ExcludeW)>thres),WStbl.EventsAr(~WStbl.ExcludeW))
%performanceX = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(WStbl.ARieF_pred_logit(~ExcludeW)>-1.2954),WStbl.EventsAr(~ExcludeW))
%performanceAR2 = PredictiveValue(1*(WStbl.ARieF2(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),WStbl.ARieF2(~ExcludeW))

%% contrast with use of original WS function

performanceAR_WS = PredictiveValue(1*(WStbl.EventsAr(~WStbl.ExcludeW)>0.5),1*(WStbl.ARieF_pred_logit(~WStbl.ExcludeW)>logit(thresOrig)),WStbl.EventsAr(~WStbl.ExcludeW))
%performanceX = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(WStbl.ARieF_pred_logit(~ExcludeW)>-1.2954),WStbl.EventsAr(~ExcludeW))
%performanceAR2 = PredictiveValue(1*(WStbl.ARieF2(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),WStbl.ARieF2(~ExcludeW))

%% Plot AR constants themselves
figure(210); clf(210);
set(gcf,'color',[1 1 1]);

mdlARtable = mdlAR.Coefficients;
mdlARtable = mdlARtable([1 7:-1:3 2 8:12],:)
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
newARpred_ = predict(mdlAR,WStblTest);
WStblTest.newARpredF = newARpred_;
WStblTest.newARpred_logit = logit(newARpred_);
WStblTest.newARpred = 1*(newARpred_>0.5);
WStblTest.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceARTest = PredictiveValue(1*(WStblTest.EventsAr(~WStblTest.ExcludeW)>0.5),1*(WStblTest.newARpredF(~WStblTest.ExcludeW)>thres),WStblTest.EventsAr(~WStblTest.ExcludeW))

performanceX = PredictiveValue(1*(WStblTest.EventsAr(~WStblTest.ExcludeW)>0.5),1*(WStblTest.ARieF_pred_logit(~WStblTest.ExcludeW)>logit(thresOrig)),WStblTest.EventsAr(~WStblTest.ExcludeW))
%performanceAR2 = PredictiveValue(1*(WStbl.ARieF2(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),WStbl.ARieF2(~ExcludeW))

%%
% I1 = find(WStbl.newARpredF>0.5&WStbl.EventsAr<0.5&~ExcludeW);
% I2 = find(WStbl.newARpredF<0.5&WStbl.EventsAr>0.5&~ExcludeW);
% I3 = find(WStbl.newARpredF<0.5&WStbl.EventsAr<0.5&~ExcludeW);
% I4 = find(WStbl.newARpredF>0.5&WStbl.EventsAr>0.5&~ExcludeW);
% I5 = find(WStbl.newARpredF>0.5&~ExcludeW);
% I6 = find(WStbl.EventsAr>0.5&~ExcludeW);

[x,y,t,AUC_AR,~] = perfcurve(1*(WStbl.EventsAr(~WStbl.ExcludeW)>0.5),1*(WStbl.newARpredF(~WStbl.ExcludeW)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_AR


%% combined WakeSleepArousal

WStbl.WAS = max([WStbl.newARpred_logit WStbl.ARieF_pred_logit]')';
WStbl.WASPr = logitinverse(WStbl.WAS);
thres=0.5;

performanceWAS = PredictiveValue(1*(WStbl.EventsAr(~WStbl.Exclude)>0.5),1*(WStbl.WASPr(~WStbl.Exclude)>thres),WStbl.EventsAr(~WStbl.Exclude))
%performanceX = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(WStbl.ARieF_pred_logit(~ExcludeW)>-1.2954),WStbl.EventsAr(~ExcludeW))
%performanceAR2 = PredictiveValue(1*(WStbl.ARieF2(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),WStbl.ARieF2(~ExcludeW))

WStblTest.WAS = max([WStblTest.newARpred_logit WStblTest.ARieF_pred_logit]')';
WStblTest.WASPr = logitinverse(WStblTest.WAS);
thres=0.5;

performanceWASTest = PredictiveValue(1*(WStblTest.EventsAr(~WStblTest.Exclude)>0.5),1*(WStblTest.WASPr(~WStblTest.Exclude)>thres),WStblTest.EventsAr(~WStblTest.Exclude))


%% AR hist
% mdl6=mdlAR_rebuild;
% %[mdl5x.Coefficients.Estimate(2:end) mdlA.Coefficients.Estimate(2:end)]./sum(abs([mdl5x.Coefficients.Estimate(2:end) mdlA.Coefficients.Estimate(2:end)]))
%
%
% Imain = find(string(mdl5a2.Coefficients.Properties.RowNames)=="ARieF_pred_logit");
% coefmain = mdl6.Coefficients.Estimate(Imain);
%
%
% newARpred_ = predict(mdl6,WStbl);
% %newARpred_ = predict(mdl5a_beta,WStbl);
%
% WStbl.newARpredF = newARpred_;
% WStbl.newARpred_logit = logit(newARpred_);
% WStbl.newARpred = 1*(newARpred_>0.5);
%     WStbl.newARpred(isnan(newARpred_)) = NaN;
%
% thres=0.5;
% performance = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),WStbl.EventsAr(~ExcludeW))
% %performanceX = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(WStbl.ARieF_pred_logit(~ExcludeW)>-1.2954),WStbl.EventsAr(~ExcludeW))

%plot histogram
figure(126); clf(126);
dStep=0.1;
Centers=-10:dStep:10;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
Edges(end)=Inf;
ax12(1)=subplot(2,1,1);
[h7,edges] = histcounts(logit(WStbl.newARpredF(WStbl.EventsAr>0.99 & ~WStbl.ExcludeW)),Edges);
[h8,edges] = histcounts(logit(WStbl.newARpredF(WStbl.EventsAr<0.01 & ~WStbl.ExcludeW)),Edges);
if 0
area = sum([h7,h8])*dStep;
h7=h7/area;
h8=h8/area;
else
%area = sum([h7,h8])*dStep;
h7=h7/(sum([h7])*dStep);
h8=h8/(sum([h8])*dStep);  
   
end

bar(Centers,h7,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);
box('off');


if 0
ax12(1)=subplot(2,1,2);
[h7,edges] = histcounts(WStbl.ARieF_pred_logit(WStbl.EventsAr>0.99 & ~WStbl.ExcludeW),Edges);
[h8,edges] = histcounts(WStbl.ARieF_pred_logit(WStbl.EventsAr<0.01 & ~WStbl.ExcludeW),Edges);
if 0
area = sum([h7,h8])*dStep;
h7=h7/area;
h8=h8/area;
else
%area = sum([h7,h8])*dStep;
h7=h7/(sum([h7])*dStep);
h8=h8/(sum([h8])*dStep);  
   
end

bar(Centers,h7,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);
box('off');
end


%% PLOT TRACES
%figure(1); clf(1);
figure(); clf(1);
set(gcf,'color',[1 1 1]);
jj=51000; %952000 shows concordance with flow, better than ar; 456856 shows clearer diffs vss ORP (noweightsversion); 457110 false arousal
%517110 arousals weak per method; local change detection may add
% information better. Same with 641910 661910
% 651110,652110 low values in wake detection, ORP similar (depending on technique)
% 10000 is a problem for MY method;
% 221000 works well 709910 722910
% 9000 nice correlation with flow
% 34000 highlights advantage of knowing underlying model, not just Pr; also MYID fail
% concordance with flow 42000
% 170000 is great %335000
% 123000 great AR but not great W 1100500 24000 1185000
% 1204000 179500 (over)scored arousals/overshoots with minimal EEG changes
% 3000000 new (underscored)arousals 3120000, no arousals scored in REM 3140000
% , 3179000, 3182000, 129000, 151000, 34000
% difficult section: 1203950
% 179500
clear ax1

plotW=1;
plotARdelta=1;
plotMY=0;
Nsampplot = 300;
RefPlot=mode(WStbl.ref(jj:jj+Nsampplot));

ax1(1)=subplot(3,1,1);
stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.EventsAr(jj:jj+Nsampplot),'k');
hold('on');
stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.Epochs(jj:jj+Nsampplot)*0.25+1.25,'k');
box('off');
set(gca,'xticklabels',[],'xcolor',[1 1 1])
if plotMY
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.ARieF_predID(jj:jj+Nsampplot),'g');
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.ARieF_predID_MY(jj:jj+Nsampplot),'b');
end
if plotW
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.ARieF_predModel(jj:jj+Nsampplot),'r');
end
if plotARdelta
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.newARpredF(jj:jj+Nsampplot),'g');
    stairs(WStbl.Time(jj:jj+Nsampplot),(WStbl.WASPr(jj:jj+Nsampplot)),'b');
    % stairs(WStbl.Time(jj:jj+Nsampplot),logitinverse(WStbl.SWSpred_logit(jj:jj+Nsampplot)),'r:');
end
% stairs(WStbl.Time(1:jj),ARieF_predCI(1:jj,1),'color',[0.5 0.2 0.2]);
% stairs(WStbl.Time(1:jj),ARieF_predCI(1:jj,2),'color',[0.5 0.2 0.2]);
hold('off')
ylim([-0.1 2.3])
logit = @(p) log(p./(1-p));
ax1(2)=subplot(3,1,2);
stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.EventsAr(jj:jj+Nsampplot),'k');
hold('on');
box('off')
set(gca,'xticklabels',[],'xcolor',[1 1 1])
if plotMY
    stairs(WStbl.Time(jj:jj+Nsampplot),logit(WStbl.ARieF_predID(jj:jj+Nsampplot)),'g');
    stairs(WStbl.Time(jj:jj+Nsampplot),logit(WStbl.ARieF_predID_MY(jj:jj+Nsampplot)),'b');
end
if plotW
    stairs(WStbl.Time(jj:jj+Nsampplot),logit(WStbl.ARieF_predModel(jj:jj+Nsampplot)),'r');
end
if plotARdelta
    
    stairs(WStbl.Time(jj:jj+Nsampplot),logit(WStbl.newARpredF(jj:jj+Nsampplot)),'g');
    stairs(WStbl.Time(jj:jj+Nsampplot),logit(WStbl.WASPr(jj:jj+Nsampplot)),'b');
    % stairs(WStbl.Time(jj:jj+Nsampplot),(WStbl.SWSpred_logit(jj:jj+Nsampplot)),'r:');
end

%ax1(3)=subplot(3,1,3);
%stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.VI(jj:jj+Nsampplot),'g-','linewidth',1.5);
hold('on');
if 0
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.Pdelta_ref(jj:jj+Nsampplot),'k-');
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.Ptheta_ref(jj:jj+Nsampplot),'b-');
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.Palpha_ref(jj:jj+Nsampplot),'r-');
    stairs(WStbl.Time(jj:jj+Nsampplot),WStbl.Pbeta_ref(jj:jj+Nsampplot),'k-');
end
linkaxes(ax1,'x');
box('off')
%%
save workspace20190610c -v7.3



%% FrWakeSleep (~2 min)

newcol_ = {'FWakeNoAR','FWake','FARinSleep','FARieF_pred_logit','meantemp_ARieF_predModel','meantemp_ARieF_pred_logit'};
subjrowrangeI=[];
subjrowrangeN=[];
for n=1:length(UniqueSubjList)
    n
    
    
    li = UniqueSubjListi(n);
    if n<length(UniqueSubjList)
        ri = UniqueSubjListi(n+1)-1;
    else
        ri = length(UniqueSubjList)
    end
    subjrowrange = li:ri;
    
    temp_WakeNoAR = WStbl.WakeNoAR(subjrowrange);
    temp_ExcludeAR = WStbl.ExcludeAR(subjrowrange);
    temp_Wake = 1*(WStbl.Epochs(subjrowrange)==4);
    temp_Wake(isnan(WStbl.Epochs(subjrowrange)))=NaN;
    temp_Wake(WStbl.Epochs(subjrowrange)>4)=NaN;
    FWakeNoAR = nansum(temp_WakeNoAR(~temp_ExcludeAR)>0.5)/nansum(~isnan(temp_WakeNoAR(~temp_ExcludeAR)));
    temp_AR = WStbl.EventsAr(subjrowrange);
    temp_ARinSleep = temp_AR>0.5&temp_Wake==0;
    FARinSleep = nansum(temp_ARinSleep==1)/nansum(~isnan(temp_ARinSleep));
    FWake = nansum(temp_Wake>0.5)/nansum(~isnan(temp_Wake));
    
    temp_WakeModel = 1*(WStbl.Epochs(subjrowrange)==4);
    temp_ARieF_pred_logit = WStbl.ARieF_pred_logit(subjrowrange);
    FARieF_pred_logit = nansum(temp_ARieF_pred_logit>0)/nansum(~isnan(temp_ARieF_pred_logit));
    
    temp_ARieF_predModel = WStbl.ARieF_predModel(subjrowrange);
    meantemp_ARieF_predModel = nanmean(temp_ARieF_predModel);
    
    temp_ARieF_pred_logit = WStbl.ARieF_pred_logit(subjrowrange);
    meantemp_ARieF_pred_logit = nanmean(temp_ARieF_pred_logit);
    %rowrangesubset = strcmp(WStbl.Subj,SubjTemp)==1;
    %temp = prctile(tempdata,centilex);
    for j=1:length(newcol_)
        if any(strcmp(newcol_{j},WStbl.Properties.VariableNames))==0
            nrow = size(WStbl,1);
            eval(['WStbl.' newcol_{j} '=zeros(nrow,1)+NaN;']);
        end
        eval(['WStbl.' newcol_{j} '(subjrowrange)=' newcol_{j} ';']);
    end
end

%% Compare to ref

figure(121)
scatter(WStbl.FWakeNoAR(UniqueSubjListi),WStbl.Pbeta5p(UniqueSubjListi),10,'filled','markerfacealpha',0.5)
figure(122)
scatter(WStbl.FWakeNoAR(UniqueSubjListi),WStbl.ref(UniqueSubjListi),10,'filled','markerfacealpha',0.5)

figure(13)
scatter(WStbl.Pbeta5p(UniqueSubjListi),WStbl.ref(UniqueSubjListi),10,'filled','markerfacealpha',0.5)
figure(132)
scatter(WStbl.Palpha25p(UniqueSubjListi),WStbl.ref(UniqueSubjListi),10,'filled','markerfacealpha',0.5)

%%
figure(14)
%{'FWakeNoAR','FWake','FARinSleep','FARieF_pred_logit','meantemp_ARieF_predModel'};
subplot(2,2,1);
scatter(WStbl.FWake(UniqueSubjListi),WStbl.FARieF_pred_logit(UniqueSubjListi),10,'filled','markerfacealpha',0.5)

subplot(2,2,2);
scatter(WStbl.FWake(UniqueSubjListi) + WStbl.FARinSleep(UniqueSubjListi),WStbl.FARieF_pred_logit(UniqueSubjListi),10,'filled','markerfacealpha',0.5)

subplot(2,2,3);
scatter(WStbl.FWake(UniqueSubjListi) + WStbl.FARinSleep(UniqueSubjListi),WStbl.meantemp_ARieF_pred_logit(UniqueSubjListi),10,'filled','markerfacealpha',0.5)

%% Predicting the future? Works within N2
figure(16)


delta_=[14];
NbreathsInFuture_= 2 + 0*delta_;

VarX_ = WStbl.ARieF_pred_logit;
% VarX_ = logitinverse(WStbl.ARieF_predID);
% VarX_ = WStbl.ARieF_predID;
% VarX_ = WStbl.ARieF_predModel;
% VarX_ = WStbl.Epochs;

for j=1:length(NbreathsInFuture_)
    NbreathsInFuture=NbreathsInFuture_(j); %15
    delta = delta_(j);
    %WStbl.EventsAr
    
    I = find((WStbl.Epochs==1)&(WStbl.EventsAr==0)==1);
    
    %remove breaths out of range
    I(find((I + NbreathsInFuture + delta)>size(WStbl,1)))=[];
    
    %remove breaths if different subjects
    I((find(WStbl.Subjn(I)~=WStbl.Subjn(I+NbreathsInFuture))))=[];
    
    if 1 %also remove breaths where the next 1 is arousal from assessment
        I((find(WStbl.EventsAr(I+1)>0)))=[];
        I((find(WStbl.EventsAr(I+2)>0)))=[];
        %I((find(WStbl.EventsAr(I+3)==0)))=[];
    end
    
    %remove breaths if clearly beyond expected level
    I((find(WStbl.Time(I+NbreathsInFuture)>WStbl.Time(I)+300)))=[];
    
    VarX = VarX_(I);
    
    currentAR = WStbl.EventsAr(I);
    futureAR = WStbl.EventsAr(I+NbreathsInFuture);
    futureAR2 = 1*(sum(WStbl.EventsAr(I+[NbreathsInFuture:NbreathsInFuture+delta])>0.5,2)>0);
    
    PlogitCut = prctile(VarX,5:5:95); %5:5:95
    PlogitCutScore = 1+sum(VarX>PlogitCut,2);
    
    num=NaN*ones(length(PlogitCut)+1,1);
    xval=NaN*ones(length(PlogitCut)+1,1);
    den=NaN*ones(length(PlogitCut)+1,1);
    for i=1:(length(PlogitCut)+1)
        Ix=find(PlogitCutScore==i);
        den(i)=sum(~isnan(futureAR(Ix)));
        num(i)=sum(futureAR2(Ix)>0.5);
        xval(i)=median(VarX(Ix));
        %num2(i)=sum(WStbl.ARieF_predModel(I));
    end
    
    normalci_ = normalci(num./den,den);
    
    
    plot(xval,num./den,'k.-','markersize',12);
    hold('on')
    plot(xval,num./den-normalci_,'-','color',[0.5 0.5 0.5]);
    plot(xval,num./den+normalci_,'-','color',[0.5 0.5 0.5]);
    %plot(logitinverse(xval),num./den,'.-')
end


%%

if 0
    mdl1test = compact(fitglm(WStbl,['WakeNoAR ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeAR)); %'weights',weights
mdl1 = fitglme(WStbl,['WakeNoAR ~ ARieF_pred_logit + (1|Subj)'],'Link','logit','Exclude',WStbl.ExcludeAR);
mdl1test = compact(fitglm(WStbl,['WakeNoAR ~ Pbeta_ref'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeAR)); %'weights',weights
mdl1 = fitglme(WStbl,['WakeNoAR ~ Pbeta_ref + (1|Subj)'],'Link','logit','Exclude',WStbl.ExcludeAR);

%mdlx = fitglme(WStbl,'EventsAr ~ Pbeta + (1|Subj)','Link','logit')

%WStbl.ExcludeAR = ExcludeAR;

%%
%tempTable1 = WStbl(~ExcludeAR,:);
tabletemp = tempTable1(1:10:100000,{'WakeNoAR','ARieF_pred_logit','Subj','Pbeta_ref'});
tabletemp.Subj = nominal(tabletemp.Subj);
mdl1_ = fitglme(tabletemp,['WakeNoAR ~ Pbeta_ref + (1|Subj)'],'Link','logit')

mdl1 = fitglme(tempTable1,['WakeNoAR ~ ARieF_pred_logit + (1|Subj)'],'Link','logit')

%%

% mdl1 = compact(mdl1);
REs = randomEffects(mdl1);

temp = WStbl.Subj;
temp2 = ['NaN';WStbl.Subj(1:end-1)];
I=find(strcmp(temp,temp2)==0);

WStbl.Subj(I);
BreathDataFullTable_ = WStbl(I,:);
BreathDataFullTable_.REs = REs;

mdlthres = fitglm(BreathDataFullTable_,'REs ~ Pbeta10p + Pbeta5p + Palpha25p + Palpha5p')
mdlthres.Rsquared.Ordinary
WStbl.PowerRef=predict(mdlthres,WStbl);

end
%% IndividualSubject WS Data plots
clear temp
%newcol_ = 'FWakeNoAR';
%subjrowrangeI=[];
%subjrowrangeN=[];
nrows=4;
ncols=8;
clear cutoff_NW_NS
plothist1=1
for n=1:length(UniqueSubjList)
    n
    SubjTemp = UniqueSubjList(n);
    subjrowrange = WStbl.Subjn==SubjTemp;
    %temp_WakeNoAR = WStbl.WakeNoAR(strcmp(WStbl.Subj,SubjTemp)==1);
    %temp_ExcludeAR = ExcludeAR(strcmp(WStbl.Subj,SubjTemp)==1);
    balanceInd = nanmean(WStbl.WakeNoAR(~WStbl.ExcludeAR&subjrowrange==1));
    
    weightsInd = 0*WStbl.WakeNoAR(subjrowrange==1);
    IW = WStbl.WakeNoAR(subjrowrange==1)>0.5 & WStbl.ExcludeAR(subjrowrange==1)==0;
    IS = WStbl.WakeNoAR(subjrowrange==1)<=0.5 & WStbl.ExcludeAR(subjrowrange==1)==0;
    %check
    
    %     IW = logical([1 0 0])';
    %     IS = logical([0 1 1])';
    %balanceInd_ = sum(IW)/(sum(IW)+sum(IS));
    
    weightsInd(IW)=1-balanceInd;
    weightsInd(IS)=balanceInd;
    
    %check
    sum(weightsInd(IW));
    sum(weightsInd(IS));
    
    mdltemp = compact(fitglm(WStbl(subjrowrange==1,:),['WakeNoAR ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeAR(subjrowrange==1),'weights',weightsInd));
    
    cutoff_NW_NS(n,:) = [-mdltemp.Coefficients.Estimate(1)/mdltemp.Coefficients.Estimate(2) sum(IW) sum(IS) balanceInd];
    if plothist1
        %figure(99 + floor((n-1)/(nrows*ncols)))
        figure(99)
        subplot(nrows,ncols,mod(n-1,nrows*ncols)+1)
        dStep=0.25;
        Centers=-10:dStep:10;
        Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
        
        [h1,edges] = histcounts(logit(WStbl.ARieF_predModel(subjrowrange==1&WStbl.Epochs==4&WStbl.EventsAr>0.95)),Edges);
        [h2,edges] = histcounts(logit(WStbl.ARieF_predModel(subjrowrange==1&WStbl.Epochs<4&WStbl.EventsAr<0.05)),Edges);
        
        %stairs(Centers-dStep/2,h1);
        %stairs(Centers-dStep/2,h2);
        bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
        hold('on');
        bar(Centers,h2,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
        plot(cutoff_NW_NS(n,1)*[1 1],get(gca,'ylim'),'k')
        hold('off');
    end
end
%Compare to ref
%
% mdltemp1 = compact(fitglm(WStbl,['WakeNoAR ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',ExcludeAR,'weights',weights))
% cutoffalldata = -mdltemp1.Coefficients.Estimate(1)/mdltemp1.Coefficients.Estimate(2)

%% No bias due to non-inclusion of EEG amplitude; but, was there before? would need to rerun analysis using "absolute power" approach

figure(100)
%cutoff_NW_NS(:,1)

temp = WStbl.Subj;
temp2 = ['NaN';WStbl.Subj(1:end-1)];


plot(RefSubj-mean(RefSubj),cutoff_NW_NS(:,1),'.')

%% Time course of Pwake after sleep onset


VarX_ = WStbl.WSBalance;
VarX2_ = WStbl.WakeIntensity;
VarX3_ = WStbl.SleepIntensity;
% VarX_ = logitinverse(WStbl.ARieF_predID);
% VarX_ = WStbl.ARieF_predID;
% VarX_ = WStbl.ARieF_predModel;
% VarX_ = WStbl.Epochs;

Exp = 1;
switch Exp
    case 1
        NbreathsInFuture=round(60*25/dT);
        NbreathsInPast=round(15/dT);
        %WStbl.EventsAr
        
        I=find([NaN;diff(WStbl.EventsAr)]<-0.5); %find front edge of arousal (could be in wake or sleep)
        %isolate analysis to those in or near wake at sleep onset
        if 1
%             temp = 1*(sum(WStbl.Epochs(I+[-NbreathsInPast:0])==4,2)==0); %remove if no recent wake
%             I(temp==1)=[];
            
            temp = 1*(sum(WStbl.EventsAr(I+[-10:-1])<0.5,2)>0); %remove if any breaths Xback are in sleep
            I(temp==1)=[];
%             
%             temp = 1*(sum(WStbl.EventsAr(I+[1:5])>0.5,2)>0); %remove if any breaths Xfwd are in wake
%             I(temp==1)=[];
        end
        
        %WStbl.EventsAr(I(1:10)+[-3:-1])<0.5,@
        
        %remove breaths out of range
        I(find((I + NbreathsInFuture)>size(WStbl,1)))=[];
        
        %remove if any arousal in range?
        if 0
            futureAR2 = 1*(sum(WStbl.EventsAr(I+[0:NbreathsInFuture])>0.5,2)>(NbreathsInFuture*0.50));
            I(futureAR2==1)=[];
        end
        %remove if no SWS in range
        if 1
            futureSWS = 1*(sum(WStbl.Epochs(I+[0:NbreathsInFuture])==0,2)<(NbreathsInFuture*0.001));
            I(futureSWS==1)=[];
        end
    case 2
        NbreathsInFuture=60;
        NbreathsInFutureRule=60;
        NbreathsInPast=30; %plot
        NbreathsInPastRule=30;
        
        I=find([NaN;diff(WStbl.EventsAr)]<-0.5); %find front edge of arousal (could be in wake or sleep)
        %isolate analysis to those in or near wake at sleep onset
        
        %remove breaths out of range
        I(find((I + NbreathsInFuture)>size(WStbl,1)))=[];
        
        if 1
            temp = 1*(sum(WStbl.Epochs(I+[-NbreathsInPastRule:0])<4,2)==0);
            I(temp==1)=[];
        end
        
        if 1
            temp = 1*(sum(WStbl.EventsAr(I+[-NbreathsInPastRule:-10])>0.5,2)>0);
            I(temp==1)=[];
        end
        
        %any arousal in range?
        if 1
            futureAR2 = 1*(sum(WStbl.EventsAr(I+[0:NbreathsInFutureRule])>0.5,2)>0);
            I(futureAR2==1)=[];
        end
end
%I(find(WStbl.Epochs(I)<4))=[];

%     I=find([NaN;diff(WStbl.Epochs)]==-1 & WStbl.Epochs==1);
%     I=find([NaN;diff(WStbl.Epochs)]==-1);
%remove data if stage is not equal to ?
if 0
    temp = 1*(sum(WStbl.Epochs(I+[0:NbreathsInFuture])~=1,2)>0);
    I(temp==1)=[];
end



%remove breaths if different subjects
I((find(WStbl.Subjn(I)~=WStbl.Subjn(I+NbreathsInFuture))))=[];
I((find(WStbl.Subjn(I)~=WStbl.Subjn(I-NbreathsInPast))))=[];

%I(find(WStbl.Epochs(I)==3)==1)=[];
%I(find(WStbl.Epochs(I)==4)==1)=[];


%remove breaths if clearly beyond expected level
I((find(WStbl.Time(I+NbreathsInFuture)>WStbl.Time(I)+NbreathsInFuture*6)))=[];
I((find(WStbl.Time(I-NbreathsInPast)<WStbl.Time(I)-NbreathsInPast*6)))=[];


Nsamples = length(I)

medians2 = nanmean(VarX2_(I+[-NbreathsInPast:NbreathsInFuture]));
upper2 = medians2 + 1.96*nanstd(VarX2_(I+[-NbreathsInPast:NbreathsInFuture]))/length(I)^0.5;
lower2 = medians2 - 1.96*nanstd(VarX2_(I+[-NbreathsInPast:NbreathsInFuture]))/length(I)^0.5;


medians3 = nanmean(VarX3_(I+[-NbreathsInPast:NbreathsInFuture]));
upper3 = medians3 + 1.96*nanstd(VarX3_(I+[-NbreathsInPast:NbreathsInFuture]))/length(I)^0.5;
lower3 = medians3 - 1.96*nanstd(VarX3_(I+[-NbreathsInPast:NbreathsInFuture]))/length(I)^0.5;



medians = nanmean(VarX_(I+[-NbreathsInPast:NbreathsInFuture]));
upper = medians + 1.96*nanstd(VarX_(I+[-NbreathsInPast:NbreathsInFuture]))/length(I)^0.5;
lower = medians - 1.96*nanstd(VarX_(I+[-NbreathsInPast:NbreathsInFuture]))/length(I)^0.5;


mediansAR = mean(WStbl.EventsAr(I+[-NbreathsInPast:NbreathsInFuture]));
modehyp = mode(WStbl.Epochs(I+[-NbreathsInPast:NbreathsInFuture]));
%    medianSWSpred = mode(WStbl.SWSpred_logit(I+[-NbreathsInPast:NbreathsInFuture]));

beta_t = median(WStbl.Pbeta_ref(I+[-NbreathsInPast:NbreathsInFuture]));
alpha_t = median(WStbl.Palpha_ref(I+[-NbreathsInPast:NbreathsInFuture]));
theta_t = median(WStbl.Ptheta_ref(I+[-NbreathsInPast:NbreathsInFuture]));
delta_t = median(WStbl.Pdelta_ref(I+[-NbreathsInPast:NbreathsInFuture]));
%    sigma_t = median(WStbl.Psigma_ref(I+[-NbreathsInPast:NbreathsInFuture]));

%    medianPredID = nanmedian(WStbl.ARieF_predID(I+[-NbreathsInPast:NbreathsInFuture]));
%    medianPredID_MY = nanmedian(WStbl.ARieF_predID_MY(I+[-NbreathsInPast:NbreathsInFuture]));

Time_ = mean(WStbl.Time(I+[-NbreathsInPast:NbreathsInFuture]) - WStbl.Time(I));

%meanVI = nanmean(WStbl.VI(I+[-NbreathsInPast:NbreathsInFuture]));
Nsamples = length(I)
figure(17); clf(17); set(gcf,'color',[1 1 1]);
plot(1/60*Time_,medians,'-','color',[0.2 0.2 0.6]);
hold('on');
h=fill(1/60*[Time_ fliplr(Time_)],[upper, fliplr(lower)],[0.2 0.2 0.6],'Edgealpha',0,'facealpha',0.2);

plot(1/60*Time_,medians2,'-','color',[0.7 0.2 0.2]);
hold('on');
h=fill(1/60*[Time_ fliplr(Time_)],[upper2, fliplr(lower2)],[0.7 0.2 0.2],'Edgealpha',0,'facealpha',0.2);


plot(1/60*Time_,medians3,'-','color',[0.2 0.7 0.2]);
hold('on');
h=fill(1/60*[Time_ fliplr(Time_)],[upper3, fliplr(lower3)],[0.2 0.7 0.2],'Edgealpha',0,'facealpha',0.2);


% plot(Time_,upper,'k-','markersize',12);
% plot(Time_,lower,'k-','markersize',12);

%plot(Time_,mediansAR,'b-','markersize',12);
%plot(Time_,logit(medianPredID),'r-','markersize',12);
%plot(Time_,logit(medianPredID_MY),'r--','markersize',12);

%plot(Time_,-medianSWSpred,'g-','markersize',12);

plot(1/60*Time_,modehyp+2,'k-','markersize',12);
%plot(Time_,meanVI*2-3,'k-','markersize',12);

%     plot(Time_,beta_t*3-5,'b-','markersize',12);
%     plot(Time_,alpha_t*3-5,'b-','markersize',12);
%     plot(Time_,theta_t*3-5,'b-','markersize',12);
%     plot(Time_,delta_t*3-5,'b-','markersize',12);
%plot(Time_,sigma_t*3-5,'b:','markersize',12);

%plot(logitinverse(xval),num./den,'.-')

set(gca,'tickdir','out','box','off','fontname','arial narrow')
% hold('off')

yyaxis left
yticks = get(gca,'ytick');
ylims = get(gca,'ylim');
set(gca,'ytick',yticks,'ylim',ylims);
yyaxis right
set(gca,'ytick',yticks,'ylim',ylims,'yticklabels',round(logitinverse(yticks),3),'ycolor',[0 0 0],'tickdir','out');

xlim(1/60*[-100 NbreathsInFuture*dT+100])
ylim([-6 6])

%% Est SWS, within sleep
WStbl.SWS = (WStbl.Epochs==0)*1;
WStbl.SWS(isnan(WStbl.Epochs))=NaN;

WStbl.ExcludeForSWS = zeros(size(WStbl,1),1);

%add in these lines to remove wake, or remove arousals
%WStbl.ExcludeForSWS(WStbl.Epochs==4)=1;
WStbl.ExcludeForSWS(WStbl.Epochs==4&WStbl.EventsAr<0.5)=1;
WStbl.ExcludeForSWS(WStbl.Epochs~=4&WStbl.EventsAr>0.5)=1;
WStbl.ExcludeForSWS = logical(WStbl.ExcludeForSWS);

balance = nanmean(WStbl.SWS(~WStbl.ExcludeForSWS))
weightsSWS = 0*weights;
weightsSWS(WStbl.SWS>0.5)=1-balance;
weightsSWS(WStbl.SWS<=0.5)=balance;
sum(weightsSWS(WStbl.SWS>0.5&~WStbl.ExcludeForSWS))
sum(weightsSWS(WStbl.SWS<=0.5&~WStbl.ExcludeForSWS))

% mdl9 = compact(fitglm(WStbl,'SWS ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref*Palpha_ref + Pdelta_ref + Pdelta_ref^2','Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeForSWS,'weights',weightsSWS))
% mdl9.Rsquared.Ordinary

mdl9 = compact(fitglm(WStbl,'SWS ~ WakeIntensity + SleepIntensity','Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeForSWS,'weights',weightsSWS))
mdl9.Rsquared.Ordinary
try
mdl9b = compact(fitglm(WStbl,'SWS ~ WakeIntensity + L9w + L7w + L5w + L3w + L1w + N1w + N3w + N5w + N7w + N9w + SleepIntensity + L9s + L7s + L5s + L3s + L1s + N1s + N3s + N5s + N7s + N9s','Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeForSWS,'weights',weightsSWS))
mdl9b.Rsquared.Ordinary
end
cutoffSWS_ = -mdl9.Coefficients.Estimate(1)/mdl9.Coefficients.Estimate(2)
%
% mdl9 = compact(fitglm(WStbl,'SWS ~ ARieF_pred_logit + Pbeta_ref + Pdelta_ref','Distribution','binomial','Link','logit','Exclude',WStbl.ExcludeForSWS,'weights',weightsSWS))
% mdl9.Rsquared.Ordinary

WStbl.SWSpred_logit = logit(predict(mdl9,WStbl));

%% WeightsWSC and Linear regression version to develop SWS->W continuum "WSC"

%StageMedians_ = [ 4.3469   -2.6775   -3.1186   -4.2288   -5.7585]    W, N1, REM, N2, N3

WStbl.WSC = NaN*zeros(size(WStbl,1),1);
WStbl.WSC(WStbl.Epochs==4) = 4.3469;
WStbl.WSC(WStbl.Epochs==3) = -3.1186;
WStbl.WSC(WStbl.Epochs==2) = -2.6775;
WStbl.WSC(WStbl.Epochs==1) = -4.2288;
WStbl.WSC(WStbl.Epochs==0) = -5.7585;
WStbl.WSC(isnan(WStbl.Epochs))=NaN;

% ExcludeAR = zeros(size(WStbl,1),1);
% ExcludeAR(WStbl.Epochs==4&WStbl.EventsAr<0.5)=1;
% ExcludeAR(WStbl.Epochs~=4&WStbl.EventsAr>0.5)=1;
% ExcludeAR = logical(ExcludeAR);

NperStage(1) = sum(WStbl.Epochs==4&~WStbl.ExcludeAR);
NperStage(2) = sum(WStbl.Epochs==3&~WStbl.ExcludeAR)
NperStage(3) = sum(WStbl.Epochs==2&~WStbl.ExcludeAR)
NperStage(4) = sum(WStbl.Epochs==1&~WStbl.ExcludeAR)
NperStage(5) = sum(WStbl.Epochs==0&~WStbl.ExcludeAR)
Ntotal = sum(~isnan(WStbl.Epochs)&~WStbl.ExcludeAR)
sum(NperStage)

NperStageF = NperStage/sum(NperStage)

weightsWSC = zeros(size(WStbl,1),1);
weightsWSC(WStbl.Epochs==4) = 1/NperStageF(1);
weightsWSC(WStbl.Epochs==3) = 1/NperStageF(2);
weightsWSC(WStbl.Epochs==2) = 1/NperStageF(3);
weightsWSC(WStbl.Epochs==1) = 2/NperStageF(4);
weightsWSC(WStbl.Epochs==0) = 4/NperStageF(5);

sum(weightsWSC(WStbl.Epochs>=0&WStbl.Epochs<4))
sum(weightsWSC(WStbl.Epochs==4))

weightsWSC(WStbl.Epochs>=0&WStbl.Epochs<4) = weightsWSC(WStbl.Epochs>=0&WStbl.Epochs<4)/nansum(weightsWSC(WStbl.Epochs>=0&WStbl.Epochs<4));
weightsWSC(WStbl.Epochs==4) = weightsWSC(WStbl.Epochs==4)/nansum(weightsWSC(WStbl.Epochs==4));


sum(weightsWSC(WStbl.Epochs>=0&WStbl.Epochs<4))
sum(weightsWSC(WStbl.Epochs==4))

mdl10 = compact(fitglm(WStbl,'WSC ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2','Exclude',WStbl.ExcludeAR,'weights',weightsWSC))
mdl10.Rsquared.Ordinary
WStbl.WSCpred = predict(mdl10,WStbl);
thres = 0.5;
performanceWSC = PredictiveValue(1*(WStbl.WakeNoAR(~WStbl.ExcludeAR)>thres),1*(WStbl.WSCpred(~WStbl.ExcludeAR)>thres),WStbl.WakeNoAR(~WStbl.ExcludeAR))

%use as AR
WStbl.WSCpred_last1 = [NaN;WStbl.WSCpred(1:end-1)];
WStbl.WSCpred_last2 = [NaN;NaN;WStbl.WSCpred(1:end-2)];
WStbl.WSCpred_last3 = [NaN;NaN;NaN;WStbl.WSCpred(1:end-3)];
WStbl.WSCpred_last4 = [NaN;NaN;NaN;NaN;WStbl.WSCpred(1:end-4)];
WStbl.WSCpred_last5 = [NaN;NaN;NaN;NaN;NaN;WStbl.WSCpred(1:end-5)];
WStbl.WSCpred_last6 = [NaN;NaN;NaN;NaN;NaN;NaN;WStbl.WSCpred(1:end-6)];
WStbl.WSCpred_last6_min = min([WStbl.WSCpred_last1';WStbl.WSCpred_last2';WStbl.WSCpred_last3';WStbl.WSCpred_last4';WStbl.WSCpred_last5';WStbl.WSCpred_last6'])';

mdlAR = compact(fitglm(WStbl,'EventsAr ~ WSCpred + WSCpred_last6_min + WSCpred_last3 + WSCpred_last4','Exclude',ExcludeW,'Distribution','binomial','Link','logit','weights',weightsAR))
mdlAR.Rsquared.Ordinary

newARpred_ = predict(mdlAR,WStbl);
WStbl.newARpredF = newARpred_;
WStbl.newARpred_logit = logit(newARpred_);
WStbl.newARpred = 1*(newARpred_>0.5);
WStbl.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceAR = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),WStbl.EventsAr(~ExcludeW))
%performanceX = PredictiveValue(1*(WStbl.EventsAr(~ExcludeW)>0.5),1*(WStbl.ARieF_pred_logit(~ExcludeW)>-1.2954),WStbl.EventsAr(~ExcludeW))

figure(121); clf(121);
dStep=0.1;
Centers=-10:dStep:10;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
subplot(1,1,1);
[h7,edges] = histcounts(logit(newARpred_(~ExcludeW&WStbl.EventsAr>0.5)),Edges);
[h8,edges] = histcounts(logit(newARpred_(~ExcludeW&WStbl.EventsAr<0.5)),Edges);
bar(Centers,h7,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);





%%
if 0
    
    save mdlA mdlA -v7.3
    save mdlAR mdlAR -v7.3
    save RefTable RefTable -v7.3
    
end