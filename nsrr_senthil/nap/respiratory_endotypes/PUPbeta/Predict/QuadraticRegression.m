% First make a version of QuadraticRegressionPrep for your study
% Run that.

%% Preliminary screening

PvalScreen = nan(1,size(xvalueslist,2));
if preliminaryscreenoutP>0
    for i=1:size(xvalueslist,2)
        [~,~,temp] = glmfit(Amatrix(:,i),criteriaR,'binomial','weights',weights); %,'weights',weights(Train)
        PvalScreen(i)=temp.p(2);
    end
    I = PvalScreen>preliminaryscreenoutP;
    Amatrix(:,I)=[];
    xvalueslist(I)=[];
end


%% multiple linear regr leave-1-out backwards elimination
Ilist = [1:size(Amatrix,2)];

Iincl=find(Exclude==0);

varlist_ = xvalueslist;
varlist = varlist_(Ilist);
shifts = [nanmean(Amatrix(Iincl,Ilist))];
shift=-shifts;
if standardizeearlyforplots
    Ain = (Amatrix(:,Ilist)-shifts)./nanstd(Amatrix(Iincl,Ilist));
else
    Ain = (Amatrix(:,Ilist)-nanmean(Amatrix(Iincl,Ilist)));
end
meansandstds = [nanmean(Amatrix(Iincl,Ilist))' nanstd(Amatrix(Iincl,Ilist))'];
Ain(Exclude==1,:)=NaN;

%%
labels = varlist;
if 1 || includequadraticterms %test that this just provides a set of zero-valued Q0
    J=length(Ilist);
    Q0=tril(ones(J,J),-ignoresquaredterms); %myf2 = @(A) K + [A]*L + sum(([A]*Q) .* [A], 2);
    %A2 contains quadratic terms e.g. a^2 + b^2 + ab (for J=2)
    
    Q0(linearonlyterms,:)=0;
    Q0(:,linearonlyterms)=0;
    
    Q0(squaredtermexcl,squaredtermexcl)=0;
    
    A2=nan(size(Ain,1),sum(sum(Q0)));
    %labels2 = cell();
    labels2 = repmat({''},sum(sum(Q0)),1)';
    count=1;
    for j=1:J
        for i=1:J
            if Q0(i,j)~=0
                A2(:,count)=[Ain(:,i).*Ain(:,j)];
                labels2{:,count} = [varlist{j} ':' varlist{i}];
                count=count+1;
            end
        end
    end
    labels = [labels labels2];
    Ain = [Ain A2];
end

%% Run the Leave1Out and Final Model Backwards Elimination Loop
clear PrPredict

M=length(Yvar)
Yvar(isnan(criteriaR))=NaN;
IlistOrig = [1:size(Ain,2)];

R(Exclude==1)=NaN;
criteriaR(isnan(Yvar)==1|Exclude==1)=NaN; % [changed from - criteriaR(isnan(DeltaAHIp)==1|Exclude==1)=NaN;]
sum(~isnan(criteriaR))
temp=nanstd(Ain(:,IlistOrig))';
Ause = Ain(:,IlistOrig);
if standardizelateforBetastd
    Ause = (Ause)./nanstd(Ause); %for fully standardized terms, but will destroy plot function
end
%VIF = diag(inv(corrcoef(Ain(Exclude==0,:))));


quadraticterm=nan(1,length(labels));
for iii=1:length(labels)
    quadraticterm(iii)=contains(labels(iii),':')==1 | contains(labels(iii),'^')==1;
end

LinearTerms = ~quadraticterm;

M=size(Ain,1);

%%
verbose = 1;
%%
clear PrPredict_
PredT_=NaN*Yvar;
for j=1:M+1 %patient leave1out loop
    if skipleaveoneout && j<M+1
        continue
    end
    disp(j)
    %xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2','Vpassive.^3','LGn.^3','VRA1.^3','Vcomp.^3','arthres0p5.^3','LG1','delay','Vpassive.^4','Vactive','BaselineAHI','Neck','BMI','Sex','Age'};%}; %'VpassiveT','VcompT','arthresT',
    pvals=Inf;
    Train=1:M;
    if j<M+1
        Train(j)=[];
    end
    
    %% Forwards - backwards stepwise regression
    
    stepdidnothing=[0 0];
    
    if BwdsOnly
        Ilisttest = 1:size(Ause,2);
        stepdidnothing=[1 0];  
    else
        Ilisttest = forcespecificterms; % initialize Ilisttest
    end
    
    if StartWithAllLinearTerms
        Ilisttest=Ilist;
    end
    
    for xx=1:Nloops
        
    if length(Ilisttest) < min([2 minterms]) ||  xx<minfwdsstartup || enableFullFwdsFullBwds
        Bkws = 0;
        Fwds = 1;
    else
        Bkws = 1;
        Fwds = 1;
    end
    
    if (stepdidnothing(1)==1)*enableFullFwdsFullBwds
        Bkws = 1;
        Fwds = 0;
    end
    
    if BwdsOnly
        Bkws = 1;
        Fwds = 0;
    end
    
    
    
    if Fwds
%         maxinModelP = 8;% number to be checked
%         if modifymaxP
%             if sum(AlreadyInModel)>maxinModelP 
%                 maxP = max(maxinModelP*maxP/(1.5*sum(AlreadyInModel)),0.05); % min 0.05 - can be changed later
%             else
%                 maxP = 0.1570;
%             end
%         else
%             maxP=0.1570;
%         end
            
        %Ilisttest = forcespecificterms;
       
         %for x=1:size(Ain,2)
            AlreadyInModel = zeros(size(Ain,2),1)'; AlreadyInModel(Ilisttest)=1;
            
            Nquadraticterms = sum(quadraticterm(Ilisttest));
            
            %AcceptableQuadraticTerms, ok to add based on current linear terms
            if Nquadraticterms == maxquadraticterms
                Q1 = zeros(size(Ain,2),1)';
            else
                Q1 = Q0;
                Iquadraticterms = Ilisttest(LinearTerms(Ilisttest)==1);
                for k=1:length(Iquadraticterms)
                    for jj=1:length(Iquadraticterms)
                        Q1(Iquadraticterms(k),Iquadraticterms(jj)) = Q1(Iquadraticterms(k),Iquadraticterms(jj))*2;
                    end
                end
                Q1=Q1(:);
                Q1(Q1==0)=[];
                Q1=Q1-1;
                Q1 = [zeros(size(Q0,1),1);Q1]';
            end
            
            if (sum(LinearTerms(Ilisttest))>=maxlinearterms)
                NewLinearTermsEnabled=0;
            else
                NewLinearTermsEnabled=1;
            end
            
            if (sum(quadraticterm(Ilisttest))>=maxquadraticterms)
                NewQuadraticTermsEnabled=0;
            else
                NewQuadraticTermsEnabled=1;
            end
            
            Candidates = AlreadyInModel==0 & (NewLinearTermsEnabled*(LinearTerms==1) | NewQuadraticTermsEnabled*(Q1==1));
            CandidatesList = find(Candidates==1);
            pvalsFWD = nan(length(CandidatesList),1);
            for i=1:length(CandidatesList)
                Ilisttest2 = [Ilisttest CandidatesList(i)];
                if linregnotlogreg
                    [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest2),Yvar(Train),'normal','weights',weights(Train)); %,'weights',weights(Train)
                else
                    [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest2),criteriaR(Train),'binomial','weights',weights(Train)); %,'weights',weights(Train)
                end
                pvalsFWD(i)=temp.p(end);
                
            end
            if PpenaltyQuadTerms~=1
                pvalsFWD(LinearTerms(Candidates)==0)=pvalsFWD(LinearTerms(Candidates)==0)*PpenaltyQuadTerms;
            end
            [minPval,ii]=min(pvalsFWD);
            
    
            if ~isempty(CandidatesList) && (minPval<maxP || length(Ilisttest)<minterms || sum(LinearTerms(Ilisttest))<minlinearterms) || sum((~LinearTerms(Ilisttest)))<minquadraticterms && (sum(LinearTerms(Ilisttest))<=maxlinearterms)
                Ilisttest = [Ilisttest CandidatesList(ii)]; %add terms
                if linregnotlogreg
                    [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),Yvar(Train),'normal','weights',weights(Train)); %,'weights',weights(Train)
                else
                    [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),criteriaR(Train),'binomial','weights',weights(Train)); %,'weights',weights(Train)
                end
                stepdidnothing(1)=0;
                if verbose==1
                    disp(['Loop=' num2str(xx) ', Fwds added: ' num2str(CandidatesList(ii)) ', p=' num2str(minPval)]);
                    disp(['  Features :' num2str(Ilisttest)]); 
                end
            else
                if linregnotlogreg
                    [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),Yvar(Train),'normal','weights',weights(Train)); %,'weights',weights(Train)
                else
                    [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),criteriaR(Train),'binomial','weights',weights(Train)); %,'weights',weights(Train)
                end
                stepdidnothing(1)=1;
                if verbose==1
                disp(['Loop=' num2str(xx) ', Fwds no features added']);
                disp(['Features :' num2str(Ilisttest)]);
                end
            end
    end
 %end
        
        
    
   
    
    if Bkws
       % maxP = 0.1570;
       %for i=1:length(Ilisttest) %backwards feature/trait elimination loop
       if stepdidnothing(1)==1 %didn't just run this above in Fwds
            if linregnotlogreg %needless if just did Fwds
                [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),Yvar(Train),'normal','weights',weights(Train)); %,'weights',weights(Train)
            else
                [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),criteriaR(Train),'binomial','weights',weights(Train)); %,'weights',weights(Train)
            end
       end
            
            pvalsBkws = temp.p(2:end);
            %minterms=1;
            
            currentforcedterms = forcespecificterms;
            if protectbaseterms
                labelset = labels(Ilisttest);
                ProtectedTerm=zeros(1,length(labelset));
                for i=1:length(labelset)
                    ProtectedTerm(i)=contains(labelset(i),'^')==0 & contains(labelset(i),':')==0 & sum(contains(labelset,labelset(i)))>1;
                end
                currentforcedterms = unique([currentforcedterms Ilisttest(ProtectedTerm==1)]);
            end
            
            if ~isempty(currentforcedterms)
                I=find(sum(Ilisttest==currentforcedterms')==1);
                pvalsBkws(I)=-1+pvalsBkws(I); %information still available but coded as negative
                %minterms=length(currentforcedterms);
            end
            
            if PpenaltyQuadTerms~=1
                pvalsBkws(LinearTerms(Ilisttest)==0)=pvalsBkws(LinearTerms(Ilisttest)==0)*PpenaltyQuadTerms;
            end
            
            BtempX=Btemp(2:end);
            
            if sum(LinearTerms(Ilisttest))<=minlinearterms
               %don't let the function remove linear terms
               pvalsBkws(LinearTerms(Ilisttest))=pvalsBkws(LinearTerms(Ilisttest))-1;
            end
            
            if sum((~LinearTerms(Ilisttest)))<=minquadraticterms
                pvalsBkws(~LinearTerms(Ilisttest))=pvalsBkws(~LinearTerms(Ilisttest))-1;
            end
            
            [maxPval,ii]=max(pvalsBkws);
            
            Nquadraticterms = sum(quadraticterm(Ilisttest));
            
            if (maxPval>maxP && length(Ilisttest)>minterms) || (Nquadraticterms > maxquadraticterms) || (sum(LinearTerms(Ilisttest))>maxlinearterms)
                if verbose==1, disp(['Loop=' num2str(xx) ', Bwds removed: ' num2str(Ilisttest(ii)) ', p=' num2str(maxPval)]); end
                Ilisttest(ii)=[];
                if verbose==1, disp(['  Features :' num2str(Ilisttest)]); end
                stepdidnothing(2)=0;
            else
                if verbose==1
                disp(['Loop=' num2str(xx) ', Bwds removed nothing']);
                disp(['  Features :' num2str(Ilisttest)]);
                end
                stepdidnothing(2)=1;
            end
            
        
    %end %backwards
    end
    
    if sum(stepdidnothing)==2
        break
    end
    
    end

    %% Make prediction model - cross-validation
    if linregnotlogreg
        PrPredict = Ause(:,Ilisttest)*Btemp(2:end)+Btemp(1);
    else
        PrPredict = Ause(:,Ilisttest)*Btemp(2:end)+Btemp(1);
    end
    
    [x,y,t,AUC,~] = perfcurve(criteriaR(Train)*1,PrPredict(Train),1); %need to find the threshold value that gives the OPTROCPT!
    
    [~,I]=max(y+(1-x));
    thresopt=mean(t(I:(I+1)));
    
    if j<M+1 %leave1out
        PrPredict_(j,1)=PrPredict(j);
        PredT_(j) = (PrPredict(j)>thresopt)*1;
        if isnan(PrPredict(j))
            PredT_(j) = NaN;
        end
    else %final model all patients
        PredT = (PrPredict>thresopt)*1;
        PredT(isnan(PrPredict)) = NaN;
    end
end

%% Performance
%variables coded with _all are based on all patients, final model
%continuous variable codes: P = percent reduction (untransformed) ...
%F is transformed percent reduction
%T is treatment value
%B is baseline value
%D is absolute delta
%underscore variables are cross-validated results i.e. PredT_ (vs PredT)
try
Tpredict = table(criteriaR,PredT,PredT_,T.DeltaAHIp,T.DeltaAHIpT,T.LGn,T.AHIbaseline,T.AHItreatment);
Tpredict.Properties.VariableNames = {'Responder','ModelResponder','CVModelResponder','DeltaAHIp','DeltaAHIpT','LGn','AHIbaseline','AHItreatment'};
catch 
Tpredict = table(criteriaR,PredT,PredT_);
Tpredict.Properties.VariableNames = {'Responder','ModelResponder','CVModelResponder'};
end

[performance_nonCV,raw_nonCV,Continuous_nonCV] = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),Yvar(rangevalidate))
% [~,~,ContBaselineDiff] = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),T.AHIbaseline(rangevalidate)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%backtransformed %reductions and 95%CI:
PredY_R_nonCV=[converttop(Continuous_nonCV.PredY(1)) converttop(Continuous_nonCV.PredY(1) - 1.96*Continuous_nonCV.PredY(2)) converttop(Continuous_nonCV.PredY(1) + 1.96*Continuous_nonCV.PredY(2))]
PredY_NR_nonCV=[converttop(Continuous_nonCV.PredN(1)) converttop(Continuous_nonCV.PredN(1) - 1.96*Continuous_nonCV.PredN(2)) converttop(Continuous_nonCV.PredN(1) + 1.96*Continuous_nonCV.PredN(2))]
[~,OR_nonCV_p,ORstats_nonCV] = fishertest([raw_nonCV([1 2]);raw_nonCV([4 3])])

%cross validated versions:
[performance_CV,raw_CV,Continuous_CV]  = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),Yvar(rangevalidate))

PredY_R_CV=[converttop(Continuous_CV.PredY(1)) converttop(Continuous_CV.PredY(1) - 1.96*Continuous_CV.PredY(2)) converttop(Continuous_CV.PredY(1) + 1.96*Continuous_CV.PredY(2))]
PredY_NR_CV=[converttop(Continuous_CV.PredN(1)) converttop(Continuous_CV.PredN(1) - 1.96*Continuous_CV.PredN(2)) converttop(Continuous_CV.PredN(1) + 1.96*Continuous_CV.PredN(2))]
[~,OR_CV_p,ORstats_CV] = fishertest([raw_CV([1 2]);raw_CV([4 3])])

%%

if 0
    %PrPredictB = PrPredict*2./(1+PrPredict);
    if ~linregnotlogreg
        %     PrPredict(PrPredict>100)=100;
        %     PrPredict_(PrPredict_>100)=100;
        PrPredict(PrPredict>1)=1;
        PrPredict_(PrPredict_>1)=1;
    end
    
    thresoptprob = 1/(1+exp(-thresopt));
    
    %figure(6)
    %subplot(1,2,1)
    xticksset=100*[-1 -0.5 0 1/3 0.7/1.3 0.9/1.1 1];
    xticksset2=[-Inf -200 0 50 70 90 100];
    if ~linregnotlogreg
        PrPredictX = 1./(1+exp(-PrPredict));
        PrPredictX_ = 1./(1+exp(-PrPredict_));
    else
        PrPredictX = PrPredict;
        PrPredictX_ = PrPredict_;
    end
    r=corr(PrPredictX(Iincl),Yvar(Iincl))
    
    Rsq = 1-sum((PrPredictX(Iincl)-Yvar(Iincl)).^2)/sum((Yvar(Iincl)-mean(Yvar(Iincl))).^2)
    
    VIF = diag(inv(corrcoef(Ause(Exclude==0,Ilisttest))));
    [pvals VIF];
    
    Afinal=Ause(Exclude==0,Ilisttest);
    
end
%% Final Model Table
T1=table(Btemp(2:end),temp.p(2:end),labels(Ilisttest)');
T1.Properties.VariableNames = {'Beta','p','Name'}

%% Nested model
if standardizefinalmodel
norm2SD = @(x) (x-nanmean(x))./(2*nanstd(x));
Tsd = T(:,[xvalueslist]); 
Tsd{:,:} = norm2SD(Tsd{:,:});
Tsd(:,OutcomeStr)=T(:,OutcomeStr);
else
    Tsd=T;
end
mdlFinal = fitglm(Tsd,[OutcomeStr '~' listtoeqn(labels(Ilisttest))],'Distribution',distr,'weights',weights)
mdlFinal.Rsquared.Ordinary;

if exist('Altlist','var')
    mdlFinalAlt = fitglm(T,[OutcomeStr '~' listtoeqn(Altlist)],'Distribution',distr,'weights',weights);
    mdlFinalAlt.Rsquared.Ordinary
else
    
end

mdlFinalT = mdlT(mdlFinal,1) %importance includes higher order terms
%mdlFinalT1 = mdlT(mdlFinal,0) %importance for base term based only on base term
mdlR2 = mdlFinal.Rsquared.Ordinary


%% Likelihood Ratio Test
% Only works for truly nested models
    % Get p value for likelihood ratio test that traits beat out AHI
    try
        if sum(sum(string(labels)==string(Altlist')))==length(Altlist)
            if mdlFinal.NumEstimatedCoefficients > mdlFinalAlt.NumEstimatedCoefficients
                mdl2=mdlFinal;
                mdl1=mdlFinalAlt;
            else
                mdl1=mdlFinal;
                mdl2=mdlFinalAlt;
            end
            [~,pmdlVersusAHI] = lratiotest(mdl2.LogLikelihood,mdl1.LogLikelihood,mdl2.NumEstimatedCoefficients-mdl1.NumEstimatedCoefficients) %Main H result
        end
    end
    
     %%
if 0    
   
    if 1
        M=40;
        clear out
        for i=1:50
            I = randsample(height(T),M,1);
            Tx = T(I,:);
            mdlFinalA = fitglm(Tx,[OutcomeStr '~' listtoeqn(labels(Ilisttest))],'Distribution',distr);
            %mdlFinalAT = mdlT(mdlFinalA);
            [predY,predYCI] = predict(mdlFinal,Tx);
            % predY=converttop(predY);
            % predYCI=converttop(predYCI);
            nanmedian(predY);
            nanmedian(diff(predYCI')'); % +/-30% (59.6817/2)
            
            [predY,predYCI] = predict(mdlFinalA,Tx);
            out(i)=nanmedian(diff(predYCI')');
        end
        nanmean(out)/2
    end
    
end

%% Random permutation
if 0
RandomPermOrder = randperm(height(T));
T_randPerm = T(RandomPermOrder,:);

PrPredict_RP = Ause(:,RandomPermOrder)*Btemp(2:end)+Btemp(1);

[x_RP,y_RP,t_RP,AUC_RP,~] = perfcurve(criteriaR(RandomPermOrder)*1,PrPredict(RandomPermOrder),1); %need to find the threshold value that gives the OPTROCPT!

[~,I_RP]=max(y_RP+(1-x_RP));
thresopt=mean(t_RP(I_RP:(I_RP+1)));


PredT_RP = (PrPredict_RP>thresopt_RP)*1;
PredT_RP(isnan(PrPredic_RPt)) = NaN;
end

%%
if ~any(Tpredict.Properties.VariableNames=="AHItreatment")
    return
end
%%
%%
figure(99); clf(99);
set(gcf,'color',[1 1 1]);
subplot(1,2,1)
I = find(Tpredict.ModelResponder==1 & Tpredict.CVModelResponder==0);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','-','color',[0.1 0.8 0.1].^0.33,'linewidth',1);
hold on
I = find(Tpredict.ModelResponder==1 & Tpredict.CVModelResponder==1);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','-','color',[0.1 0.8 0.1],'linewidth',1);

I = find(Tpredict.ModelResponder==1);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
Ymean = nanmean(Y.^0.5).^2;
plot(.3+0.2*[-1 1],Ymean(1)*[1 1],'-','color',[0 0 0],'linewidth',2);
plot(2-.3+0.2*[-1 1],Ymean(2)*[1 1],'-','color',[0 0 0],'linewidth',2);
hold off
set(gca,'xtick',[],'tickdir','out');
xlim([0 2]);
box off
ylim([0 100]);

I = find(Tpredict.ModelResponder==0);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
NRahis = nanmean(Y.^0.5).^2;

subplot(1,2,2)
I = find(Tpredict.ModelResponder==0 & Tpredict.CVModelResponder==1);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','-','color',[0.8 0.1 0.1].^0.33,'linewidth',1)
hold on
I = find(Tpredict.ModelResponder==0 & Tpredict.CVModelResponder==0);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','-','color',[0.8 0.1 0.1],'linewidth',1);

I = find(Tpredict.ModelResponder==0);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
Ymean = nanmean(Y.^0.5).^2;
plot(.3+0.2*[-1 1],Ymean(1)*[1 1],'-','color',[0 0 0],'linewidth',2);
plot(2-.3+0.2*[-1 1],Ymean(2)*[1 1],'-','color',[0 0 0],'linewidth',2);
hold off
set(gca,'xtick',[],'tickdir','out');
xlim([0 2]);
box off
ylim([0 100]);

%%
figure(101); clf(101);
set(gcf,'color',[1 1 1]);

I = find(Tpredict.Responder==0& Tpredict.CVModelResponder==0);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','-','color',[0.8 0.1 0.1],'linewidth',1);
hold on

I = find(Tpredict.Responder==0& Tpredict.CVModelResponder==1);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','--','color',[0.8 0.1 0.1],'linewidth',1);
hold on

% I = find(Tpredict.Responder==1);
% Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
% Ymean = nanmean(Y.^0.5).^2;
% plot(.3+0.2*[-1 1],Ymean(1)*[1 1],'-','color',[0 0 0],'linewidth',2);
% plot(2-.3+0.2*[-1 1],Ymean(2)*[1 1],'-','color',[0 0 0],'linewidth',2);


%I = find(Tpredict.Responder==0);
%Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
%nanmean(Y.^0.5).^2
I = find(Tpredict.Responder==1& Tpredict.CVModelResponder==1);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','-','color',[0.1 0.8 0.1],'linewidth',1);

hold on
I = find(Tpredict.ModelResponder==1 & Tpredict.CVModelResponder==0);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
X = repmat([0.3 1.7],length(I),1);
plot(X',Y','--','color',[0.1 0.8 0.1],'linewidth',1);


I = 1:height(Tpredict);
Y = [Tpredict.AHIbaseline(I) Tpredict.AHItreatment(I)];
Ymean = nanmean(Y.^0.5).^2;
plot(.3+0.2*[-1 1],Ymean(1)*[1 1],'-','color',[0 0 0],'linewidth',2);
plot(2-.3+0.2*[-1 1],Ymean(2)*[1 1],'-','color',[0 0 0],'linewidth',2);

set(gca,'xtick',[],'tickdir','out');
xlim([0 2]);
box off
ylim([0 100]);

%% Display this
ORstats_CV
 