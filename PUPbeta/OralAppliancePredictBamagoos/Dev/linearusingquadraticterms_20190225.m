%% Before running this
%1. Prepare excluded variables "Exclude"
%2. Prepare Amatrix with  {'Vpassive0p5','LGn','Vcomp','arthres0p5','VRA1'}

if 0
Yvariable_=DeltaAHIf;
%Yvariable_=DeltaAHIfnrem;
%Yvariable_=
endotypesubset=[1 2 3 4 5];

Amatrix_ = Amatrix(:,endotypesubset);
end
%xvalueslist = endotypelist(endotypesubset);
%Exclude = 0*Yvariable_;
forcespecificterms=[];
linregnotlogreg=1;   

criteriaPlot = [Yvariable>70&TreatmentAHI<10 Yvariable>70 Yvariable<=70&Yvariable>50 Yvariable<=50];
Coptions = [0 0.5 0;0.35 0.7 0;0.9 0.5 0;0.75 0.1 0.1];

% criteriaPlot = [Yvariable>70 Yvariable<=70];
% Coptions = [0 0.5 0;0.75 0.1 0.1];

plotrange1 = 0*Yvariable; 
plotrange2 = 0*Yvariable; 
    plotrange1(1:108)=1; 
YvariablePlot = Yvariable_;
    plotfinal=1;
    
    
%% multiple linear regr leave-1-out backwards elimination
Ilist = [1:length(endotypesubset)];
standardizeearlyforplots=0; %1
standardizelateforBetastd=0; %0
% A_=SVMModel.ScaleData.scaleFactor.*(A+SVMModel.ScaleData.shift);
% fZsvm = @(A) (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,SVMModel.ScaleData.scaleFactor.*(A+SVMModel.ScaleData.shift),SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias; 
% fZsvm_prescaled = @(A) (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,A,SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias; 
% Z = -fZsvm(A); %Z<0 is responder
% Z_ = -fZsvm_prescaled(A_);

Iincl=find(Exclude==0);

varlist_ = xvalueslist;
%Ilist=[1:4]; %change to LGn for residual event info ARI/CAI etc
varlist = varlist_(Ilist);
%includequadraticterms=1;
    %ignoresquaredterms=0;
    %maxP=0.15; %0.15-0.2
    
    
Amatrix2=Amatrix_;
%Amatrix2(Amatrix2(:,1)<0,1)=0; %option to put lower/upper limits on specific traits
shifts = [nanmean(Amatrix2(Iincl,Ilist))]
%shifts = [nanmedian(Amatrix2(Iincl,Ilist))]
%     shifts(1)=0.9;
%     shifts(2)=0.5;
%     shifts(3)=1.25;
%     shifts(4)=0;
    shift=-shifts;
if standardizeearlyforplots
Ain = (Amatrix2(:,Ilist)-shifts)./nanstd(Amatrix2(Iincl,Ilist));
else
Ain = (Amatrix2(:,Ilist)-nanmean(Amatrix2(Iincl,Ilist)));    
end
meansandstds = [nanmean(Amatrix2(Iincl,Ilist))' nanstd(Amatrix2(Iincl,Ilist))'];

%Ain = Amatrix(:,Ilist);
% 
% scaleFactors = 1./nanstd(Amatrix(Iincl,Ilist));
% shiftFactors = -nanmean(Amatrix(Iincl,Ilist),1);
% Ain=scaleFactor.*(A+shiftFactor);
Ain(Exclude==1,:)=NaN;

%%
labels = varlist;
if includequadraticterms
J=length(Ilist);
Q0=tril(ones(J,J),-ignoresquaredterms); %myf2 = @(A) K + [A]*L + sum(([A]*Q) .* [A], 2);
%A2 contains quadratic terms e.g. a^2 + b^2 + ab (for J=2)
A2=[];

for j=1:J
    for i=1:J
        if Q0(i,j)~=0
            A2=[A2,Ain(:,i).*Ain(:,j)];
            labels = [labels [varlist{j} '.' varlist{i}]];
        end
    end
end
Ain = [Ain A2]; 
end 

%%
clear PrPredict

Yvar = DeltaAHIf; %nrem
criteriaR = NaN*DeltaAHIp;
criteriaR = 1*(DeltaAHIp>50);
M=length(Yvar)
Yvar(isnan(criteriaR))=NaN;
IlistOrig = [1:size(Ain,2)];

% Irem = find(sum(IlistOrig==disallowlist')==1);
% IlistOrig(Irem)=[];


% Yvar = TreatmentAHI.^0.5;
% criteriaR = NaN*TreatmentAHI;
% criteriaR = 1*(TreatmentAHI>10);

R(Exclude==1)=NaN;
criteriaR(isnan(DeltaAHIp)==1|Exclude==1)=NaN;
sum(~isnan(criteriaR))
%criteriaR(81)=NaN; %obtained data only by forcing software to analyze noise
temp=nanstd(Ain(:,IlistOrig))';

%IlistOrig([6 8 11 15 20])=[]; %4+[0 5 10 15]
Ause = Ain(:,IlistOrig);

%Ause = [BaselineAHI BMI Age Sex Neck];
if standardizelateforBetastd
    %Ause = (Ause-nanmean(Ause))./(nanstd(Ause)); %never need this
    Ause = (Ause)./nanstd(Ause); %for fully standardized terms, but will destroy plot function
end

if 0
Yvar = (Yvar-nanmean(Yvar)./nanstd(Yvar));
end
weights = 1*ones(M,1);
if 0
    weights(criteriaR==0)=nansum(criteriaR==1);
    weights(criteriaR==1)=nansum(criteriaR==0);
    weights=weights/nanmean(weights);
elseif 0
    weights(Yvariable_>66.7)=2/nansum(Yvariable_>66.7);
    weights(Yvariable_>50&Yvariable_<=66.7)=1/nansum(Yvariable_>50&Yvariable_<=66.7);
    weights(Yvariable_>0&Yvariable_<=50)=1/nansum(Yvariable_>0&Yvariable_<=50);
    weights(Yvariable_<=0)=1/nansum(Yvariable_<=0);
    weights=weights/nanmean(weights);
elseif 0
    weights(Yvariable_>66.7)=2/nansum(Yvariable_>66.7);
    weights(Yvariable_>50&Yvariable_<=66.7)=1/nansum(Yvariable_>50&Yvariable_<=66.7);
    weights(Yvariable_>0&Yvariable_<=50)=1/nansum(Yvariable_>0&Yvariable_<=50);
    weights(Yvariable_<=0)=1/nansum(Yvariable_<=0);
    weights=weights/nanmean(weights);
else
    weights=0*weights+1;
end


VIF = diag(inv(corrcoef(Ain(Exclude==0,:))));

M=size(Ain,1);

    clear PrPredict_
    PredT_=NaN*BaselineAHI;
    for j=1:M+1
    Ilisttest = 1:size(Ause,2);
        %xvalueslist={'Vpassive','LGn','VRA1','Vcomp','arthres0p5','Vpassive.^2','LGn.^2','VRA1.^2','Vcomp.^2','arthres0p5.^2','Vpassive.^3','LGn.^3','VRA1.^3','Vcomp.^3','arthres0p5.^3','LG1','delay','Vpassive.^4','Vactive','BaselineAHI','Neck','BMI','Sex','Age'};%}; %'VpassiveT','VcompT','arthresT',
    pvals=Inf;
    Train=1:M;
    if j<M+1
        Train(j)=[];
    end
    
    for i=1:length(Ilisttest)
    %xvalueslist(Ilisttest);
        if linregnotlogreg
        %[Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),criteriaR(Train),'binomial'); %,'weights',weights
            [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),Yvar(Train),'normal','weights',weights(Train)); %,'weights',weights(Train)
        else
            [Btemp,~,temp] = glmfit(Ause(Train,Ilisttest),criteriaR(Train),'binomial','weights',weights(Train)); %,'weights',weights(Train)
        end
    
    pvals = temp.p(2:end);
    minterms=1;
    if ~isempty(forcespecificterms)
        I=find(sum(Ilisttest==forcespecificterms')==1);
        pvals(I)=-1+pvals(I); %information still available but coded as negative
        minterms=length(forcespecificterms);
    end
    BtempX=Btemp(2:end);
    [~,ii]=max(pvals);
    
    %i=i+5;.
    if max(pvals)>maxP&&length(Ilisttest)>minterms%max(pvals)>maxP&&length(Ilisttest)>4 %2
        %xvalueslist(Ilisttest(ii));
        Ilisttest(ii)=[];
    else
        break
    end
    end
    %xvalueslist(Ilisttest)';
    if linregnotlogreg
        PrPredict = Ause(:,Ilisttest)*Btemp(2:end)+Btemp(1);
    else
        PrPredict = Ause(:,Ilisttest)*Btemp(2:end)+Btemp(1);
        %PrPredict = 1./(1+exp(-PrPredict));
    end
    
    if 0 %re-run regression with same parameters but now use regularization
        warning('off')
       [Mdl] = fitrlinear(Ause(Train',Ilisttest),Yvar(Train),'Regularization','ridge','weights',weights(Train),'learner','leastsquares'); %ridge
       %[Mdl] = fitrlinear(Ause(:,Ilisttest),Yvar,'Regularization','ridge','solver','bfgs','weights',weights,'learner','svm'); %ridge
       %[Mdl] = fitrlinear(Ause(:,Ilisttest),Yvar,'Regularization','lasso','solver','svm'); %lasso
       PrPredict = predict(Mdl,Ause(:,Ilisttest)); 
    end
    %figure(10)
    %scatter(PrPredict(Train),DeltaAHIf(Train),20,[0.4 0.4 1],'filled','markerfacealpha',0.7); hold('on')
    %scatter(PrPredict(Test),DeltaAHIf(Test),20,[1 0.4 0.4],'filled','markerfacealpha',0.7); hold('off')
    [x,y,t,AUC,~] = perfcurve(criteriaR(Train)*1,PrPredict(Train),1); %need to find the threshold value that gives the OPTROCPT!
    
             [~,I]=max(y+(1-x));
             thresopt=mean(t(I:(I+1)));
             
             if j<M+1
                 PrPredict_(j,1)=PrPredict(j);
                PredT_(j) = (PrPredict(j)>thresopt)*1;
                if isnan(PrPredict(j))
                    PredT_(j) = NaN;
                end
             else
             PredT = (PrPredict>thresopt)*1;  
             PredT(isnan(PrPredict)) = NaN;
             %xvalueslist(Ilisttest)'
             end
    end
    
    performanceP_all = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),DeltaAHIp(rangevalidate));
    performanceT_all = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),TreatmentAHI(rangevalidate));
    performanceP = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),DeltaAHIp(rangevalidate));      
    performanceF = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),DeltaAHIf(rangevalidate));
    performanceF_all = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),DeltaAHIf(rangevalidate)); 
    converttop = @(x) 2*x/(x+1);
    [converttop(performanceF_all.predY_mean) converttop(performanceF_all.predY_mean - 1.96*performanceF_all.predY_sem) converttop(performanceF_all.predY_mean + 1.96*performanceF_all.predY_sem)]
    [converttop(performanceF_all.predN_mean) converttop(performanceF_all.predN_mean - 1.96*performanceF_all.predN_sem) converttop(performanceF_all.predN_mean + 1.96*performanceF_all.predN_sem)]
    
    [converttop(performanceF.predY_mean) converttop(performanceF.predY_mean - 1.96*performanceF.predY_sem) converttop(performanceF.predY_mean + 1.96*performanceF.predY_sem)]
    [converttop(performanceF.predN_mean) converttop(performanceF.predN_mean - 1.96*performanceF.predN_sem) converttop(performanceF.predN_mean + 1.96*performanceF.predN_sem)]
    
    
    performanceT = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentAHI(rangevalidate));
    performanceB = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHI(rangevalidate));
    performanceD = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHI(rangevalidate)-TreatmentAHI(rangevalidate));
    
    performanceTp33 = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentAHI(rangevalidate).^0.33);
    
    converttop2 = @(x) x.^3;
    [converttop2(performanceTp33.predY_mean) converttop2(performanceTp33.predY_mean - 1.96*performanceTp33.predY_sem) converttop2(performanceTp33.predY_mean + 1.96*performanceTp33.predY_sem)]
    [converttop2(performanceTp33.predN_mean) converttop2(performanceTp33.predN_mean - 1.96*performanceTp33.predN_sem) converttop2(performanceTp33.predN_mean + 1.96*performanceTp33.predN_sem)]
    
    performanceBp33 = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHI(rangevalidate).^0.33);
    [converttop2(performanceBp33.predY_mean) converttop2(performanceBp33.predY_mean - 1.96*performanceBp33.predY_sem) converttop2(performanceBp33.predY_mean + 1.96*performanceBp33.predY_sem)]
    [converttop2(performanceBp33.predN_mean) converttop2(performanceBp33.predN_mean - 1.96*performanceBp33.predN_sem) converttop2(performanceBp33.predN_mean + 1.96*performanceBp33.predN_sem)]
    
    performanceTp33all = PredictiveValue(1*criteriaR(rangevalidate),PredT(rangevalidate),TreatmentAHI(rangevalidate).^0.33);
    
    [converttop2(performanceTp33all.predY_mean) converttop2(performanceTp33all.predY_mean - 1.96*performanceTp33all.predY_sem) converttop2(performanceTp33all.predY_mean + 1.96*performanceTp33all.predY_sem)]
    [converttop2(performanceTp33all.predN_mean) converttop2(performanceTp33all.predN_mean - 1.96*performanceTp33all.predN_sem) converttop2(performanceTp33all.predN_mean + 1.96*performanceTp33all.predN_sem)]
    
    %
    
    %PrPredictB = PrPredict*2./(1+PrPredict);
    if linregnotlogreg
%     PrPredict(PrPredict>100)=100;
%     PrPredict_(PrPredict_>100)=100;
    PrPredict(PrPredict>1)=1;
    PrPredict_(PrPredict_>1)=1;
    end
    
    thresopt_ = thresopt*2./(1+thresopt)
    
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
    performanceP_all
    performanceP
    Afinal=Ause(Exclude==0,Ilisttest);
    
    T=table(Btemp(2:end),temp.p(2:end),labels(Ilisttest)')
    T_=table(Btemp(2:end),temp.p(2:end),labels(Ilisttest)')
    
    %% all data for table
    temp1=BaselineAHI(Iincl).^0.33;
    temp2=1/sqrt(93);
    [converttop2(mean(temp1)) converttop2(mean(temp1) - 1.96*std(temp1)*temp2) converttop2(mean(temp1) + 1.96*std(temp1)*temp2)]
    temp1=TreatmentAHI(Iincl).^0.33;
    [converttop2(mean(temp1)) converttop2(mean(temp1) - 1.96*std(temp1)*temp2) converttop2(mean(temp1) + 1.96*std(temp1)*temp2)]
    temp1=DeltaAHIf(Iincl);
    [converttop(mean(temp1)) converttop(mean(temp1) - 1.96*std(temp1)*temp2) converttop(mean(temp1) + 1.96*std(temp1)*temp2)]
    
    
    %%
    DeltaAIf = (BaselineAI-TreatmentAI)./(BaselineAI+TreatmentAI);
    performanceP_AI = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),DeltaAIf(rangevalidate));
    performanceT_AI = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentAI(rangevalidate).^0.33);
    
    performanceT_MinSpO2 = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentMinSpO2(rangevalidate));
    performanceB_MinSpO2 = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineMinSpO2(rangevalidate));
    
    temp = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHInrem(rangevalidate));
    temp = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentAHInrem(rangevalidate));
    temp = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),BaselineAHIrem(rangevalidate));
    temp = PredictiveValue(1*criteriaR(rangevalidate),PredT_(rangevalidate),TreatmentAHIrem(rangevalidate));
    
    DeltaAHIfnrem = (BaselineAHInrem-TreatmentAHInrem)./(BaselineAHInrem+TreatmentAHInrem);
    
    
    
    
    %%
    
    %scatter(PrPredictX,Yvar,20,[0.2 0.8 0.5],'filled','markerfacealpha',0.7);
%     set(gca,'Ytick',xticksset,'yticklabels',xticksset2)
%     if linregnotlogreg
%         set(gca,'Xtick',xticksset,'xticklabels',xticksset2)
%     end
%     subplot(1,2,2)
%     scatter(PrPredictX_,Yvar,20,[0.2 0.8 0.5],'filled','markerfacealpha',0.7);
%     if linregnotlogreg
%         set(gca,'Xtick',xticksset,'xticklabels',xticksset2)
%     end
%     set(gca,'Ytick',xticksset,'yticklabels',xticksset2)
%     
   
%     subplot(2,2,3)
%     warning('off')
%        [Mdl] = fitrlinear(Ause(:,Ilisttest),Yvar,'learner','leastsquares'); %ridge
%        %[Mdl] = fitrlinear(Ause(:,Ilisttest),Yvar,'Regularization','ridge','solver','bfgs','weights',weights,'learner','svm'); %ridge
%        %[Mdl] = fitrlinear(Ause(:,Ilisttest),Yvar,'Regularization','lasso','solver','svm'); %lasso
%        PrPredictX = predict(Mdl,Ause(:,Ilisttest)); 
%        PrPredictX(PrPredictX>1)=1;
%        scatter(PrPredictX,Yvar,20,[0.2 0.8 0.5],'filled','markerfacealpha',0.7);
%     set(gca,'Xtick',xticksset,'xticklabels',xticksset2)
%     set(gca,'Ytick',xticksset,'yticklabels',xticksset2)
       
    %checkPredict = Ause(:,Ilisttest)*Mdl.Beta+Mdl.Bias; %good
        %[checkPredict PrPredict]
    
    
    
    %%
    [r,p]=corr(Vpassive0p5(Iincl),DeltaAHIf(Iincl))
    [r,p]=corr(LGn(Iincl),DeltaAHIf(Iincl))
    [r,p]=corr(Vcomp(Iincl),DeltaAHIf(Iincl))
    [r,p]=corr(arthres0p5(Iincl),DeltaAHIf(Iincl))
    [r,p]=corr(VRA1(Iincl),DeltaAHIf(Iincl))
    %%
    [r,p]=corr(Vpassive0p5(Iincl),DeltaAHI(Iincl))
    [r,p]=corr(LGn(Iincl),DeltaAHI(Iincl))
    [r,p]=corr(Vcomp(Iincl),DeltaAHI(Iincl))
    [r,p]=corr(arthres0p5(Iincl),DeltaAHI(Iincl))
    [r,p]=corr(VRA1(Iincl),DeltaAHI(Iincl))
    
    
    
    