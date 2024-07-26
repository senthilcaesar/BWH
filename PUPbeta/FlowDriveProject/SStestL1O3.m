clear all
close all
clc

%edited to keep test data at each step, slower

%load('SA_FD_16.mat');
load('SSblend_FS_FD_56.mat')
%%
clearvars -except PtData FeatureNames Amatrix


%% Just in case, look at NaN again

allnanrows = sum(isnan(Amatrix),2)==size(Amatrix,2);
sum(allnanrows) %still some rows with no data in here
Amatrix(allnanrows,:)=[];
PtData(allnanrows,:)=[];

Fnan=sum(isnan(Amatrix)|isinf(Amatrix))/size(Amatrix,1);
%checked=ok

%% Rem fluttering
if 1
% FtrsToExclude = {...
%     'InspFlutPow_Vpeak2_21_O',...
%     'ExpFlutPow_VpeakE2_22_O',...
%     'InspFlutPow_Vpeak2_20_T',...
%     'ExpFlutPow_VpeakE2_21_T',...
%     'InspExpFlutPow_22_T',...
%     'AA_PeaksRatio_25_T',...
%     'MedianFlowChange_26_T',...
%     'STDFlowChange_27_T',...
%     'AA_PeaksRatio_24_O',...
%     'MinFlowChange_25_O',...
%     'MedianFlowChange_26_O'...
%     }
FtrsToExclude = {...
    'AA_PeaksRatio',...
    }
Ind = [];
for i=1:length(FtrsToExclude) %needs checking
    temp=find(startsWith(FeatureNames.Name,FtrsToExclude(i)));
    if ~isempty(temp)
        Ind(end+1)=temp;
    end
end

FeatureNames(Ind,:)=[];
Amatrix(:,Ind)=[];
end


%% Names
temp = FeatureNames.Name;
clear temp1 temp2 temp3
for i=1:length(temp)
    temp1{i}=[temp{i} '_1p0'];
    temp2{i}=[temp{i} '_0p5'];
    temp3{i}=[temp{i} '_2p0'];
    tempX{i}=[temp{i} '_0pX'];
    tempY{i}=[temp{i} '_0pY'];
end

%% Make large matrix including square-root and square transformed data where possible

I1=Amatrix>=0;
I2=Amatrix<0;

AmatrixSQRT=NaN*Amatrix;
AmatrixSQRT(I1) = Amatrix(I1).^0.5;
AmatrixSQRT(I2) = -(-Amatrix(I2)).^0.5;

AmatrixSQ=NaN*Amatrix;
AmatrixSQ(I1) = Amatrix(I1).^2;
AmatrixSQ(I2) = -(-Amatrix(I2)).^2;


addextratransform=0;
AmatrixX=NaN*Amatrix;
AmatrixX = exp(-Amatrix);

AmatrixY=NaN*Amatrix;
AmatrixY = exp(Amatrix);

%AmatrixX(I2) = exp(-(-Amatrix(I2)));


% AmatrixSQ = Amatrix.^2; %removing if SQ/SQRT not appropriate etc
%     AmatrixSQ(Amatrix<0)=0;
%     AmatrixSQ(:,crit)=[];
%
%     temp2(crit)=[];
%     temp3(crit)=[];

%% Preparation
Gtest_ = PtData{:,13}; %Keep

squish=@(data)(100+((data(data>1)-1)*100).^0.5885)/100;
if 0
     Gtest_(Gtest_>1)=squish(Gtest_(Gtest_>1));
end
colofones = ones(length(Gtest_),1);
predyL1O = NaN*Gtest_;

NfeatureSteps=100;
predyL1O_array = NaN*ones(length(Gtest_),NfeatureSteps);

%% Weights 
maxG=1.5;
clear Ndata
dx=0.2; xbins=[0 0.3:dx:0.9 maxG];

%xbins=[0 0.1:0.05:1.1 1.5];
for i=1:length(xbins)-1
    Ix=Gtest_>xbins(i)&Gtest_<=xbins(i+1);
    Ndata(i)=sum(Ix);
end
%Ndata = Ndata.^0.5;
weightsbins = 1./(Ndata);
%weightsbins = [2 1 1 1 0.5];
weightsbins = weightsbins/mean(weightsbins);
%weightsbins = weightsbins/weightsbins(end);

weights = NaN*Gtest_;
for i=1:length(xbins)-1
    Ix=Gtest_>=xbins(i)&Gtest_<=xbins(i+1);
    weights(Ix)=weightsbins(i);
end
weights = weights/nanmean(weights);

useweights=1;
if ~useweights %overwrite, use no weights
    weights = ones(length(weights),1);
end

%% Other, off
if 0
Gtest = Gtest_;
Amatrix2 = [Amatrix AmatrixSQ AmatrixSQRT];
Labels = [temp1 temp2 temp3]';
nanx = sum(isnan(Amatrix2),2)>0;
    nany = isnan(Gtest);
    Amatrix2(nanx|nany,:)=[];
    Gtest(nanx|nany)=[];
    
colofones = 1 + 0*Gtest;
clear Rsq

for i=1:size(Amatrix2,2)
    [b,bint,r,rint,stats] = regress(Gtest(:),[colofones(:) Amatrix2(:,i)]);
    Rsq(i)=stats(1);
end

%large correlation array
clear Rsq
for i=1:size(Amatrix,2)
    for j=1:i-1
        
        %[b,bint,r,rint,stats] = regress(Amatrix2(:,j),[colofones(:) Amatrix2(:,i)]);
        %Rsq(i,j)=stats(1);
        I = ~isnan(Amatrix(:,j))&~isnan(Amatrix(:,i));
       Rsq(i,j)= corr(Amatrix(I,j),Amatrix(I,i)).^2;
    end
end

Inanrsq = find(isnan(Rsq))
Labels(Inanrsq)

temp = flipud(sortrows([Rsq;1:length(Rsq)]'));

FirstN=50;
ItopN = temp(1:FirstN,2);

if 0
figure(1); clf(1);
for i=1:FirstN
    subplot(6,10,i);
    scatter(Amatrix2(I,ItopN(i)),Gtest(I),2,'filled','markerfacealpha',0.1)
end
end
Amatrix2 = Amatrix2(:,ItopN);

Labels = Labels(ItopN);
end

%% Leave one out loop
clear RvalueTrain Err ErrRms Nfeatures
clear RsqTrain_array ErrTrain_array ErrRmsTrain_array              
labels_Step_Subj =[];               
alpha=0.05;
MaxNfeatures=1;
neverbreak=1;
   tic
    
Fwd=0; %backwards is crashing, perhaps duplicate data included
for subj=1:max(unique(PtData.PT))
    disp(['Starting Subj=' num2str(subj)]);
    if 1
        Amatrix2 = [Amatrix AmatrixSQ AmatrixSQRT];
    else
        Amatrix2 = [Amatrix];
    end
    Labels = [temp1 temp2 temp3]';
    if addextratransform>0
        Amatrix2 = [Amatrix2 AmatrixX];
        Labels = [Labels;tempX'];
    end
    if addextratransform>1
        Amatrix2 = [Amatrix2 AmatrixY];
        Labels = [Labels;tempY'];
    end
    Labels_ = Labels;
    Gtest = PtData{:,13}; %Edit per loop
    Gtest(Gtest>maxG)=maxG;
    Isubj=(PtData.PT==subj);
    
    Gtest(Isubj)=[];
    
    w = weights;
    w(Isubj)=[];
    w=w/nanmean(w);
    
    colofones_ = colofones;
    colofones_(Isubj)=[];
    
    %isnanF=sum(isnan(Amatrix2))/size(Amatrix2,1);
    %find(isnanF>0.001)
    
    Amatrix2_ = Amatrix2;
    Amatrix2_(Isubj,:)=[];
    
    %remove all NAN from all rows
    nanx = sum(isnan(Amatrix2_),2)>0;
    nany = isnan(Gtest);
    Amatrix2_(nanx|nany,:)=[];
    Gtest(nanx|nany)=[];
    w(nanx|nany)=[];
    w=w/nanmean(w);
    colofones_(nanx|nany)=[];
    warning('off')
    maxp=Inf;
    if Fwd==0 %backwards elimination
        If = 1:size(Amatrix2,2);
        while length(If)>=MaxNfeatures || maxp>alpha
            if 1
                [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2_(:,If),Gtest,w,0); %faster
                Pvals(1)=[];
            elseif 0
                [b,dev,stats]=glmfit(Amatrix2_(:,If),Gtest,'normal','weights',w);
                Pvals = stats.p(2:end);
            end
            if length(If)<=NfeatureSteps
                predyL1O_array(Isubj,length(If)) = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
                
                predytrain = [colofones_ Amatrix2_(:,If)]*b;
                predytrain(predytrain>maxG)=maxG;
                predytrain(predytrain<0)=0; %will be overwritten
                RsqTrain_array(subj,length(If)) = 1-nansum(w.*(Gtest-predytrain).^2)/nansum(w.*(Gtest-nanmean(Gtest)).^2);
                ErrTrain_array(subj,length(If)) = nanmean(w.*abs(predytrain-Gtest));
                ErrRmsTrain_array(subj,length(If)) = nanmean((w.*(predytrain-Gtest)).^2).^0.5;
                labels_Step_Subj{subj,length(If)} = If;
            end
            %remove least important
            if 0&length(If)>200 %speed up at higher levels
                temp=[Pvals';1:length(Pvals)]';
                temp=sortrows(temp);
                remi = temp(end-20+1:end,2);
                maxp = temp(end-20+1,1);
            else
                [maxp,remi]=max(Pvals);
            end
            if length(If)==1 || (~neverbreak && maxp<alpha && (length(If)-1)<=MaxNfeatures)
                break
            end
            disp(['Removing:p=' num2str(maxp)])
            If(remi)=[];
            Labels(remi)=[];
        end
    else %forwards stepwise inclusion
        if 1 %Matlab stock
            [b,se,Pvals,inmodel,stats,nextstep,history] = stepwisefit([Amatrix2_],Gtest,'penter',alpha,'premove',2*alpha,'MaxIter',MaxNfeatures); %no weighting possible here
            If = [];
            for i=1:sum(inmodel)
                temp = find(history.in(i,:));
                x=sum(temp'==If,2);
                temp(x==1)=[];
                if ~isempty(temp)
                    If(end+1)=temp;
                end
            end
        else %robust
            If = [];
            notinc = 1:size(Amatrix2_,2);
            while length(If)<=MaxNfeatures %&& maxp<alpha
                RMSE=NaN*notinc;
                P=NaN*notinc;
                tic
                for i=1:length(notinc)
                    %[~,Pvals,RMSE(i),~]=glmfitFast(Amatrix2_(:,[If notinc(i)]),Gtest,w);
                    %P(i) = Pvals(end);
                    RMSE(i)=glmfitFastRMSE(Amatrix2_(:,[If notinc(i)]),Gtest,w);
                end
                toc
                [~,mini]=min(RMSE);
                [~,Pvals,~,~]=glmfitFast(Amatrix2_(:,[If notinc(mini)]),Gtest,w,0);
                if Pvals(end)>alpha
                    disp('end FSR');
                    break
                end
                If = [If notinc(mini)];
                notinc(mini)=[];
                Labels_(If)
                %check for remove
                [~,Pvals,RMSE,~]=glmfitFast(Amatrix2_(:,If),Gtest,w,0);
                Pvals(1)=[];
                while max(Pvals)>2*alpha
                [~,remi]=max(Pvals);
                notinc = [notinc If(remi)];
                If(remi)=[];
                [~,Pvals,RMSE,~]=glmfitFast(Amatrix2_(:,If),Gtest,w,0);
                Pvals(1)=[];
                end
            end
        end
        [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2_(:,If),Gtest,w,0);
        Pvals(1)=[];
        Labels = Labels_(If);
    end
    
    Bvals = b(2:end);
    Bvals_ = abs(Bvals./std(Amatrix2_(:,If))');
    
    temp = std(Amatrix2_(:,If));
    if 1
    Data=[Pvals,(1:length(Pvals))'];
    Data2 = sortrows(Data);
    If2 = Data2(:,2);
    else
    Data=[Bvals_,(1:length(Bvals_))'];
    Data2 = sortrows(Data,'descend');
    If2 = Data2(:,2);
    end
    Labels2 = Labels(If2);
    P = Pvals(If2);
    B_ = Bvals_(If2);
    
    predytrain = [colofones_ Amatrix2_(:,If)]*b;
    predytrain(predytrain>maxG)=maxG;
    predytrain(predytrain<0)=0;
    Rsq = 1-nansum(w.*(Gtest-predytrain).^2)/nansum(w.*(Gtest-nanmean(Gtest)).^2);
    RvalueTrain(subj) = Rsq^0.5;
    Err(subj) = nanmean(abs(predyL1O(Isubj)-Gtest_(Isubj)));
    ErrRms(subj) = nanmean((predyL1O(Isubj)-Gtest_(Isubj)).^2).^0.5;
    
    predyL1O(Isubj) = [colofones(Isubj) Amatrix2(Isubj,If)]*b;

    Nfeatures(subj) = length(Labels2);
    %keyboard
end
toc
 
%% Most important features (penalty scores = square of distance from first position)
% DM comments:
% this method rewards those that appear early, but... if one feature usually appears
% early (and is accumulating an overall fairly low (i.e good) score), and then shows
% up last in just a couple of patients, it gets badly damaged by this scoring.
% An alternative would be to score as (x^-1), and pick the features with the highest 
% score. This way, features that regularly show up early, get the biggest scores. 
% The features that never show up early, never get big scores. The scenario described 
% above (a feature that is usually early, but on a couple occasions is late) is ok. 

score=zeros(1,size(Amatrix2,2));
maxscore = (size(labels_Step_Subj,2))^2;
for i=1:size(labels_Step_Subj,1)
    score1 = maxscore + 0*score;
    If = [];
            for j=1:size(labels_Step_Subj,2)
                temp = labels_Step_Subj{i,j};
                x=sum(temp'==If,2);
                temp(x==1)=[];
                if ~isempty(temp)
                    If(end+1)=temp;
                end
            end
    for j=1:length(If)        
        score1(If(j)) = (j-1).^0.5;
    end
    score=score+score1;
end

scoredata = [score;1:length(score)]';
    scoredata = sortrows(scoredata);
    LabelsOrdered = Labels_(scoredata(:,2));
    
    
    %% Process test data
    
    for i=1:size(predyL1O_array,2)
        RsqL1O_(i) = 1-nansum(weights.*(Gtest_-predyL1O_array(:,i)).^2)/nansum(weights.*(Gtest_-nanmean(Gtest_)).^2);
        ErrL1O_(i) = nanmean(weights.*abs(predyL1O_array(:,i)-Gtest_));
        ErrL1Orms_(i) = nanmean((weights.*(predyL1O_array(:,i)-Gtest_)).^2).^0.5;
    end
    figure(21)
    plot([RsqL1O_;ErrL1O_;ErrL1Orms_]')
    
   
    %% Manually pick optimal number of features
    
    NfeaturesOpt=28; 
    LabelsOrderedOpt = LabelsOrdered(1:NfeaturesOpt);
    
     %% Final model using all data
    If = 1:size(Amatrix2,2);
    Labels = Labels_;
    MaxNfeatures=NfeaturesOpt;
        while length(If)>=MaxNfeatures
               [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2(:,If),Gtest_,weights,1); %faster
                Pvals(1)=[];
            [maxp,maxpi]=max(Pvals);
            if length(If)==MaxNfeatures
                break
            end
            disp(['Removing:p=' num2str(maxp)])
            If(maxpi)=[];
            Labels(maxpi)=[];
        end
        Labels = Labels_(If)
  a=sortrows( [Pvals (1:NfeaturesOpt)'])
     Labels(a(:,2))
     direction = b(2:end)>0;
     
     
    
%% What data to use
if 0 %training, last loop data
    %Gtest = Gtest;
    %predy = predyL1O;
    predy = predytrain;% predyL1O_fsmr;
    %Yval = Gtest(NotI);
    Yval = Gtest;
    Rvalue = RvalueTrain;
    Err_ = nanmean(abs(predy-Gtest));
    STD_ = nanstd(predy-Gtest);
else %test data
    if 0
    RsqL1O = 1-nansum(weights.*(Gtest_-predyL1O).^2)/nansum(weights.*(Gtest_-nanmean(Gtest_)).^2);
    Rvalue = RsqL1O^0.5;
    ErrL1O = nanmean(abs(predyL1O-Gtest_));
    ErrL1Orms = nanmean((predyL1O-Gtest_).^2).^0.5;
    
    predy = predyL1O;% predyL1O_fsmr;
    Yval = Gtest_;
    else
       
       Rvalue = RsqL1O_(NfeaturesOpt)^0.5;
       predy = predyL1O_array(:,NfeaturesOpt);% predyL1O_fsmr;
       ErrL1O = ErrL1O_(NfeaturesOpt);
       ErrL1Orms = ErrL1Orms_(NfeaturesOpt);
       Yval = Gtest_;
       
    end  
end




%% Custom colormap
customcmap = ...
    [0.0416666666666667,0,0;0.0833333333333333,0,0;0.125000000000000,0,0;0.166666666666667,0,0;0.208333333333333,0,0;0.250000000000000,0,0;0.291666666666667,0,0;0.333333333333333,0,0;0.375000000000000,0,0;0.416666666666667,0,0;0.458333333333333,0,0;0.500000000000000,0,0;0.541666666666667,0,0;0.583333333333333,0,0;0.625000000000000,0,0;0.666666666666667,0,0;0.708333333333333,0,0;0.750000000000000,0,0;0.791666666666667,0,0;0.833333333333333,0,0;0.875000000000000,0,0;0.916666666666667,0,0;0.958333333333333,0,0;1,0,0;1,0.0416666666666667,0;1,0.0833333333333333,0;1,0.125000000000000,0;1,0.166666666666667,0;1,0.208333333333333,0;1,0.250000000000000,0;1,0.291666666666667,0;1,0.333333333333333,0;1,0.375000000000000,0;1,0.416666666666667,0;1,0.458333333333333,0;1,0.500000000000000,0;1,0.541666666666667,0;1,0.583333333333333,0;1,0.625000000000000,0;1,0.666666666666667,0;1,0.708333333333333,0;1,0.750000000000000,0;1,0.791666666666667,0;1,0.833333333333333,0;1,0.875000000000000,0;1,0.916666666666667,0;1,0.958333333333333,0;1,1,0;1,1,0.0625000000000000;1,1,0.125000000000000;1,1,0.187500000000000;1,1,0.250000000000000;1,1,0.312500000000000;1,1,0.375000000000000;1,1,0.437500000000000;1,1,0.500000000000000;1,1,0.562500000000000;1,1,0.625000000000000;1,1,0.687500000000000;1,1,0.750000000000000;1,1,0.812500000000000;1,1,0.875000000000000;1,1,0.937500000000000;1,1,1];

%% Analysis and plots

figure(2); clf(2);
subplot(1,3,1);
scatter(100*predy,100*Yval,2,'filled','markerfacealpha',0.4)
xlim([0 150]);
hold('on')

%dx=0.1; xbins=[0 0.15:dx:1.05 1.5];
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
title(['Rvalue: ' num2str(Rvalue(end),2)])
%plot binned data
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy>xbins(i)&predy<xbins(i+1);
    medianX(i)=prctile(predy(Ix),50);
    medianY(i)=prctile(Yval(Ix),50);
    upperIQRY(i)=prctile(Yval(Ix),75);
    lowerIQRY(i)=prctile(Yval(Ix),25);
    upperIQRY2(i)=prctile(Yval(Ix),90);
    lowerIQRY2(i)=prctile(Yval(Ix),10);
end
%medianX = [0.2:dx:1]; % overwrite
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01)
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01)

%% moving time median
if 0
xdata = [0 medianX];
ydata = [0 medianY];
%plot(100*xdata,100*ydata,'r');

x0 = [1 1];
upper = [1.1 2];
lower = [0.9 0.5];
fun = @(x,xdata)x(1)*(xdata.^x(2));    
[Xmodel]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper); 
xin = 0:0.01:1.5;
yout = fun(Xmodel,xin);
plot(100*xin,100*yout,'r');

predybackup = predy; 
if 0
    predy = predybackup;
end
% replot (overwrite)

predy = fun(Xmodel,predy);

figure(2); clf(2);
subplot(1,3,1);
scatter(100*predy,100*Yval,2,'filled','markerfacealpha',0.4)
xlim([0 150]);
hold('on')

%dx=0.1; xbins=[0 0.15:dx:1.05 1.5];
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
title(['Rvalue: ' num2str(Rvalue(end),2)])
%plot binned data
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy>xbins(i)&predy<xbins(i+1);
    medianX(i)=prctile(predy(Ix),50);
    medianY(i)=prctile(Yval(Ix),50);
    upperIQRY(i)=prctile(Yval(Ix),75);
    lowerIQRY(i)=prctile(Yval(Ix),25);
    upperIQRY2(i)=prctile(Yval(Ix),90);
    lowerIQRY2(i)=prctile(Yval(Ix),10);
end
%medianX = [0.2:dx:1]; % overwrite
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01)
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01)
end

%% create classes from regression SVM results
subplot(1,3,2)

%classcutoffs = [0.7 0.4];
classcutoffs = [0.9 0.7 0.5 0.3]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
Nclasses=length(classcutoffs)+1;

g = NaN*Yval;
ghat = NaN*predy;
g(Yval>classcutoffs(1))=1;
ghat(predy>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval<classcutoffs(i-1))=i;
    ghat(predy<classcutoffs(i-1))=i;
end

[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100 %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100

% plotconfusion(g',ghat') %needs neuralnet tools, and
%                           data must be in columns not rows, and
%                           works for binary classification only
AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));

AccAactual = mean(sum(AccA_C_Factual,2));
AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));

AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));

AccAestimated = mean(sum(AccA_C_Festimated));
AccBestimated = mean(sum(AccB_C_Festimated)) + mean(sum(AccB_C_Festimated));

ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
%second col gives accuracy if accepting next-category error as correct

x = order';                  % Random data
y = order';
if 0
    C1 = C_Festimated;
else
    C1 = C_Factual;
end
%C = [0 1 2 3 ; 1 2 3 4 ; 2 3 4 5];
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
set(gcf,'colormap',customcmap);
pcolor(XGrid,YGrid,C2)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
%colormap(gcf,'hot');
%cmap = get(gcf,'colormap');
h=colorbar();
set(gcf,'colormap',customcmap);

set(h,'FontName','Arial Narrow','Limits',[0 max(max(C2))]);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<18%((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
xlabel('Flow Shape Classification');
ylabel('Actual');
labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));

    
%% create classes from regression SVM results
if 1
    subplot(1,3,3)
    classcutoffs = [0.9 0.5];
    labeltext={'Normal','Mild-Moderate','Severe'};
    %classcutoffs = [0.9 0.7];
    %classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
    Nclasses=length(classcutoffs)+1;
    
    g = NaN*Yval;
    ghat = NaN*predy;
    g(Yval>classcutoffs(1))=1;
    ghat(predy>classcutoffs(1))=1;
    for i=2:Nclasses
        g(Yval<classcutoffs(i-1))=i;
        ghat(predy<classcutoffs(i-1))=i;
    end
    
    [C,order] = confusionmat(g,ghat);
    %sumactual=sum(C')';
    sumactual=sum(C,2);
    sumestimated=sum(C);
    %C_Factual=C./sum(C')'
    C_Factual=C./sum(C,2)*100 %rows are actual, cols are estimated
    C_Festimated=C./sum(C)*100
    
    % plotconfusion(g',ghat') %needs neuralnet tools, and
    %                           data must be in columns not rows, and
    %                           works for binary classification only
    %         AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    %         AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    %
    %         AccAactual = mean(sum(AccA_C_Factual,2));
    %         AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    %
    %         AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    %         AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    %
    %         AccAestimated = mean(sum(AccA_C_Festimated));
    %         AccBestimated = mean(sum(AccB_C_Festimated)) + mean(sum(AccB_C_Festimated));
    %
    %         ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    
    
    x = order';                  % Random data
    y = order';
    if 1
        C1 = C_Festimated;
    else
        C1 = C_Factual;
    end
    %C = [0 1 2 3 ; 1 2 3 4 ; 2 3 4 5];
    xSplit = diff(x)/2;                 % Find edge points
    ySplit = diff(y)/2;
    xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
    yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
    [XGrid, YGrid] = meshgrid(xEdges,yEdges);
    YGrid = flipud(YGrid);              % To match expected behavior
    XGrid = fliplr(XGrid);
    C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
    pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
    hold on                             % Plot original data points
    [X,Y] = meshgrid(x,y);
    colormap(gcf,'hot')
    colorbar();
    C1 = flipud(C1);
    C1 = fliplr(C1);
    C1 = C1';
    for i=1:size(C1,1)
        for j=1:size(C1,2)
            if C1(i,j)<((max(max(C1))-min(min(C1)))/2+min(min(C1)))
                textcolor=[1 1 1];
            else
                textcolor=[0 0 0];
            end
            text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
        end
    end
    xlabel('Flow Shape Classification');
    ylabel('Actual');
    
    yticks(y); yticklabels(gca,fliplr(labeltext));
    xticks(x); xticklabels(gca,fliplr(labeltext));
    
end


%%
AAsummary = [alpha MaxNfeatures median(Nfeatures) Fwd useweights median(RvalueTrain) Rvalue ACCs(2,:) median(Err) median(ErrRms) ErrL1O ErrL1Orms];
%%
return
%% Patient level summary data

%Need to put back breaths that are thrown out because of low VE here (apneas)
for i=1:max(PtData.PT)   
    I=PtData.PT==i&PtData.NotAr==1;
    medianG(i) = nanmedian(Gtest_(I));
    medianGest(i) = nanmedian(predy(I));
    
    denominator = sum(~isnan(Gtest_(I)));
    
    Fmildmod(i) = sum(Gtest_(I)<0.9&Gtest_(I)>0.5)/denominator;
    Fmildmodest(i) = sum(predy(I)<0.9&Gtest_(I)>0.5)/denominator;
    
    Fmodsev(i) = sum(Gtest_(I)<0.7)/denominator;
    Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
    
    Fsev(i) = sum(Gtest_(I)<0.5)/denominator;
    Fsevest(i) = sum(predy(I)<0.5)/denominator;
end
corr(medianGest',medianG')
corr(Fmodsevest',Fmodsev')
corr(Fsevest',Fsev')
corr(Fmildmodest',Fmildmod')
figure(9); clf(9);
subplot(1,5,1);
[Rtemp,Ptemp]=plotregressionwithSEM(medianGest,medianG); title(num2str(Rtemp));
subplot(1,5,2);
[Rtemp,Ptemp]=plotregressionwithSEM(Fmodsevest,Fmodsev); title(num2str(Rtemp));
subplot(1,5,3);
[Rtemp,Ptemp]=plotregressionwithSEM(Fsevest,Fsev); title(num2str(Rtemp));
subplot(1,5,4);
[Rtemp,Ptemp]=plotregressionwithSEM(Fmildmodest,Fmildmod); title(num2str(Rtemp));

%% Correlate with AHI
figure(10); clf(10);
subplot(1,1,1);
[Rtemp,Ptemp]=plotregressionwithSEM(AHItotal',medianGest); title(num2str(Rtemp));

[b,dev,stats]=glmfit([AHItotal medianGest'],medianG')
[b,dev,stats]=glmfit([AHItotal Fsevest'],Fsev')


%% Repeat, sensitivity for Nfeatures


for n=1:100
    clear medianG medianGest Fmildmod Fmildmodest Fmodsev Fmodsevest Fsev Fsevest
    predy = predyL1O_array(:,n);
for i=1:max(PtData.PT)   
    I=PtData.PT==i&PtData.NotAr==1;
    medianG(i) = nanmedian(Gtest_(I));
    medianGest(i) = nanmedian(predy(I));
    
    denominator = sum(~isnan(Gtest_(I)));
    
    Fmildmod(i) = sum(Gtest_(I)<0.9&Gtest_(I)>0.5)/denominator;
    Fmildmodest(i) = sum(predy(I)<0.9&Gtest_(I)>0.5)/denominator;
    
    Fmodsev(i) = sum(Gtest_(I)<0.7)/denominator;
    Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
    
    Fsev(i) = sum(Gtest_(I)<0.5)/denominator;
    Fsevest(i) = sum(predy(I)<0.5)/denominator;
end
n
Rvals(n,:)=[...
corr(medianGest',medianG'),...
corr(Fmodsevest',Fmodsev'),...
corr(Fsevest',Fsev'),...
corr(Fmildmodest',Fmildmod'),...
];
end



%%
predy(predy<0)=0;

%Plot FL values over time for an example subject (e.g. compare with Spike data). Random check (n=8) Looks great.
n=8; % 1313?
%for n=1:30
%ypredAll = predict(SVMModel_reg,Amatrix2(:,If));
%ypred = predict(SVMModel_reg,Amatrix2(:,If));
Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) predy(PtData.PT==n) Gtest_(PtData.PT==n)];

addNaNgaps=1;
cols = [NaN 1 2];
if addNaNgaps
    tol2=0.1;
    i=1;
    M=size(Data1,2);
    while i<(size(Data1,1)-1)
        if (Data1(i,cols(2))+Data1(i,cols(3))+tol2)<Data1(i+1,cols(2))
            Data1 = [Data1(1:i,:); NaN*ones(1,M); Data1((i+1):size(Data1,1),:)];
            %keyboard
            i=i+1;
        end
        i=i+1;
    end
end


figure(8); clf(8) %ax(1) = subplot(2,1,1);
ax(1)=subplot(3,1,1);
stairs(Data1(:,1),100*Data1(:,4),'k');
hold('on')
stairs(Data1(:,1),100*Data1(:,3),'r');
plot([Data1(1,1) Data1(end,1)],[0 0],'k:')
plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:')

set(gca,'xtick',[],'box','off')


load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1313.mat', 'Edi', 'Flow');
Time = 0:0.008:(length(Flow.values)-1)*0.008;


ax(2)=subplot(3,1,2); plot(Time,Flow.values,'g'); hold('on');
dsf=5; dt=Flow.interval;
FlowF=Flow.values;
if 1
    filter_HFcutoff_butter0 = 12.5;
    filter_order0 = 1;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
end
ax(2)=subplot(3,1,2); plot(downsample(Time,dsf),downsample(FlowF,dsf),'k');
set(gca,'xtick',[],'box','off')
ax(3)=subplot(3,1,3); plot(Time,Edi.values);
box('off');
linkaxes(ax,'x');

xlim([16200 16900])

%%

figure()
subplot(1,2,1);
mfl=100*Data1(:,4);
mfl_categories2_ = 100*[0:0.05:1.51];
for jj=2:length(mfl_categories2_)
    mfl_categories2(n,jj-1)=sum((mfl>=mfl_categories2_(jj-1)&mfl<mfl_categories2_(jj)))/length(mfl);
end

bar((mfl_categories2_(1:end-1)+0.025),mfl_categories2(n,:));
box('off');
xlim([0 150]);
ylim([0 0.14]);

subplot(1,2,2);
mfl=100*Data1(:,3);
mfl_categories2_ = 100*[0:0.05:1.51];
for jj=2:length(mfl_categories2_)
    mfl_categories2(n,jj-1)=sum((mfl>=mfl_categories2_(jj-1)&mfl<mfl_categories2_(jj)))/length(mfl);
end

bar((mfl_categories2_(1:end-1)+0.025),mfl_categories2(n,:));
box('off');
xlim([0 150]);
ylim([0 0.14]);
% a very simplistic way to see the tow trend lines
%plot(PtData.BB_time(PtData.PT==n),smooth(Gtest(PtData.PT==n),99),'c-');
%plot(PtData.BB_time(PtData.PT==n),smooth(Gpredicted(PtData.PT==n),99),'m-');
%end