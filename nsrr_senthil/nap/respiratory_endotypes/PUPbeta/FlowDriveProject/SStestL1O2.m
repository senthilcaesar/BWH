clear all
close all
clc

load('SA_FD_16.mat')
load('SSblend_FS_FD_56.mat')
%%
clearvars -except PtData FeatureNames Amatrix

  
%% Rem fluttering

if 0
FtrsToExclude = {'InspFlutPow_Vpeak2_21_O','ExpFlutPow_VpeakE2_22_O','InspFlutPow_Vpeak2_20_T','ExpFlutPow_VpeakE2_21_T','InspExpFlutPow_22_T'};
Ind = [];
for i=1:length(FtrsToExclude)
    Ind(i)=find(strcmp(FeatureNames.Name,FtrsToExclude(i)));
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

%% Find and sort based on univariate linear correlations
Gtest_ = PtData{:,13}; %Keep

maxG = 1.5;
    Gtest_(Gtest_>maxG)=maxG;

colofones = ones(length(Gtest_),1);
predyL1O = NaN*Gtest_;

%% Weights

clear Ndata
dx=0.2; xbins=[0 0.3:dx:0.9 maxG];

%xbins=[0 0.1:0.05:1.1 1.5];
for i=1:length(xbins)-1
    Ix=Gtest_>xbins(i)&Gtest_<=xbins(i+1);
    Ndata(i)=sum(Ix);
end
%Ndata = Ndata.^0.5;
weightsbins = 1./Ndata;
weightsbins = [2 1 1 1 0.5];
weightsbins = weightsbins/mean(weightsbins);
%weightsbins = weightsbins/weightsbins(end);


weights = NaN*Gtest_;
for i=1:length(xbins)-1
    Ix=Gtest_>xbins(i)&Gtest_<=xbins(i+1);
    weights(Ix)=weightsbins(i);
end
weights = weights/nanmean(weights);

useweights=0;
if ~useweights %overwrite, use no weights
    weights = ones(length(weights),1);
end

%%
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

clear Rsq
for i=1:size(Amatrix2,2)
    for j=1:i-1
        
        %[b,bint,r,rint,stats] = regress(Amatrix2(:,j),[colofones(:) Amatrix2(:,i)]);
        %Rsq(i,j)=stats(1);
       Rsq(i,j)= corr(Amatrix2(:,j),Amatrix2(:,i)).^2;
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
alpha=0.0001; 
MaxNfeatures=25;
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
    w=w/mean(w);
    
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
    w=w/mean(w);
    colofones_(nanx|nany)=[];
    warning('off')
    maxp=Inf;
    if Fwd==0 %backwards elimination
        If = 1:size(Amatrix2,2);
        while length(If)>=MaxNfeatures || maxp>alpha
            if 1
                [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2_(:,If),Gtest,w); %faster
                Pvals(1)=[];
            elseif 0
                [b,dev,stats]=glmfit(Amatrix2_(:,If),Gtest,'normal','weights',w);
                Pvals = stats.p(2:end);
            end
            [maxp,maxpi]=max(Pvals);
            if maxp<alpha && (length(If)-1)<=MaxNfeatures
                break
            end
            disp(['Removing:p=' num2str(maxp)])
            If(maxpi)=[];
            Labels(maxpi)=[];
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
                [~,Pvals,~,~]=glmfitFast(Amatrix2_(:,[If notinc(mini)]),Gtest,w);
                if Pvals(end)>alpha
                    disp('end FSR');
                    break
                end
                If = [If notinc(mini)];
                notinc(mini)=[];
                Labels_(If)
                %check for remove
                [~,Pvals,RMSE,~]=glmfitFast(Amatrix2_(:,If),Gtest,w);
                Pvals(1)=[];
                while max(Pvals)>2*alpha
                [~,remi]=max(Pvals);
                notinc = [notinc If(remi)];
                If(remi)=[];
                [~,Pvals,RMSE,~]=glmfitFast(Amatrix2_(:,If),Gtest,w);
                Pvals(1)=[];
                end
            end
        end
        [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2_(:,If),Gtest,w);
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
    
    predyL1O(Isubj) = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
    Err(subj) = nanmean(abs(predyL1O(Isubj)-Gtest_(Isubj)));
    ErrRms(subj) = nanmean((predyL1O(Isubj)-Gtest_(Isubj)).^2).^0.5;
    
    Nfeatures(subj) = length(Labels2);
    %keyboard
end



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
    
    RsqL1O = 1-nansum(weights.*(Gtest_-predyL1O).^2)/nansum(weights.*(Gtest_-nanmean(Gtest_)).^2);
    Rvalue = RsqL1O^0.5;
    ErrL1O = nanmean(abs(predyL1O-Gtest_));
    ErrL1Orms = nanmean((predyL1O-Gtest_).^2).^0.5;
    
    predy = predyL1O;% predyL1O_fsmr;
    Yval = Gtest_;
    
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
%%

%Plot FL values over time for an example subject (e.g. compare with Spike data). Random check (n=8) Looks great.
n=8; % 1313?
%for n=1:30
ypredAll = predict(SVMModel_reg,Amatrix2(:,If));
Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) ypredAll(PtData.PT==n) Gtest(PtData.PT==n)];

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


ax(2)=subplot(3,1,2); plot(Time,Flow.values);
set(gca,'xtick',[],'box','off')
ax(3)=subplot(3,1,3); plot(Time,Edi.values);
box('off');
linkaxes(ax,'x');



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