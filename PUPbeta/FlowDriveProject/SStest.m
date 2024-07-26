clear all
close all
clc

load('SA_FD_16.mat')

%%

Gtest = PtData{:,13};

downsamplefactor=3; %downsample 10 takes just 2s [similar for 1-20 features],5->11s, 4->15s, 3->25s, none->3-5 min
    %I=ceil(rand(1,1)*downsamplefactor):downsamplefactor:length(Gtest); 
    I=1:downsamplefactor:length(Gtest); 
    
    NotI = 1:size(Amatrix,1);
        NotI(I)=[];
        
        

%% Names
temp = FeatureNames.Name;
clear temp1 temp2 temp3
for i=1:length(temp)
    temp1{i}=[temp{i} '_1p0'];
    temp2{i}=[temp{i} '_0p5'];
    temp3{i}=[temp{i} '_2p0'];
end

%% Make large matrix including square-root and square transformed data where possible

Nonzero = sum(Amatrix>0)./sum(~isnan(Amatrix));
crit = Nonzero<0.999;
AmatrixSQRT = Amatrix.^0.5;
    AmatrixSQRT(Amatrix<0)=0;
    AmatrixSQRT(:,crit)=[];
    
AmatrixSQ = Amatrix.^2;
    AmatrixSQ(Amatrix<0)=0;
    AmatrixSQ(:,crit)=[];
    
    Amatrix2 = [Amatrix AmatrixSQ AmatrixSQRT];
    
    temp2(crit)=[];
    temp3(crit)=[];
    Labels = [temp1 temp2 temp3]';
    
%% Find and sort based on univariate linear correlations
subj=1
Gtest_ = Gtest;
I=(PtData.PT==1);

colofones = 1 + 0*Gtest;
clear Rsq
    
for i=1:size(Amatrix2,2)
    [b,bint,r,rint,stats] = regress(Gtest(I),[colofones(I) Amatrix2(I,i)]);
    Rsq(i)=stats(1);
end

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

%% Forward Stepwise Linear Feature Preselection

%dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
xbins=[0 0.1:0.05:1.5];
for i=1:length(xbins)-1
    Ix=Gtest>xbins(i)&Gtest<xbins(i+1);
    Ndata(i)=sum(Ix);
end
weightsbins = 1./Ndata;
for i=1:length(xbins)-1
    Ix=Gtest>xbins(i)&Gtest<xbins(i+1);
    weights(Ix)=weightsbins(i);
end    
weights = weights/sum(weights);
    %weights = 1+0*weights;
    
%Amatrix_ = Amatrix(:,ItopN);
%[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit([Amatrix],Gtest,'maxiter',15);

%replace this with a stepwise algorithm that takes weights
ChooseN=40;
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit([Amatrix2(I,:)],Gtest(I),'maxiter',ChooseN);
Nvarsselected = sum(inmodel)

predy = Amatrix2(I,inmodel==1)*b(inmodel==1) + stats.intercept;
Rsq = 1-nansum(weights(I).*(Gtest(I)-predy).^2)/nansum(weights(I).*(Gtest(I)-nanmean(Gtest(I))).^2);
Rvalue = Rsq^0.5

if 1
predy = Amatrix2(NotI,inmodel==1)*b(inmodel==1) + stats.intercept;
Rsq = 1-nansum(weights(NotI).*(Gtest(NotI)-predy).^2)/nansum(weights(NotI).*(Gtest(NotI)-nanmean(Gtest(NotI))).^2);
Rvalue_ = Rsq^0.5
end
%%
If = [];
for i=1:sum(inmodel)
    temp = find(history.in(i,:));
    x=sum(temp==If');
    temp(x==1)=[];
    If(i)=temp;
end

Labels = Labels(If);

clear Rvalue Rvalue_
figure(1); clf(1);
for i=1:sum(inmodel) %no need to do this one by one but makes it easy to assess progressive improvement
    %subplot(5,10,i);
    [b,dev,stats]=glmfit(Amatrix2(I,If(1:i)),Gtest(I),'normal','weights',weights(I));
    predy = [colofones(I) Amatrix2(I,If(1:i))]*b;
    Rsq = 1-nansum(weights(I).*(Gtest(I)-predy).^2)/nansum(weights(I).*(Gtest(I)-nanmean(Gtest(I))).^2);
    Rvalue(i) = Rsq^0.5
    predy = [colofones(NotI) Amatrix2(NotI,If(1:i))]*b;
    Rsq = 1-nansum(weights(NotI).*(Gtest(NotI)-predy).^2)/nansum(weights(NotI).*(Gtest(NotI)-nanmean(Gtest(NotI))).^2);
    Rvalue_(i) = Rsq^0.5
    %scatter(Amatrix2(:,I(i)),Gtest,2,'filled','markerfacealpha',0.1)
    %scatter(predy,Gtest(NotI),2,'filled','markerfacealpha',0.1)
    %xlim([0 1.25]);
end

if 0 %add additional features based on univariate
II=find(inmodel==0);
inmodel(II(1:3))=1;
end

%%
%Remove the transforms for the SVM model

for n=1:3
i=length(Labels);
while i>1
    for j=1:length(Labels)
        if i==j
            continue
        end
        if strcmp(Labels{i}(1:end-4),Labels{j}(1:end-4))
           Labels(i)=[];
           If(i)=[];
           i=i-1;
           break
        end
    end
    i=i-1;
end
end

%% Build SVM model based on set features (skip to assess linear regression model)
%weights = 1+0*weights;

    %Run SVM regress
    tic
    svpmethodstr = 'SMO';
    kernelf = 'rbf';%'polynomial', rbf
    sigma=3;
    SVMModel_reg = fitrsvm(Amatrix2(I,If),Gtest(I),'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelf,'KernelScale',sigma,'Weights',weights(I));%'OptimizeHyperparameters','auto' ,'Weights',weights
    toc
    
    tic
    predyTrain = predict(SVMModel_reg,Amatrix2(I,If));
    %predy = NaN*Gtest;
    predy = predict(SVMModel_reg,Amatrix2(NotI,If));
    Nnan = sum(isnan(predy))
    toc
    
    RsqTrain = 1-nansum(weights(I).*(Gtest(I)-predyTrain).^2)/nansum(weights(I).*(Gtest(I)-nanmean(Gtest(I))).^2);
RvalueTrain = RsqTrain^0.5

Rsq = 1-nansum(weights(NotI).*(Gtest(NotI)-predy).^2)/nansum(weights(NotI).*(Gtest(NotI)-nanmean(Gtest(NotI))).^2);
Rvalue = Rsq^0.5
%test set, training values are NaN

%% Custom colormap

if ~exist('customcmap')
    customcmap = [...
    0.2422    0.1504    0.6603 ; ...
    0.2504    0.1650    0.7076 ; ...
    0.2578    0.1818    0.7511 ; ...
    0.2647    0.1978    0.7952 ; ...
    0.2706    0.2147    0.8364 ; ...
    0.2751    0.2342    0.8710 ; ...
    0.2783    0.2559    0.8991 ; ...
    0.2803    0.2782    0.9221 ; ...
    0.2813    0.3006    0.9414 ; ...
    0.2810    0.3228    0.9579 ; ...
    0.2795    0.3447    0.9717 ; ...
    0.2760    0.3667    0.9829 ; ...
    0.2699    0.3892    0.9906 ; ...
    0.2602    0.4123    0.9952 ; ...
    0.2440    0.4358    0.9988 ; ...
    0.2206    0.4603    0.9973 ; ...
    0.1963    0.4847    0.9892 ; ...
    0.1834    0.5074    0.9798 ; ...
    0.1786    0.5289    0.9682 ; ...
    0.1764    0.5499    0.9520 ; ...
    0.1687    0.5703    0.9359 ; ...
    0.1540    0.5902    0.9218 ; ...
    0.1460    0.6091    0.9079 ; ...
    0.1380    0.6276    0.8973 ; ...
    0.1248    0.6459    0.8883 ; ...
    0.1113    0.6635    0.8763 ; ...
    0.0952    0.6798    0.8598 ; ...
    0.0689    0.6948    0.8394 ; ...
    0.0297    0.7082    0.8163 ; ...
    0.0036    0.7203    0.7917 ; ...
    0.0067    0.7312    0.7660 ; ...
    0.0433    0.7411    0.7394 ; ...
    0.0964    0.7500    0.7120 ; ...
    0.1408    0.7584    0.6842 ; ...
    0.1717    0.7670    0.6554 ; ...
    0.1938    0.7758    0.6251 ; ...
    0.2161    0.7843    0.5923 ; ...
    0.2470    0.7918    0.5567 ; ...
    0.2906    0.7973    0.5188 ; ...
    0.3406    0.8008    0.4789 ; ...
    0.3909    0.8029    0.4354 ; ...
    0.4456    0.8024    0.3909 ; ...
    0.5044    0.7993    0.3480 ; ...
    0.5616    0.7942    0.3045 ; ...
    0.6174    0.7876    0.2612 ; ...
    0.6720    0.7793    0.2227 ; ...
    0.7242    0.7698    0.1910 ; ...
    0.7738    0.7598    0.1646 ; ...
    0.8203    0.7498    0.1535 ; ...
    0.8634    0.7406    0.1596 ; ...
    0.9035    0.7330    0.1774 ; ...
    0.9393    0.7288    0.2100 ; ...
    0.9728    0.7298    0.2394 ; ...
    0.9956    0.7434    0.2371 ; ...
    0.9970    0.7659    0.2199 ; ...
    0.9952    0.7893    0.2028 ; ...
    0.9892    0.8136    0.1885 ; ...
    0.9786    0.8386    0.1766 ; ...
    0.9676    0.8639    0.1643 ; ...
    0.9610    0.8890    0.1537 ; ...
    0.9597    0.9135    0.1423 ; ...
    0.9628    0.9373    0.1265 ; ...
    0.9691    0.9606    0.1064 ; ...
    0.9769    0.9839    0.0805 ; ... 
    ];
end
    
    %% Analysis and plots
    
    Yval = Gtest(NotI);
    
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
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01)
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01)

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
    pcolor(XGrid,YGrid,C2)
    hold on                             % Plot original data points
    [X,Y] = meshgrid(x,y);
    colormap(gcf,customcmap)
    h=colorbar();
        set(h,'FontName','Arial Narrow','Limits',[0 max(max(C2))]);
    C1 = flipud(C1); 
    C1 = fliplr(C1);
    C1 = C1'; 
    for i=1:size(C1,1)
        for j=1:size(C1,2)
            if C1(i,j)<25%((max(max(C1))-min(min(C1)))/2+min(min(C1)))
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
    %classcutoffs = [0.7 0.4];
    classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
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
    labeltext={'Normal','Mild','FL'};
    yticks(y); yticklabels(gca,fliplr(labeltext));
    xticks(x); xticklabels(gca,fliplr(labeltext));
    
    end
    
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