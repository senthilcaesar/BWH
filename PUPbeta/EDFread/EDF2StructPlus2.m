
useregEDFread=0; %0=blockedfread
Mrange=1:2060;
%%
OutTall=table();

%%
for n=Mrange

tic
        if useregEDFread

            [Label,Transducer,Fs] = EDFChannelLabels(settings.SourceDir{n}); %now fails in some MESA EDFs [2/2021]

            %fix labels to make them legal variable names, may need to do more here.
            LabelTemp = erase(Label,' ');
            LabelTemp = replace(LabelTemp,'-','_');
            LabelTemp = replace(LabelTemp,'/','_');
            LabelTemp = replace(LabelTemp,'(','_');
            LabelTemp = erase(LabelTemp,')');
            Label = LabelTemp;
            for i=1:length(Label)
                n*1000 + i
                temp = readedfrev3(settings.SourceDir{n},i-1,0,Inf);
                SigS=setfield(SigS,Label{i},temp);
            end

        else
            SigS = struct();
            n
            clear signalCellX signalHeader
            [Header,signalHeader] = blockEdfLoad(settings.SourceDir{n});
            signalHeader=struct2table(signalHeader); %used if blockEDFLoad methods are used
            [~,~,signalCellX] = blockEdfLoad(settings.SourceDir{n},signalHeader.signal_labels);
            LabelTemp = erase(signalHeader.signal_labels,' ');
            LabelTemp = replace(LabelTemp,'-','_');
            LabelTemp = replace(LabelTemp,'/','_');
            LabelTemp = replace(LabelTemp,'(','_');
            LabelTemp = erase(LabelTemp,')');
            % caution: will handle only one duplicated channel here
            %fixed at line 59
            [~, ~, ic] = unique(LabelTemp,'stable');
            findDuplicateLabel=find(hist(ic,unique(ic,'stable'))>1);
            if ~isempty(findDuplicateLabel)
                %             LabelTemp{findDuplicateLabel}=strcat(LabelTemp{findDuplicateLabel},'_2');
                LabelTemp(findDuplicateLabel)=cellfun(@(c)[c '_2'],LabelTemp(findDuplicateLabel),'uni',false);
            end
            Label = LabelTemp';

            SigS = cell2struct(signalCellX, Label,2);
        end


OutT = DealNaNIntoTable(1,fieldnames(SigS));
OutT(1,:)=[];

OutT = table();
OutT.index = n*1000 + [1:length(fieldnames(SigS))]';
OutT.name = fieldnames(SigS);

for i=1:height(OutT)
  data = getfield(SigS,OutT.name{i});
  
  temp = data;
        temp(temp<0.4)=NaN;
        if nanmedian(temp)<1 && prctile(temp,95)<1 && prctile(temp,95)>0.5 && nanmedian(temp)>0.5
            data=data*100;
        end
        
  data(data<40)=NaN;
  OutT.Fnotisnan(i) = mean(~isnan(data));
  OutT.p05(i) = prctile(data,5);
  OutT.p25(i) = prctile(data,25);
  OutT.p50(i) = prctile(data,50);
  OutT.p75(i) = prctile(data,75);
  OutT.p95(i) = prctile(data,95);
end
%OutTall.SpO2=[];
OutTall = [OutTall;OutT];
toc
end

%% Build and Test model
OutTall2 = OutTall;
OutTall2.SpO2 = OutTall2.name=="SpO2";
OutTall2.p05(isnan(OutTall2.p05))=100*rand(sum(isnan(OutTall2.p05)),1);
OutTall2.p25(isnan(OutTall2.p25))=100*rand(sum(isnan(OutTall2.p25)),1);
OutTall2.p50(isnan(OutTall2.p50))=100*rand(sum(isnan(OutTall2.p50)),1);
OutTall2.p75(isnan(OutTall2.p75))=100*rand(sum(isnan(OutTall2.p75)),1);
OutTall2.p95(isnan(OutTall2.p95))=100*rand(sum(isnan(OutTall2.p95)),1);

%mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p50)*Fnotisnan','distribution','binomial');
%mdlPrSpO2L = fitglm(OutTall2,'SpO2 ~ (p05+p25+p50+p75+p95)*Fnotisnan','distribution','binomial');

Iremove = find(OutTall2.name=="HR" & OutTall2.p95==100);
OutTall2(Iremove,:)=[];
Iremove = find(OutTall2.name=="SpO2" & OutTall2.Fnotisnan==0);
OutTall2(Iremove,:)=[];

weights = ones(height(OutTall2),1);
W1 = sum(OutTall2.SpO2==1)/height(OutTall2)
weights(OutTall2.SpO2==1)=1./W1;
weights(OutTall2.SpO2==0)=1./(1-W1);
weights = weights./mean(weights);

OutTall2.p95under100 = 1.*(OutTall2.p95<=100);
%mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p25+p75):(1+Fnotisnan)','distribution','binomial','weights',weights); %optimized by SS, 6/8/2023
%mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p05+p25+p50+p75+p95):(1+Fnotisnan)','distribution','binomial','weights',weights); %optimized by SS, 6/8/2023
%mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p05+p25+p50+p75+p95):(1+Fnotisnan+Fnotisnan^2)','distribution','binomial','weights',weights); %optimized by SS, 6/8/2023

%mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p25*(p50+p75+p95) + p50*(p75+p95) + p75*p95):(Fnotisnan)','distribution','binomial','weights',weights); %optimized by SS, 6/8/2023
mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p05*(p25+p50+p75+p95) + p25*(p50+p75+p95) + p50*(p75+p95) + p75*p95):(Fnotisnan)','distribution','binomial','weights',weights); %optimized by SS, 6/8/2023

mdlPrSpO2Use = mdlPrSpO2; 

OutTall2.PrSpO2 = predict(mdlPrSpO2Use,OutTall2); 

[thresX,AUC,SEM,~,~,SensSpec]=ROCAUCSEM(OutTall2.SpO2,OutTall2.PrSpO2,1); 
[PerfT,Raw,Cont] = PredictiveValue(OutTall2.SpO2,OutTall2.PrSpO2>=thresX,OutTall2.PrSpO2); 
Indwrong = find(OutTall2.SpO2~=(OutTall2.PrSpO2>=thresX));
IndFN = find(OutTall2.SpO2==1 & OutTall2.PrSpO2<thresX);

OutTall2errors = OutTall2(Indwrong,:);

%%
save PrSpO2Workspace2
save PrSpO2Models mdlPrSpO2 OutTall2 thresX



