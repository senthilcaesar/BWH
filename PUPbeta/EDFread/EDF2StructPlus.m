

Mrange=6:50;
%%
OutTall=table();

%%
for n=Mrange
[Label,Transducer,Fs] = EDFChannelLabels(settings.SourceDir{n}); %now fails in some MESA EDFs [2/2021]

%fix labels to make them legal variable names, may need to do more here.
LabelTemp = erase(Label,' ');
LabelTemp = replace(LabelTemp,'-','_');
Label = LabelTemp;

SigS = struct();

for i=1:length(Label)
    n*1000 + i
temp = readedfrev3(settings.SourceDir{n},i-1,0,Inf);
SigS=setfield(SigS,Label{i},temp);
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
mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p25+p75):(1+Fnotisnan)','distribution','binomial'); %optimized by SS, 6/8/2023

mdlPrSpO2Use = mdlPrSpO2; 

OutTall2.PrSpO2 = predict(mdlPrSpO2Use,OutTall2); 

[thresX,AUC,SEM,~,~,SensSpec]=ROCAUCSEM(OutTall2.SpO2,OutTall2.PrSpO2,1); 
[PerfT,Raw,Cont] = PredictiveValue(OutTall2.SpO2,OutTall2.PrSpO2>=thresX,OutTall2.PrSpO2); 
Indwrong = find(OutTall2.SpO2~=(OutTall2.PrSpO2>=thresX)) 

OutTall2errors = OutTall2(Indwrong,:);

%%

save PrSpO2Models mdlPrSpO2 OutTall2 thresX