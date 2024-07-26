

Mrange=1578%1:size(settings.Filenames,1);
LabelSpO2(Mrange,1)={NaN(1,1)};
ChannelNumberSpO2(Mrange,1)=NaN;


load('PrSpO2Models.mat');

HChannelSpreadsheet=settings.HChannelSpreadsheet;
HarmonizedChannelT=readtable(HChannelSpreadsheet);
SigList=HarmonizedChannelT.Properties.VariableNames;
ii=find(SigList=="SpO2");
SigToLoad=HarmonizedChannelT{:,ii};

%%
for n=Mrange

    try
        OutTall=table();
        clear Indright
        useregEDFread=0; %0=blockedfread

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
        % tying to remove other signal data that may qualify of PrSpo2
        possibleOtherSigInd=OutT.p95>105;
        OutT(possibleOtherSigInd,4:end)= array2table(nan(sum(possibleOtherSigInd==1),length(size(OutT,2)-4:size(OutT,2))));

        OutTall = [OutTall;OutT];
        OutTall.p05(isnan(OutTall.p05))=100*rand(sum(isnan(OutTall.p05)),1);
        OutTall.p25(isnan(OutTall.p25))=100*rand(sum(isnan(OutTall.p25)),1);
        OutTall.p50(isnan(OutTall.p50))=100*rand(sum(isnan(OutTall.p50)),1);
        OutTall.p75(isnan(OutTall.p75))=100*rand(sum(isnan(OutTall.p75)),1);
        OutTall.p95(isnan(OutTall.p95))=100*rand(sum(isnan(OutTall.p95)),1);

        mdlPrSpO2Use = mdlPrSpO2;
        OutTall.PrSpO2 = predict(mdlPrSpO2Use,OutTall);

        Indright = find(OutTall.PrSpO2>=thresX);
        if ~isempty(Indright)
            LabelPosible=OutTall.name(Indright);
            if sum(ismember(LabelPosible,SigToLoad))>0
                ind=Indright(ismember(LabelPosible,SigToLoad));
                if length(ind)==1
                    % only one signal detected
                    LabelSpO2{n,1}=OutTall.name(ind);
                    ChannelNumberSpO2(n,1)=find(strcmp(OutTall.name(ind),Label));
                    disp(['SpO2 signal found:' num2str(ChannelNumberSpO2(n,1)) LabelSpO2{n,1}])
                end
            end
        end
    end
end
%% save to excel sheet
Spo2col='AM';
xlswrite(AMasterSpreadsheet,ChannelNumberSpO2(Mrange),1,[Spo2col num2str(3+1+ Mrange(1)-1)]);

Spo2label='P';
xlswrite(AMasterSpreadsheet,(LabelSpO2{Mrange,1}),1,[Spo2label num2str(3+1+ Mrange(1)-1)]);



%% Build and Test model
buildmodel=0;
if buildmodel
    OutTall2 = OutTall;
    OutTall2.SpO2 = OutTall2.name=="SaO2";
    OutTall2.p05(isnan(OutTall2.p05))=100*rand(sum(isnan(OutTall2.p05)),1);
    OutTall2.p25(isnan(OutTall2.p25))=100*rand(sum(isnan(OutTall2.p25)),1);
    OutTall2.p50(isnan(OutTall2.p50))=100*rand(sum(isnan(OutTall2.p50)),1);
    OutTall2.p75(isnan(OutTall2.p75))=100*rand(sum(isnan(OutTall2.p75)),1);
    OutTall2.p95(isnan(OutTall2.p95))=100*rand(sum(isnan(OutTall2.p95)),1);

    %mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p50)*Fnotisnan','distribution','binomial');
    %mdlPrSpO2L = fitglm(OutTall2,'SpO2 ~ (p05+p25+p50+p75+p95)*Fnotisnan','distribution','binomial');
    mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p25+p75):(1+Fnotisnan)','distribution','binomial'); %optimized by SS, 6/8/2023
%     mdlPrSpO2 = fitglm(OutTall2,'SpO2 ~ (p25+p75+p95):(1+Fnotisnan)','distribution','binomial'); %optimized by SS, 6/8/2023

    mdlPrSpO2Use = mdlPrSpO2;

    OutTall2.PrSpO2 = predict(mdlPrSpO2Use,OutTall2);

    [thresX,AUC,SEM,~,~,SensSpec]=ROCAUCSEM(OutTall2.SpO2,OutTall2.PrSpO2,1);
    [PerfT,Raw,Cont] = PredictiveValue(OutTall2.SpO2,OutTall2.PrSpO2>=thresX,OutTall2.PrSpO2);
    Indwrong = find(OutTall2.SpO2~=(OutTall2.PrSpO2>=thresX))

    OutTall2errors = OutTall2(Indwrong,:);
    %%

    % save PrSpO2Models mdlPrSpO2 OutTall2
end
