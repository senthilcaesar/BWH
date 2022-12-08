function OutcomeT = getFLandSnore(MrangeOverride)
global settings
%errors: SnoreDB95thcentile is incorrect, Etype is missing
% %
% settings=ImportSettings(settings,AMasterSpreadsheet);
% 
% [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
% NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
% analyzelist = logical(num(:,2));
% settings.invertflowlist = logical(num(:,1));
% settings.protocol = patients(:,3);
% settings.patients = patients;

M=size(settings.patients,1);
Mrange=1:M;
if exist('MrangeOverride')
    Mrange=MrangeOverride;
end

dir=[settings.workdir 'Analyzed\'];
%dir=[settings.OutputDataDirectory 'EventAnalyzed\'];
OutcomeT = table([]);

minflowdrive = 0.1;

for i=1:length(Mrange) %1:11 13:19 fails: 20 23 24 33
    clear BreathDataTable BreathFLDataTable
    try
        success(i)=1;
        filen = [dir settings.savename '_' num2str(Mrange(i))];
        load(filen,'BreathDataTable','BreathFLDataTable');
        
        if 1
            [~,~,~,~,T]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
        else
            T = BreathDataTableFulls;
        end
        
        %QC on FlowDrive
        T.FlowDrive(T.FlowDrive<minflowdrive)=minflowdrive;
        T.FlowDrive(T.FlowDrive>1.2)=1.2;
        T.FlowDrive(T.ApneaB==1)=minflowdrive;
        T.FlowDrive(T.VI<minflowdrive & T.Ecentralapnea==0)=minflowdrive;
        
        %    a = 1*(T.AR==1 & T.hypnog_B==4)
        %    win=zeros(height(T),1);
        %    T.distfromsleepB = howfarawayfromsleepInf(a,win);
        %
        criteriaW = T.AR==1 & T.hypnog_B==4;
        criteriaW = T.AR==1 & T.hypnog_B==4 & ~isnan(T.spo2);% & T.distfromsleepB<10;
        sum(criteriaW);
        OutcomeT.ID(i)=Mrange(i);
        OutcomeT.FlowDriveMedianW(i) = nanmedian(T.FlowDrive(criteriaW));
        
        criteriaFL = T.AR==0 & T.hypnog_B<4 & T.hypnog_B>=0;
        sum(criteriaFL);
        
        disp(strcat("Median FlowDrive = ",string(nanmedian(T.FlowDrive(criteriaFL))*100)));
        
        OutcomeT.FlowDriveMedianSleepNoAr(i) = nanmedian(T.FlowDrive(criteriaFL));
        OutcomeT.FLfrequency50SleepNoAr(i) = nanmean(T.FlowDrive(criteriaFL)<0.5)*100;
        OutcomeT.FLfrequency70SleepNoAr(i) = nanmean(T.FlowDrive(criteriaFL)<0.7)*100;
        OutcomeT.FLfrequency30SleepNoAr(i) = nanmean(T.FlowDrive(criteriaFL)<0.3)*100;
        OutcomeT.FLfrequency40SleepNoAr(i) = nanmean(T.FlowDrive(criteriaFL)<0.4)*100;
        OutcomeT.FLfrequency60SleepNoAr(i) = nanmean(T.FlowDrive(criteriaFL)<0.6)*100;
        
        criteriaFL2 = T.hypnog_B<4 & T.hypnog_B>=0; %All Sleep
            OutcomeT.FlowDriveMedianSleepAll(i) = nanmedian(T.FlowDrive(criteriaFL2));
            OutcomeT.FLfrequency50SleepAll(i) = nanmean(T.FlowDrive(criteriaFL2)<0.5)*100;
            OutcomeT.FLfrequency70SleepAll(i) = nanmean(T.FlowDrive(criteriaFL2)<0.7)*100;
            OutcomeT.FLfrequency30SleepAll(i) = nanmean(T.FlowDrive(criteriaFL2)<0.3)*100;
            OutcomeT.FLfrequency40SleepAll(i) = nanmean(T.FlowDrive(criteriaFL2)<0.4)*100;
            OutcomeT.FLfrequency60SleepAll(i) = nanmean(T.FlowDrive(criteriaFL2)<0.6)*100;
            
        OutcomeT.RIPcorrMedian(i) = nanmedian(T.RIPcorr(criteriaFL));
        OutcomeT.RIPcorrMean(i) = nanmean(T.RIPcorr(criteriaFL));
        
        if sum(T.Properties.VariableNames=="SnoreDBmean")==1
            criteriaSn = T.hypnog_B<4 & T.hypnog_B>=0 & ~isnan(T.SnoreDBmean) & T.SnoreDBmean>50;
            sum(criteriaSn);
            %All sleep epoch data, regardless of arousals
            SnoreRMSsq = 10.^((T.SnoreDBmean(criteriaSn)/10))*(0.00002^2);
            T.Ttot = T.Time_end - T.Time_start;
            %time weighted average:
            SnoreRMSsq_mean = nansum(SnoreRMSsq.*T.Ttot(criteriaSn))./nansum(T.Ttot(criteriaSn));
            SnoreDB_mean=10*(log10(SnoreRMSsq_mean/(0.00002^2))); %convert back to dB after averaging
            OutcomeT.SnoreDB_mean(i) = SnoreDB_mean;
            
            disp(strcat("Snore Mean dB = ",string(SnoreDB_mean)));
            
            OutcomeT.FSnore80DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>80)*100;
            OutcomeT.FSnore85DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>85)*100;
            OutcomeT.FSnore90DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>90)*100;
            OutcomeT.FSnore95DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>95)*100;
            OutcomeT.FSnore100DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>100)*100;
            
            OutcomeT.FSnore110DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>110)*100;
            
            OutcomeT.FSnore120DB_90c(i) = mean(T.SnoreDB90thcentile(criteriaSn)>120)*100;
            
            OutcomeT.FSnore80DB_95c(i) = mean(T.SnoreDB95thcentile(criteriaSn)>80)*100;
            OutcomeT.FSnore85DB_95c(i) = mean(T.SnoreDB95thcentile(criteriaSn)>85)*100;
            OutcomeT.FSnore90DB_95c(i) = mean(T.SnoreDB95thcentile(criteriaSn)>90)*100;
            OutcomeT.FSnore95DB_95c(i) = mean(T.SnoreDB95thcentile(criteriaSn)>95)*100;
            OutcomeT.FSnore100DB_95c(i) = mean(T.SnoreDB95thcentile(criteriaSn)>100)*100;
            
            OutcomeT.FSnore110DB_95c(i) = mean(T.SnoreDB95thcentile(criteriaSn)>110)*100;
            
            OutcomeT.SnoringMedianW(i) = nanmedian(T.SnoreDBmean(criteriaW));
            OutcomeT.SnoreDBminMedian(i) = nanmedian(T.SnoreDBmin(criteriaSn));
            
            
            OutcomeT.FSnore90DB_Leq(i) = mean(T.SnoreDBmean(criteriaSn)>90)*100;
            OutcomeT.FSnore85DB_Leq(i) = mean(T.SnoreDBmean(criteriaSn)>85)*100;
            OutcomeT.FSnore80DB_Leq(i) = mean(T.SnoreDBmean(criteriaSn)>80)*100;
            OutcomeT.FSnore75DB_Leq(i) = mean(T.SnoreDBmean(criteriaSn)>75)*100;
            
        else
            OutcomeT.SnoreDB_mean(i) = NaN;
            OutcomeT.FSnore80DB_90c(i) = NaN;
            OutcomeT.FSnore85DB_90c(i) = NaN;
            OutcomeT.FSnore90DB_90c(i) = NaN;
            OutcomeT.FSnore95DB_90c(i) = NaN;
            OutcomeT.FSnore100DB_90c(i) = NaN;

            OutcomeT.FSnore110DB_90c(i) = NaN;

            OutcomeT.FSnore120DB_90c(i) = NaN;

            OutcomeT.FSnore80DB_95c(i) = NaN;
            OutcomeT.FSnore85DB_95c(i) = NaN;
            OutcomeT.FSnore90DB_95c(i) = NaN;
            OutcomeT.FSnore95DB_95c(i) = NaN;
            OutcomeT.FSnore100DB_95c(i) = NaN;

            OutcomeT.FSnore110DB_95c(i) = NaN;

            OutcomeT.SnoringMedianW(i) = NaN;
            OutcomeT.SnoreDBminMedian(i) = NaN;


            OutcomeT.FSnore90DB_Leq(i) = NaN;
            OutcomeT.FSnore85DB_Leq(i) = NaN;
            OutcomeT.FSnore80DB_Leq(i) = NaN;
            OutcomeT.FSnore75DB_Leq(i) = NaN;

            %I = GiantTable;
           
        end
        
        if sum(T.Properties.VariableNames=="SnoreDBNoxmean")==1
            criteriaSn = T.hypnog_B<4 & T.hypnog_B>=0 & ~isnan(T.SnoreDBNoxmean) & T.SnoreDBNoxmean>50;
            sum(criteriaSn);
            SnoreRMSNoxsq = 10.^((T.SnoreDBNox95thcentileI(criteriaSn)/10))*(0.00002^2);
            T.Ttot = T.Time_end - T.Time_start;
            %time weighted average:
            I = ~isnan(SnoreRMSNoxsq);
            SnoreRMSNoxsq_mean = nansum(SnoreRMSNoxsq(I).*T.Ttot(I))./nansum(T.Ttot(I));
            SnoreDBNox_mean=10*(log10(SnoreRMSNoxsq_mean/(0.00002^2))); %convert back to dB after averaging

            OutcomeT.SnoreDBNox_mean(i) = SnoreDBNox_mean;

            disp(strcat("Snore Mean dB = ",string(SnoreDBNox_mean)));

            OutcomeT.FSnoreNox70DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>70)*100;
            OutcomeT.FSnoreNox75DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>75)*100;
            OutcomeT.FSnoreNox80DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>80)*100;
            OutcomeT.FSnoreNox85DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>85)*100;
            OutcomeT.FSnoreNox90DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>90)*100;

            OutcomeT.FSnoreNox100DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>100)*100;

            OutcomeT.FSnoreNox110DB_90c(i) = mean(T.SnoreDBNox90thcentileI(criteriaSn)>110)*100;

            OutcomeT.FSnoreNox70DB_95c(i) = mean(T.SnoreDBNox95thcentileI(criteriaSn)>70)*100;
            OutcomeT.FSnoreNox75DB_95c(i) = mean(T.SnoreDBNox95thcentileI(criteriaSn)>75)*100;
            OutcomeT.FSnoreNox80DB_95c(i) = mean(T.SnoreDBNox95thcentileI(criteriaSn)>80)*100;
            OutcomeT.FSnoreNox85DB_95c(i) = mean(T.SnoreDBNox95thcentileI(criteriaSn)>85)*100;
            OutcomeT.FSnoreNox90DB_95c(i) = mean(T.SnoreDBNox95thcentileI(criteriaSn)>90)*100;

            OutcomeT.FSnoreNox100DB_95c(i) = mean(T.SnoreDBNox95thcentileI(criteriaSn)>100)*100;

            OutcomeT.SnoringNoxMedianW(i) = nanmedian(T.SnoreDBNoxmeanI(criteriaW));
            OutcomeT.SnoreDBNoxminMedian(i) = nanmedian(T.SnoreDBNoxminI(criteriaSn));


            OutcomeT.FSnoreNox85DB_Leq(i) = mean(T.SnoreDBNoxmeanI(criteriaSn)>85)*100;
            OutcomeT.FSnoreNox80DB_Leq(i) = mean(T.SnoreDBNoxmeanI(criteriaSn)>80)*100;
            OutcomeT.FSnoreNox75DB_Leq(i) = mean(T.SnoreDBNoxmeanI(criteriaSn)>75)*100;
            OutcomeT.FSnoreNox70DB_Leq(i) = mean(T.SnoreDBNoxmeanI(criteriaSn)>70)*100;

        else
            criteriaSn = T.hypnog_B<4 & T.hypnog_B>=0 & ~isnan(T.SnoreDBNoxmean) & T.SnoreDBNoxmean>50;
            sum(criteriaSn);
            SnoreRMSNoxsq = 10.^((T.SnoreDBNox95thcentileI(criteriaSn)/10))*(0.00002^2);
            T.Ttot = T.Time_end - T.Time_start;
            %time weighted average:
            I = ~isnan(SnoreRMSNoxsq);
            SnoreRMSNoxsq_mean = nansum(SnoreRMSNoxsq(I).*T.Ttot(I))./nansum(T.Ttot(I));
            SnoreDBNox_mean=10*(log10(SnoreRMSNoxsq_mean/(0.00002^2))); %convert back to dB after averaging

            OutcomeT.SnoreDBNox_mean(i) = SnoreDBNox_mean;

            disp(strcat("Snore Mean dB = ",string(SnoreDBNox_mean)));

            OutcomeT.FSnoreNox70DB_90c(i) = NaN;
            OutcomeT.FSnoreNox75DB_90c(i) = NaN;
            OutcomeT.FSnoreNox80DB_90c(i) = NaN;
            OutcomeT.FSnoreNox85DB_90c(i) = NaN;
            OutcomeT.FSnoreNox90DB_90c(i) = NaN;

            OutcomeT.FSnoreNox100DB_90c(i) = NaN;

            OutcomeT.FSnoreNox110DB_90c(i) = NaN;

            OutcomeT.FSnoreNox70DB_95c(i) = NaN;
            OutcomeT.FSnoreNox75DB_95c(i) = NaN;
            OutcomeT.FSnoreNox80DB_95c(i) = NaN;
            OutcomeT.FSnoreNox85DB_95c(i) = NaN;
            OutcomeT.FSnoreNox90DB_95c(i) = NaN;

            OutcomeT.FSnoreNox100DB_95c(i) = NaN;

            OutcomeT.SnoringNoxMedianW(i) = NaN;
            OutcomeT.SnoreDBNoxminMedian(i) = NaN;


            OutcomeT.FSnoreNox85DB_Leq(i) = NaN;
            OutcomeT.FSnoreNox80DB_Leq(i) = NaN;
            OutcomeT.FSnoreNox75DB_Leq(i) = NaN;
            OutcomeT.FSnoreNox70DB_Leq(i) = NaN;
        end
        disp(strcat("Subject ",string(i)," completed"));
    catch me
        success(i)=0;
        disp(strcat("Subject ",string(i)," failed"));
    end
end
OutcomeT.Var1=[];
OutcomeT{success==0,2:end}=NaN;


%%
if 0
    
    cd(settings.AMasterdir)
    if 0
        OutcomeT2 = OutcomeT
        save OutcomeT2 OutcomeT2
    end
    
    load Data
    save OutcomeT OutcomeT
    load OutcomeT OutcomeT
    
    %Ifail = find(success==0);
    %OutcomeT{Ifail,:}=NaN;
    
    %OutcomeT_ = OutcomeT;
    %OutcomeT = OutcomeT_;
end

