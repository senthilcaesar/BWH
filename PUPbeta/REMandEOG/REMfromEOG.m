function [PrREM,Tbl,EOGrmsConj_p50_180sb,ymtaAr180b,t_,rmssumbyrbyar,rmssumbyr,rmssum,Thres,LOCfilt,ROCfilt] = REMfromEOG(LOC,ROC,Time,EventsAr,Exclude,Options,mdlREM)

if isempty(Options) 
        filter_LFcutoff_butter1 = 0.25 %0.25, 0.67 at training
        filter_HFcutoff_butter1 = 7; %7, 4,0.3        
        widthrms = 6;
        widthmta = 120;
        FAr = 0.9999;
        ThresSig = 0.5; %R value threshold for defining conjugate movements
else
        filter_LFcutoff_butter1 = Options(1); %4,0.3
        filter_HFcutoff_butter1 = Options(2);
        widthrms = Options(3);
        widthmta = Options(4);
        FAr = Options(5);
        ThresSig = Options(6);
end        
        logitinverse = @(p) 1./(1+exp(-p));
        
        dT = Time(2) - Time(1);
        
        filter_order0 = 2;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dT/2));
        LOCfilt = nanfilter(B_butter0,A_butter0,LOC,1);
        ROCfilt = nanfilter(B_butter0,A_butter0,ROC,1);
        
        StepSize1 = 0.25;
        BufferI = buffer(1:length(LOC),round(widthrms*1/dT),round((widthrms-StepSize1)*1/dT),'nodelay');
                BufferI(:,sum(BufferI==0)>0)=[];
        
        t = Time(BufferI);% buffer(DataEventHypnog_Mat_ds.Time,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
        t_ = t(1,:)' + widthrms/2;
        
        LOCbuffer = LOCfilt(BufferI); %buffer(LOCfilt,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay'); 
        ROCbuffer = ROCfilt(BufferI); 
        if 0
            LOCrms = (nanmean((LOCbuffer.^2)).^0.5)';
            ROCrms = (nanmean((ROCbuffer.^2)).^0.5)'; %rms of yR
        else
            LOCrms = nanmedian(abs(LOCbuffer))'; %abs is instantaneous rms, median will track dynamic state shifts better
            ROCrms = nanmedian(abs(ROCbuffer))'; 
        end
        
        ar = EventsAr(BufferI) | Exclude(BufferI); 
        ar(isnan(ar))=1;
            ar = (nanmean(ar)')>0.5; %0.05

        rval=corrfast(LOCbuffer,ROCbuffer)';  

        rmssum=(LOCrms + ROCrms);
        rmssumbyr=rmssum.*logitinverse((-rval-ThresSig)*10);
        rmssumbyrbyar=rmssum.*logitinverse((-rval-ThresSig)*10).*(1-FAr*(ar));
        
        %BufferI2 = buffer(1:length(t_),round(widthmta*1/(t_(2)-t_(1))),round(widthmta*1/(t_(2)-t_(1)))-1,'nodelay');
        BufferI2 = buffer(1:length(t_),round(widthmta*1/(t_(2)-t_(1))),round((widthmta-3)/(t_(2)-t_(1))),'nodelay');
            BufferI2(:,sum(BufferI2==0)>0)=[];

        ymtaAr180b = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssumbyrbyar(BufferI2)),Time,'nearest'));
        
        
        %replace for plots and TonicREM/PhasicREM
        rmssum=interp1(t_,rmssum,Time,'nearest');
        rmssumbyr=interp1(t_,rmssumbyr,Time,'nearest');
        rmssumbyrbyar=interp1(t_,rmssumbyrbyar,Time,'nearest');
       
        refsig = log10(interp1(nanmean(t_(BufferI2)),nanmean(rmssum(BufferI2)),Time,'nearest'));
        
        p50 = prctile(refsig(~Exclude),50) * ones(length(LOC),1);
        
        EOGrmsConj_p50_180sb = ymtaAr180b - p50;

        p5 = prctile(refsig(~Exclude),5) * ones(length(LOC),1);
            p5_p50  = (p5 - p50);
        p95 = prctile(refsig(~Exclude),95) * ones(length(LOC),1);
            p95_p50 = (p95 - p50);
        p25 = prctile(refsig(~Exclude),25) * ones(length(LOC),1);
            p25_p50 = (p25 - p50);
        p75 = prctile(refsig(~Exclude),75) * ones(length(LOC),1);
            p75_p50 = (p75 - p50);
        p99 = prctile(refsig(~Exclude),99) * ones(length(LOC),1);
            p99_p50 = (p99 - p50);
        
        Tbl = table(EOGrmsConj_p50_180sb,p50,p5_p50,p25_p50,p75_p50,p95_p50,p99_p50);
        
        if exist('mdlREM')
            try
            PrREM = predict(mdlREM,Tbl);
            
            B = mdlREM.Coefficients.Estimate;
            Bintercept = B(1);
            
            Imain = find(strcmp('EOGrmsConj_p50_180sb',mdlREM.Coefficients.Properties.RowNames));
            
                
            Bmain = B(Imain);
            Rownames = mdlREM.Coefficients.Properties.RowNames;
            Rownames(Imain,:)=[]; B(Imain)=[];
            Rownames(1,:)=[]; B(1)=[];
            sumterms=0;
            Nskipterms = 0;
            for i=1:(length(B)-Nskipterms)
                sumterms = sumterms + eval(Rownames{i})*B(i);
            end
            Thres = p50 + -(Bintercept + sumterms)/Bmain;
            catch me
                PrREM = NaN*Time;
                Thres = NaN*Time;
             me.message
            end
        else
            PrREM = NaN*Time;
            Thres = NaN*Time;
        end


