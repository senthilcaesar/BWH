function [ARpred,WASPr]=RunWASar2(DataEventHypnog_Mat_ds,mdlAR,dT)
    
    logitinverse = @(p) 1./(1+exp(-p));
    
    N = height(DataEventHypnog_Mat_ds);
    %DataEventHypnog_Mat_ds = table(ARieF_pred_logit);
    Time = [0:dT:((N-1)*dT)]';

    dTtrain=3; %must stay 3
    mtalength = 1; %seconds
    makelist = mdlAR.Coefficients.Properties.RowNames(2:end-2);

    %make Xs (see dTtrain) moving time average to use for reference inputs
    mtalengthi = round(mtalength/dT);
    wts = [0.5;repmat(1,mtalengthi-1,1);0.5];
    wts=wts/sum(wts);
    yS = conv(DataEventHypnog_Mat_ds.WakeIntensity,wts);
    yS(1:round(mtalengthi/2))=[];
    yS(length(DataEventHypnog_Mat_ds.WakeIntensity)+1:end)=[];
    WakeIntensity_mta = yS;
    
    yS = conv(DataEventHypnog_Mat_ds.SleepIntensity,wts);
    yS(1:round(mtalengthi/2))=[];
    yS(length(DataEventHypnog_Mat_ds.SleepIntensity)+1:end)=[];
    SleepIntensity_mta = yS;
    clear yS

    for k=1:length(makelist)
        %makelist{k}
        if strcmp(makelist{k}(1),'L')
            shiftdir=-1;
        else %case 'N'
            shiftdir=+1;
        end
        shiftval = dTtrain*str2num(makelist{k}(2));
        
        if strcmp(makelist{k}(3),'w')
        %using mta for reference
        temp = interp1(Time,WakeIntensity_mta,Time + shiftdir*shiftval,'nearest');
        eval(['DataEventHypnog_Mat_ds.' makelist{k} '=temp;']);
        elseif strcmp(makelist{k}(3),'s')
        temp = interp1(Time,SleepIntensity_mta,Time + shiftdir*shiftval,'nearest');
        eval(['DataEventHypnog_Mat_ds.' makelist{k} '=temp;']);
        end
    end

    ARpred = predict(mdlAR,DataEventHypnog_Mat_ds);
    WASPr = max([ARpred logitinverse(DataEventHypnog_Mat_ds.WSBalance)]')';
    
    
