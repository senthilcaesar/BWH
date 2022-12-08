function [ArThresData,Narousals,Nbreaths,ArThresDataAll,Iqualify] = ArThresNew(Ar_M,drive,drive_eupnea,win,settings,varargin)
    
    
    if ~isempty(varargin)
        if strcmp(varargin{1},'minwake')
            minwake=varargin{2};
            minwakearthres=varargin{3};
            clear minwakeB
            for i=1:length(win)
                if isnan(win(i))
                    minwakeB(i)=NaN;
                else
                    minwakeB(i)=minwake(win(i)+1);
                end
            end
            Ar_M(minwakeB>minwakearthres)=[];
            drive(minwakeB>minwakearthres)=[];
            win(minwakeB>minwakearthres)=[];
        end
    end

    %backupAr = Ar_M;
    BrY = settings(1); %ignore first X breaths after sleep onset
    swapbreathsforhigherdrive=settings(2);
    Nbreathsattributabletoarousal = settings(3);
    increasingdriveonly=settings(4);
    deleteifbeloweupnea=settings(5);
    setaseupneaifbeloweupnea=settings(6);
    
    %ArOnset = [[NaN diff(Ar_M)]]; %1 means first breath of arousal
    %ArOnset(ArOnset==-1)=0;
    
    ArOnsetNext = [[diff(Ar_M) NaN]]; %1 means next breath is arousal, -1 means next breath is sleep
    ArOnsetNext(Ar_M==1)=NaN;
    
%     if length(settings>8)&&settings(9)==1 %use drive AT arousal rather than BEFORE
%         ArOnsetNext=ArOnset; %overwrite
%     end %needs testing
    
    %handle discontinuous windows of data 
    Iendwin = [find(diff(win(:))>0)];
    Istartwin = [1;Iendwin(1:end-1)+1];
    for i=1:length(Istartwin)
        drive(Istartwin(i):Istartwin(i)+BrY-1)=NaN;
    end
    Ar_M(Istartwin)=NaN;
    ArOnsetNext(Iendwin)=NaN;
    
    %find length of time during sleep
    Iarousalonset=find(ArOnsetNext==1);
    for i=1:length(Iarousalonset)
        Imostdistantsleep(i)=1+find(Ar_M(1:Iarousalonset(i))~=0,1,'last');
    end
    ArThres_Dur = Iarousalonset-Imostdistantsleep+1;
    
    %remove short sleep periods
    Iarousalonset(ArThres_Dur<(BrY))=[];
    Imostdistantsleep(ArThres_Dur<(BrY))=[];
    ArThres_Dur(ArThres_Dur<(BrY))=[];
    
    %ArThres_Dur = ArThres_Dur-BrY;
    %I_ArThresL = I_ArThresL+BrY;
    
    I_ArThresInclBr=[];
    for i=1:length(ArThres_Dur)
        for j=1:(ArThres_Dur(i)-BrY)
            I_ArThresInclBr(end+1) = Iarousalonset(i)-ArThres_Dur(i)+BrY+j;
        end
    end
    
    
    
    %if usePmusinstead
        scores = drive(I_ArThresInclBr);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gradientdrive = gradient(drive);
        gradientdrive = gradientdrive(I_ArThresInclBr);%doesn't correct for a window change
%     elseif 1
%         scores = F_FlowtoEdiRatio_used(P)*Vdot_edi_M(I_ArThresInclBr)/Veupnea_M*60;
%         gradientdrive = gradient(F_FlowtoEdiRatio_used(P)*Vdot_edi_M/Veupnea_M*60); %doesn't correct for a window change
%     else
%         scores = Vchem_M(I_ArThresInclBr)'+1;
%         gradientdrive = gradient(Vchem_M+1); %doesn't correct for a window change
%     end
    
    labels = ArOnsetNext(I_ArThresInclBr); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    IqualifyA = zeros(length(ArOnsetNext),1);
    IqualifyA(I_ArThresInclBr)=1;
    
    Iqualify = 1*(ArOnsetNext(:)==1 & IqualifyA & ~isnan(drive(:)));
    %Iuse = find(Iqualify==1);
    %nanmedian(drive(Iuse))
    
    if swapbreathsforhigherdrive
        I = find(ArOnsetNext(I_ArThresInclBr)==1);
        for i=1:length(I)
            rangeleft = I(i)-Nbreathsattributabletoarousal+1;
            if rangeleft<1; rangeleft=1; end
            rangeright = I(i);
            range = rangeleft:rangeright;
            [tempval,tempI]=max(scores(rangeleft:rangeright));
            labels(I(i))=0;
            labels(rangeleft+tempI-1)=1;
        end
    end
    
    %default
    criteria = scores*0+1;    
    
    if increasingdriveonly&&deleteifbeloweupnea
        criteria = (scores>drive_eupnea)&(gradientdrive>0);
    elseif increasingdriveonly
        criteria = gradientdrive>0;
    elseif deleteifbeloweupnea
        criteria = scores>drive_eupnea;
    end
    
    scores(criteria==0)=[];
    labels(criteria==0)=[];
    
    Narousals = sum(labels==1);
    Nbreaths = length(labels);
    
    %diffa = [diff(nota);NaN];
    %prea = diffa==-1;
    %temp = [nota a diffa];
        %minarthres=1*x_;
        %arthresarray=x(prea==1&x>0.9*x_);
        %arthresarray(arthresarray<minarthres)=minarthres; %or = 
    if settings(7) 
        if setaseupneaifbeloweupnea
            scores(scores<drive_eupnea)=drive_eupnea;
        end
        ArThresData=nanmedian(scores(labels==1));
        ArThresDataAll=scores(labels==1&~isnan(scores));
    else
    posclass = 1;
    [X,Y,T,AUC_Arthres,OPTROCPT] = perfcurve(labels,scores,posclass); %need to find the threshold value that gives the OPTROCPT!
    %X= 1- specificity, OPTROCPT = FPR TPR
    [~,tempi] = max(-X + Y);
    if settings(8)==-1
        ArThresData = T(tempi);
    else
        ArThresData = T(find(X>settings(8),1));
    end
        ArThresDataAll=NaN;
    end
%     ArThres_ROC_(4,P) = T(find(X>0.1,1));
%     ArThres_ROC_(5,P) = T(find(X>0.05,1));
%     
%     xdata=T;
%     ydata=X;
%     lsqoptions=optimset('display','off');
%     x0 = [0.5 7];
%     upper = [10 50];
%     lower = [0.5 0.001];
%     fun = @(x,xdata)1./(1+exp(x(2)*(xdata-x(1))));
%     %fun = @(x,xdata)1./(1+exp(-x(2)*(xdata-x(1))));
%     [Xmodel]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper,lsqoptions);
%     xdatamodel = 0:0.001:9;
%     ydatamodel = fun(Xmodel,xdata);
%     ydatamodel2 = fun(Xmodel,xdatamodel);
%     findthres = [0.05 0.1 0.2 0.25 0.3 0.33 0.5];
%     for i=1:length(findthres)
%         ArThres_ROC_(i+5,P)=xdatamodel(find(ydatamodel2<=findthres(i),1));
%     end
%     
%     figure(434);
%     subplot(2,1,1); stairs(X,Y); ylabel('TPR'); xlabel('FPR');
%     subplot(2,1,2); stairs(T,X); hold('on'); stairs(T,Y,'r'); h=legend('FPR','TPR'); set(h,'EdgeColor',[1 1 1]);
%     plot([ArThres_ROC_(4,P) ArThres_ROC_(4,P)],[0 1],'r:');
%     plot([ArThresEdi_M(P) ArThresEdi_M(P)],[0 1],'k:');
%     ylabel('TPR/FPR'); xlabel('VdriveEdi, Feupnea');