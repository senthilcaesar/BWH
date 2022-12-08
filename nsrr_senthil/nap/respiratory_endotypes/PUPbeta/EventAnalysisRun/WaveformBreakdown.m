function [FsinEst,FexpEst,FsquEst,HarmonicPeaks,T1] = WaveformBreakdown(v,Time,Models,ploton,T1only)
    
    N=length(v);
    T0=Time(end)-Time(1);
    v = v - mean(v);
    w = hann(N); w = w/rms(w);
    v2=v.*w;
    
    if ploton    
        figure(1); clf(1);
        subplot(1,2,1);
        plot(Time,[v v2]);
    end
    %% analyze waveform v
    
    %padding improves frequency resolution for peak detection
    paddingratio=29; %oddnumberinteger, %19 works well
    z=zeros(round(N*(paddingratio-1)/2),1);
    X = fft([z;v2;z])/(N*paddingratio)*2; %faster without windowing
    Pxx = (conj(X).*X)*(paddingratio^2);
    
    %Phase = phase(X); %unused phase info; %SS removed this because phase
    %function is obsolete; try angle() in future if needed.
    
        df = 1/(T0*paddingratio);
        df_ = 1/(T0);
        Freq = 0:df:(N*paddingratio-1)*df;
        trange=[10 100];
        frange=fliplr(1./trange);
        frangei=round(frange./df)+1;
        IsearchmaxF=frangei(1):frangei(2);
        Harmonics = (1:8)'; %uses first H harmonics
    
    totalpower=sum(Pxx)*df;
    
    [Pxx1,f1_iround]=max(Pxx(IsearchmaxF));
    f1_iround = f1_iround + frangei(1)-1;
    FWaveEst=Freq(f1_iround);
    
    f1to5 = FWaveEst*Harmonics;
    T1 = 1/f1to5(1);
    if exist('T1only') && T1only==1
        FsinEst=[];
        FexpEst=[];
        FsquEst=[];
        HarmonicPeaks=[];
        return
    end
    
    f1to5_iround = round(f1to5/df) + 1;
    
    %harmonics are all normalzied to first
    HarmonicPeaks = Pxx(f1to5_iround)/Pxx(f1_iround); %shuold be same as: Pxx(f1to5_iround)/Pxx(f1to5_iround(1));
    HarmonicPeaksAbs = Pxx(f1to5_iround);
    H2 = HarmonicPeaks(2);
    H3 = HarmonicPeaks(3);
    H4 = HarmonicPeaks(4);
    H5 = HarmonicPeaks(5);
    H6 = HarmonicPeaks(6);
    
    if ~isempty(Models)
        T = HarmonicPeaks(2:6)'; %linear is better for squ and exp, less good for sine
        Tlog = log10(HarmonicPeaks(2:6)');
        T2=[T Tlog];
        VarNames = {'H2','H3','H4','H5','H6','logH2','logH3','logH4','logH5','logH6'};
        T2=array2table(T2); T2.Properties.VariableNames=VarNames;
        FsinEst = predict(Models.MdlSin,T2);
            FsinEst(FsinEst>1)=1; FsinEst(FsinEst<0)=0;
        FexpEst = predict(Models.MdlExp,T2);
            FexpEst(FexpEst>1)=1; FexpEst(FexpEst<0)=0;
        FsquEst = predict(Models.MdlSqu,T2);
            FsquEst(FsquEst>1)=1; FsquEst(FsquEst<0)=0;
    else
        FsinEst=NaN;
        FexpEst=NaN;
        FsquEst=NaN;
    end
    
    
    if ploton
        subplot(1,2,2);
        semilogy(Freq,Pxx);
        
        hold('on');
        I=[1:6];
        semilogy(Freq(f1to5_iround(I)),Pxx(f1to5_iround(I)),'b.','markersize',10);
        
        xlim([0 6.5*FWaveEst]);
        ylim([Pxx(f1to5_iround(1))/10000 1.1*Pxx(f1to5_iround(1))]);
    end        
