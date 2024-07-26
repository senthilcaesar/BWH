if ~exist('mdlSmooth')
    save FilterSmoothCutoffDetector
end

%%

for Exp=31:33%26:30%1:20
    switch Exp
        case 1
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5229_XHz.mat'); %RICCADSA
        case 2
            load('5032_XHz.mat');
        case 3
            load('5012_XHz.mat');
        case 4
            load('5447_XHz.mat');
        case 5
            load('5478_XHz.mat');
        case 6
            cd('G:\Dropbox (Personal)\SaraODB OAT\Converted2019');
            load('UZA012_XHz.mat') %Antwerp OAT
        case 7
            load('UZA024_XHz.mat')
        case 8
            load('UZA039_XHz.mat')
        case 9
            load('UZA067_XHz.mat')
        case 10
            load('UZA085_XHz.mat')
        case 11
            cd('G:\Dropbox (Partners HealthCare)\MAD-OX\Traits\Converted');
            load('1115N0_XHz.mat') %MADOX, pneumotach
        case 12
            load('1731N0_XHz.mat') %MADOX, pneumotach
        case 13
            load('1870N0_XHz.mat') %MADOX, pneumotach
        case 14
            load('1881N0_XHz.mat') %MADOX, pneumotach
        case 15
            load('1884N0_XHz.mat') %MADOX, pneumotach
        case 16
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5111_XHz.mat'); %RICCADSA
        case 17
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5112_XHz.mat'); %RICCADSA
        case 18
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5113_XHz.mat'); %RICCADSA
        case 19
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5115_XHz.mat'); %RICCADSA
        case 20
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5117_XHz.mat'); %RICCADSA
        case 21
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5119_XHz.mat'); %RICCADSA
        case 22
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5120_XHz.mat'); %RICCADSA
        case 23
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5121_XHz.mat'); %RICCADSA
        case 24
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5122_XHz.mat'); %RICCADSA
        case 25
            cd('G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Converted');
            load('5124_XHz.mat'); %RICCADSA
        case 26
            cd('G:\Dropbox (Personal)\SaraODB OAT\Converted2019');
            load('UZA002_XHz.mat') %Antwerp OAT
        case 27
            load('UZA004_XHz.mat')
        case 28
            load('UZA005_XHz.mat')
        case 29
            load('UZA010_XHz.mat')
        case 30
            load('UZA013_XHz.mat')
            
            % DLM data, case 31-33
        case 31
            DataDir = ['C:\Users\dwayne\OneDrive\OneDrive - The University of Queensland\Projects\NSRR\MESA\SampleData\Converted\'];
            load([DataDir, '3010023_20110817_XHz.mat']);
        case 32
            DataDir = ['C:\Users\dwayne\OneDrive\OneDrive - The University of Queensland\Projects\RainePilot\Converted\'];
            load([DataDir, '20170330_402_flowfiltersoff_XHz.mat']);
        case 33
            DataDir = ['C:\Users\dwayne\OneDrive\OneDrive - The University of Queensland\Projects\GMRF\Converted\Filters_allOff\'];
            load([DataDir, 'P005_XHz.mat']);
    end
    
    %% filtering
    
    flow = DataEventHypnog_Mat(:,2);
    dt = DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1);
    
    Ierr = isnan(flow)| flow==Inf | flow==-Inf;
    
    flow(Ierr)=nanmedian(flow);
    
    freqinput = [1:12 15 20];
    %clear log10Prel
    freqsweep = [1:12];
    for j=1:length(freqinput)
        j
        filter_HFcutoff_butter1 = freqinput(j);
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter1]/(1/dt/2),'low');
        flowfiltered = filter(B_butter0,A_butter0,flow);
        
        %flowfiltered(Ierr)=0;
        %filter_order0 = 1;
        
        nfft = 200*round(1/dt);
        [Flowpsd,F] = pwelch(flowfiltered-mean(flowfiltered),hann(nfft),round(nfft/2),nfft,1/dt);
        df = F(2)-F(1);
        Fmax = 20;
        I = F>Fmax;
        Flowpsd(I)=[];
        F(I)=[];
        if 0
            figure(1); clf(1);
            semilogy(F,Flowpsd)
        end
        Fl = 0.1; Fr = 12.5; Fli = round(Fl/df) + 1; Fri = round(Fr/df) ; %not inclusive
        P1to12 = sum(Flowpsd(Fli:Fri))*df;
        
        range_ = [1:12]';
        clear P
        for i=1:size(range_,1)
            Fl = range_(i)-0.5; Fr = range_(i)+0.5; Fli = round(Fl/df) + 1; Fri = round(Fr/df) ; %not inclusive
            Prel(i,1) = sum(Flowpsd(Fli:Fri))*df/P1to12;
        end
        log10Prel(j,:) = log10(Prel)';
        
        %
        % SDratiodiagnosis = interp1(freqsweep,log10Prel,filter_HFcutoff_butter1)
    end
    
    %%
    if ~exist('log10Prel_')
        log10Prel_ = log10Prel;
        freqinput_ = freqinput;
    else
        log10Prel_ = [log10Prel_;log10Prel];
        freqinput_ = [freqinput_ freqinput];
    end
    
    %%
    predrange = [1:6]
    mdlSmooth = fitglm(log10Prel_(:,[predrange]),freqinput_)
    
    % mdlSmooth = stepwiseglm(log10Prel_,freqinput_,'linear','PEnter',0.5,'PRemove',0.51,'Criterion','sse')
    %
    % %mdlSmooth = stepwiseglm(log10Prel_,freqinput_,'interactions','PEnter',0.005,'PRemove',0.010,'Criterion','sse')
    %
    % NtermsCurrent = length(mdlSmooth.Coefficients.Estimate)-1
    % Nsteps=NtermsCurrent-25%NtermsCurrent-50
    % pthres = 0.000000000000000000000000000000000000000000001
    % mdlSmooth = step(mdlSmooth,'Criterion','sse','PEnter',pthres,'PRemove',pthres*1.01,'Nsteps',Nsteps,'upper','linear')
    % NtermsCurrent = length(mdlSmooth.Coefficients.Estimate)-1
    
    %%
    [FSmoothPred]=predict(mdlSmooth,log10Prel_(:,[predrange]));
    figure(1);
    plot(FSmoothPred,freqinput_,'.');
    if savefigs
        fig = gcf;
        str = ['FSmoothPredVsFreqInput'];
        print(fig, str, '-dtiff', '-r300');
    end
    %% Run New data
    
    
    %load the _XHz data then run this cell
    flow = DataEventHypnog_Mat(:,2);
    dt = DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1);
    Ierr = isnan(flow)| flow==Inf | flow==-Inf;
    flow(Ierr)=nanmedian(flow);
    
    %or skip the above and provide "flow" and "dt" then proceed below
    %note: need to load "FilterSmoothCutoffDetector" to get mdlSmooth
    
    if 0 %simulate a filter
        freqinputtest = 3;
        %clear log10Prel
        filter_HFcutoff_butter1 = freqinputtest;
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter1]/(1/dt/2),'low');
        flowfiltered = filter(B_butter0,A_butter0,flow);
    else %test a signal
        flowfiltered = flow;
    end
    
    
    freqsweep = [1:12];
    freqsweep = freqsweep(predrange);
    nfft = 200*round(1/dt); %200 seconds per fft
    [Flowpsd,F] = pwelch(flowfiltered-mean(flowfiltered),hann(nfft),round(nfft/2),nfft,1/dt);
    df = F(2)-F(1);
    Fmax = 20;
    I = F>Fmax;
    Flowpsd(I)=[];
    F(I)=[];
    if 1
        figure(10); clf(10);
        loglog(F,Flowpsd);
    end
    Fl = 0.1; Fr = 12.5; Fli = round(Fl/df) + 1; Fri = round(Fr/df) ; %not inclusive
    P1to12 = sum(Flowpsd(Fli:Fri))*df;
    
    range_ = freqsweep';
    clear Prel log10PrelTest
    for i=1:size(range_,1)
        Fl = range_(i)-0.5; Fr = range_(i)+0.5; Fli = round(Fl/df) + 1; Fri = round(Fr/df) ; %not inclusive
        Prel(i,1) = sum(Flowpsd(Fli:Fri))*df/P1to12;
    end
    log10PrelTest = log10(Prel)';
    
    %
    % SDratiodiagnosis = interp1(freqsweep,log10Prel,filter_HFcutoff_butter1)
    
    [FSmoothEst,FSmoothEst95CI]=predict(mdlSmooth,log10PrelTest)
end


%%
if 0
    save FilterSmoothCutoffDetector mdlSmooth log10Prel_ freqinput_ -v7.3 
end
