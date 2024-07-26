clear all
close all


%set ploton=1 to view waveforms with "pause"
%set ploton=0 to train model, takes several minutes
savenewworkspace=0; %default
ploton=0;


%% Setup
T=100;
Fs = 1;
dt = 1/Fs;
Time = (-T:dt:T-dt)';
N=length(Time);
T0 = N*dt;

%% Load classification models
if ploton
load('workspaceWF','MdlSin','MdlExp','MdlSqu');
Models.MdlSin=MdlSin;
Models.MdlExp=MdlExp;
Models.MdlSqu=MdlSqu;
else
    Models=[];
end

%% Filter

    filter_HFcutoff_butter0 = 0.05;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    
    
%% Run
clear AllDataTot
for xx=1:50000
    %% Generate waveform v
    disp(num2str(xx));
    
    wavetype = round(1 + rand(1)); %hybrid sin-square or sin-exp
    duty=0.2+0.3*rand(1);
%     alpha=rand(1);
    apneaat=-0.5-2.5*rand(1);
    
    Twave=30+60*rand(1);
    
    filter_HFcutoff_butter0 = 4/Twave;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    
    tau=Twave/(0.01 + 4.99*rand(1)); %<0.1 is sawtooth
%   
%     Fsine = 1-alpha;
%     
%     Fexp_=rand(1); %how exponential versus triangular
%     tau=Twave/(0.01 + 4.99*Fexp_); %<0.1 is sawtooth
%     
%     
%     Fsqu = alpha*(wavetype==1);
%     Fexp = alpha*(wavetype==2);

    
    Fsine_ = rand(1);
    Fexp_ = rand(1);
    Fsqu_ = rand(1);

    Fsine=Fsine_/sum([Fsine_,Fsqu_,Fexp_]);
    Fexp=Fexp_/sum([Fsine_,Fsqu_,Fexp_]);
    Fsqu=Fsqu_/sum([Fsine_,Fsqu_,Fexp_]);
        
    %noise = pinknoise(N);
    %noise=0.01*noise(:)/std(noise); %correction for precise power.    
    noise=0;
    
%     switch wavetype
%         case 9 %pure sin
%             v = sin(2*pi*Time/Twave);
%         case 10 %pure square
%             v = square(2*pi*Time/Twave,duty*100);
%         case 11 %pure exponential/sawtooth
%             v = -0.5+1/(1-exp(-Twave/tau))*(1*exp(-mod(Time,Twave)/tau) - exp(-Twave/tau)) ; % Get a vector.  No noise in this Y yet.
%         case 1 %hybrid square-sin
%             v1 = square(2*pi*Time/Twave,duty*100);
%             v2 = sin(2*pi*Time/Twave);
%             v1=v1/rms(v1);
%             v2=v2/rms(v2);
%             v = alpha*v1 + (1-alpha)*v2;
%         case 2 %hybrid exp-sin
            v1 = 0.5*sin(2*pi*Time/Twave);
            v2 = -0.5+1/(1-exp(-Twave/tau))*(1*exp(-mod(Time,Twave)/tau) - exp(-Twave/tau)) ; % Get a vector.  No noise in this Y yet.
            v3 = square(2*pi*Time/Twave,duty*100);
            v1=v1/rms(v1);
            v2=v2/rms(v2);
            v3=v3/rms(v3);
            v = Fsine*v1 + Fexp*v2 + Fsqu*v3;
%     end
    v=v/rms(v);
    v=v+noise;
    v(v<apneaat)=apneaat;
    
    v = filter(B_butter0,A_butter0,v);
    
    %artificially add a pre-window
    %N2 = round(N*((0.5+0.5*rand(1)))/2)*2;
    %wA = [zeros((N-N2)/2,1);hann(N2);zeros((N-N2)/2,1)];
    %v=v.*wA;
    
    %% window data for analysis
%     [FsinEst,FexpEst,FsquEst,HarmonicPeaks,T1] = WaveformBreakdown(v,Time,Models,ploton);
%     
%     [FsinEst,FexpEst,FsquEst,HarmonicPeaks] = WaveformBreakdown(v,Time,Models,ploton);
%     
    [~,~,~,~,T1] = WaveformBreakdown(v,Time,Models,ploton);
    I = find(Time>-1.67*T1&Time<1.67*T1);
    [FsinEst,FexpEst,FsquEst,HarmonicPeaks,~] = WaveformBreakdown(v(I),Time(I),Models,ploton);


    %% analyze waveform v
    
    if ploton
        ShapeEst = table(FsinEst,FexpEst,FsquEst)
        if 1
            ShapeAct = table(Fsine,Fexp,Fsqu)
        end
    end
    
    ShapeTable = table(Fsine,Fexp,Fsqu);
    H2=HarmonicPeaks(2);
    H3=HarmonicPeaks(3);
    H4=HarmonicPeaks(4);
    H5=HarmonicPeaks(5);
    H6=HarmonicPeaks(6);
    Predictors = table(H2,H3,H4,H5,H6);
    AllData = [ShapeTable Predictors];
    if exist('AllDataTot')
        AllDataTot = [AllDataTot;AllData];
    else
        AllDataTot = AllData;
    end
    
    if ploton
        pause
    end
end

%% Make model for later use

T = (AllDataTot{:,[4:8]}); %linear is better for squ and exp, less good for sine
Tlog = log10(AllDataTot{:,[4:8]});
T2=[T Tlog];
VarNames = {'H2','H3','H4','H5','H6','logH2','logH3','logH4','logH5','logH6'};
methodA='quadratic'; %'quadratic' 'constant'
methodX='quadratic'; %'quadratic'
methodY='normal'; %binomial was less effective

T2=array2table(T2); T2.Properties.VariableNames=VarNames;

%Models are trained to detech how "sinusoidal", "exponential" or "rectangular" a pattern is 
T2.y = AllDataTot.Fsine; %H5 low, H5on2 high
MdlSin=stepwiseglm(T2,methodA,'upper',methodX,'distribution',methodY)

T2.y = AllDataTot.Fexp; %H5 low, H5on2 high
MdlExp=stepwiseglm(T2,methodA,'upper',methodX,'distribution',methodY)

T2.y = AllDataTot.Fsqu; %H5 low, H5on2 high
MdlSqu=stepwiseglm(T2,methodA,'upper',methodX,'distribution',methodY)

disp('model fits are complete');

% check Rsquared
if exist('Rsqu')
Rsqbackup=Rsq;
end

% 
clear Rsq
FsinEst = predict(MdlSin,T2);
Rsq.FsinEst=1-MdlSin.SSE/MdlSin.SST;
FexpEst = predict(MdlExp,T2);
Rsq.FexpEst=1-MdlExp.SSE/MdlExp.SST;
FsquEst = predict(MdlSqu,T2);
Rsq.FsquEst=1-MdlSqu.SSE/MdlSqu.SST;


%% save

if 0
    save('workspaceWF2','T2','AllDataTot','MdlSin','MdlExp','MdlSqu');
end


%%
