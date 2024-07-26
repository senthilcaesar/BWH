function [sig3]=FindWaveLetCoefs(sig,WaveFun,fs,mode)
%[D,t,sig2]=FindWaveLetCoefs(sig,WaveLetLevel,WaveFun,fs,plotflg,Flow)
%sig: input signal
%WaveLetLevel: wavelet decomposition level
%WaveFun: Wavelet basis function

NewFs=300;
NewSig=resample(sig,NewFs,fs);

                    
Time=(0:length(sig)-1)/fs;


if strcmp(mode,'Volume')
    WaveLetLevel=11; %%% To capture the frequency range between [0.1666 0.5]Hz, 
                    %this is the frequency range of normal breathing at rest  
    [C,L] = wavedec(NewSig,WaveLetLevel,WaveFun);
    t=(0:L(WaveLetLevel+2)-1)/NewFs;
    C2=zeros(size(C));
    nzInd=L(2)+1:sum(L(1:9));
    C2(nzInd)=C(nzInd);
    sig2 = waverec(C2,L,WaveFun);
    
elseif strcmp(mode,'FlowDC')
    WaveLetLevel=12; %%% To capture the frequency range between [0.1666 0.5]Hz, 
                    %this is the frequency range of normal breathing at rest  
    [C,L] = wavedec(NewSig,WaveLetLevel,WaveFun);
    t=(0:L(WaveLetLevel+2)-1)/NewFs;
    C2=zeros(size(C));
    nzInd=1:L(1);%L(1)+1:sum(L(1:12));
    C2(nzInd)=C(nzInd);
    sig2 = waverec(C2,L,WaveFun);
    
elseif strcmp(mode,'Flow HP')
    WaveLetLevel=10; %%% To capture the frequency range between [0.1666 0.5]Hz, 
                    %this is the frequency range of normal breathing at rest  
    [C,L] = wavedec(NewSig,WaveLetLevel,WaveFun);
    t=(0:L(WaveLetLevel+2)-1)/NewFs;
    C2=zeros(size(C));
    nzInd=[sum(L(3))+1:sum(L(1:6)) ];
    C2(nzInd)=C(nzInd);
    sig2 = waverec(C2,L,WaveFun);
elseif strcmp(mode,'Flow LP')
    WaveLetLevel=10; %%% To capture the frequency range between [0.1666 0.5]Hz, 
                    %this is the frequency range of normal breathing at rest  
    [C,L] = wavedec(NewSig,WaveLetLevel,WaveFun);
    t=(0:L(WaveLetLevel+2)-1)/NewFs;
    C2=zeros(size(C));
    nzInd=1:sum(L(1:7));
    C2(nzInd)=C(nzInd);
    sig2 = waverec(C2,L,WaveFun);
elseif strcmp(mode,'Flow LP2')
    WaveLetLevel=10; %%% To capture the frequency range between [0.1666 0.5]Hz, 
                    %this is the frequency range of normal breathing at rest  
    [C,L] = wavedec(NewSig,WaveLetLevel,WaveFun);
    t=(0:L(WaveLetLevel+2)-1)/NewFs;
    C2=zeros(size(C));
    nzInd=1:sum(L(1:4));
    C2(nzInd)=C(nzInd);
    sig2 = waverec(C2,L,WaveFun);
elseif strcmp(mode,'Flow Noise')
    WaveLetLevel=10; %%% To capture the frequency range between [0.1666 0.5]Hz, 
                    %this is the frequency range of normal breathing at rest  
    [C,L] = wavedec(NewSig,WaveLetLevel,WaveFun);
    t=(0:L(WaveLetLevel+2)-1)/NewFs;
    C2=zeros(size(C));
    nzInd=[sum(L(8))+1:sum(L(1:11)) ];
    C2(nzInd)=C(nzInd);
    sig2 = waverec(C2,L,WaveFun);
end


sig3=resample(sig2,fs,NewFs);
tt=(0:length(sig3)-1)/fs;

while Time(end)>tt(end)
    tt(end+1)=NaN;
    sig3(end+1)=NaN;
end

while Time(end)<tt(end)
    tt(end)=[];
    sig3(end)=[];
end

    

