function FlowFilterDetect = FlowFilterDetector(Flow,Time,ploton,verbose, plotdims)

%%
developmentmode=0;

if ~exist('plotdims')
    plotdims = [8 6];
end


if developmentmode
    %Commented out for safety. Saving new moedls should be intentional:
    %save('FilterDetectModels','mdlLPFadj','LPFtransmdls','HPFtransmdls','predrange','mdlHPF','predrangeH');
else
    load('FilterDetectModels','mdlLPFadj','LPFtransmdls','HPFtransmdls','predrange','mdlHPF','predrangeH');
end

%load data here

if 0 | ~developmentmode %assess raw data
    lowpassmethod=0; %set to 0
    LFinputtest=25; %set to 25, but not important
    HFinputtest=0; %set to 0
    HFmethod=0;
else %simulate problems
    lowpassmethod=2; %
    LFinputtest=10; %
    HFinputtest=0.03; %
    HFmethod=1; %
end

%Get fft data
[log10Prel,log10PrelH]=FlowFilterSimulator(Flow,Time,lowpassmethod,LFinputtest,HFmethod,HFinputtest,ploton, plotdims);

%Estimate LF
[temp,temp2] = predict(mdlLPFadj,log10Prel(:,[predrange]));
FlowFilterDetect.LowPassPredict = [10.^temp 10.^temp2];
FlowFilterDetect.LowPassPredict_1stOrderEquivalent = interp1(10.^predict(LPFtransmdls{1},log10([0.1:50]')),[0.1:50]',FlowFilterDetect.LowPassPredict);
FlowFilterDetect.LowPassPredict_2ndOrderEquivalent = interp1(10.^predict(LPFtransmdls{2},log10([0.1:50]')),[0.1:50]',FlowFilterDetect.LowPassPredict);
FlowFilterDetect.LowPassPredict_samplerateEquivalent = 2*interp1(10.^predict(LPFtransmdls{4},log10([0.1:50]')),[0.1:50]',FlowFilterDetect.LowPassPredict);

%Estimate HF
[temp,temp2] = predict(mdlHPF,log10PrelH(:,[predrangeH]));
FlowFilterDetect.HighPassPredict = [10.^temp 10.^temp2];
FlowFilterDetect.HighPassPredict_1stOrderEquivalent = interp1(10.^predict(HPFtransmdls{1},log10([0.00000001 0.001:0.001:1]')),[0.00000001 0.001:0.001:1]',FlowFilterDetect.HighPassPredict);
FlowFilterDetect.HighPassPredict_2ndOrderEquivalent = interp1(10.^predict(HPFtransmdls{2},log10([0.00000001 0.001:0.001:1]')),[0.00000001 0.001:0.001:1]',FlowFilterDetect.HighPassPredict);


%Disp to screen, debug
if 0
    disp(FlowFilterDetect.LowPassPredict_1stOrderEquivalent);
    disp(FlowFilterDetect.LowPassPredict_2ndOrderEquivalent);
    disp(FlowFilterDetect.LowPassPredict_samplerateEquivalent);
    disp(FlowFilterDetect.HighPassPredict);
end

if verbose
    
    %Disp to screen
    if FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1)<10
        disp('Warning: Flow appears smoothed or downsampled');
    else
        disp('Checked: Flow appears free of smoothing');
    end
    disp(['   2nd order cutoff equivalent ' num2str(FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1),2) ' Hz']);
    disp(['   1st order cutoff equivalent ' num2str(FlowFilterDetect.LowPassPredict_1stOrderEquivalent(1),2) ' Hz']);
    disp(['   Sampling rate equivalent ' num2str(FlowFilterDetect.LowPassPredict_samplerateEquivalent(1),2) ' Hz']);
    
    if FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1)>0.011
        disp('Warning: Flow appears distorted by baseline adjustment (high pass)');
    else
        disp('Checked: Flow appears to have a preserved DC baseline');
    end
    disp(['   2nd order cutoff equivalent  ' num2str(FlowFilterDetect.HighPassPredict_2ndOrderEquivalent(1),2) ' Hz']);
    disp(['   1st order cutoff equivalent  ' num2str(FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1),2) ' Hz']);
    
end