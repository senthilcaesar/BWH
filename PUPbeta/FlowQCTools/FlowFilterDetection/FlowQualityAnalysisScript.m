%% How to use model %%

%save('FilterDetectModels','mdlLPFadj','LPFtransmdls','HPFtransmdls','predrange','mdlHPF','predrangeH');
if 1
    Flow = DataEventHypnog_Mat(:,2);
    Time = DataEventHypnog_Mat(:,1);
else
    Flow = Pnasal;
    Time = TimeDN-TimeDN(1);
end

%copy this code into somewhere, Flow and Time must exist, output is FlowFilterDetect "structure"

ploton=1;
verbose=1;
if verbose
   disp('    ------Flow quality analysis------    ');
end
FlowFilterDetect = FlowFilterDetector(Flow,Time,ploton,verbose);

% add clip detection:
ClipThresholdFmax=0.90;
ClipFthresholdFnormal=0.002; %higher value removes more (i.e. false) clipping (0.002)
[~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[],ClipThresholdFmax,ClipFthresholdFnormal,1);

FlowFilterDetect.FclippedI=FclippedI;
FlowFilterDetect.FclippedE=FclippedE;
FlowFilterDetect.FclippedTotal=FclippedE+FclippedI;

if verbose
if FlowFilterDetect.FclippedTotal(1)>0.005
    disp('Warning: Flow appears clipped');
else
    disp('Checked: Flow appears free of clipping');
end
disp(['   Clipping fraction: ' num2str(100*FlowFilterDetect.FclippedTotal,2) ' %']);
end

if 0
    fig = gcf;
    saveas(fig, 'FrequencyAnalysis', 'png');
end
