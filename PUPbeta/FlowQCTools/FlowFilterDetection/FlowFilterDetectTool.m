clear all; close all; clc;
%First run start here to set directories
disp('Select Folder with Converted Data');
mainFolder = uigetdir(); 
%mainFolder ='F:\numom2b\Converted';  % Use absolute paths
%addpath(genpath(pwd));
Flpattern=fullfile(mainFolder, '*.mat');
dinfo = dir(Flpattern); % mat files

%% common to all files
ploton=1;
verbose=1;
if verbose
    disp('    ------Flow quality analysis------    ');
end
ClipThresholdFmax=0.90;
ClipFthresholdFnormal=0.002; %higher value removes more (i.e. false) clipping (0.002)

%% Loop for each file
clear Subject
for jj=1:size(dinfo,1)
    Subject{jj,1}=dinfo(jj).name(1:end-8);
end
Subject = table(Subject);

for jj=1:size(dinfo,1)
    if verbose
        jj
        disp(dinfo(jj).name)
        
    end
    labfile=[dinfo(jj).folder '\' dinfo(jj).name];
    temp = load(labfile);
    if 1
        Flow =temp.DataEventHypnog_Mat(:,2);
        Time = temp.DataEventHypnog_Mat(:,1);
    else
        Flow = Pnasal;
        Time = TimeDN-TimeDN(1);
    end
    
    FlowFilterDetect = FlowFilterDetector(Flow,Time,ploton,verbose);
    Ttemp=array2table([FlowFilterDetect.LowPassPredict FlowFilterDetect.LowPassPredict_1stOrderEquivalent...
        FlowFilterDetect.LowPassPredict_2ndOrderEquivalent FlowFilterDetect.LowPassPredict_samplerateEquivalent...
        FlowFilterDetect.HighPassPredict FlowFilterDetect.HighPassPredict_1stOrderEquivalent  ...
        FlowFilterDetect.HighPassPredict_2ndOrderEquivalent]);
    
    
    % add clip detection:
    [~,~,FclippedI,FclippedE,~] = ClipDetection(Flow,[],[],[],ClipThresholdFmax,ClipFthresholdFnormal,1);
    
    Ttemp.FclippedI=FclippedI;
    Ttemp.FclippedE=FclippedE;
    Ttemp.FclippedTotal=FclippedE+FclippedI;
    
    if verbose
        
        if Ttemp.FclippedTotal(1)>0.005
            disp('Warning: Flow appears clipped');
        else
            disp('Checked: Flow appears free of clipping');
        end
        disp(['   Clipping fraction: ' num2str(100*Ttemp.FclippedTotal,2) ' %']);
    end
    
    if exist ('T')
        T=[T;Ttemp];
    else
        T=Ttemp;
       
    end
    
end

 T.Properties.VariableNames={'LowPassPredict','LowPassPredictL95','LowPassPredictU95',...
            'LowPassPredict_1stOrderEquivalent','LowPassPredict_1stOrderEquivalentL95','LowPassPredict_1stOrderEquivalentU95',...
            'LowPassPredict_2ndOrderEquivalent','LowPassPredict_2ndOrderEquivalentL95','LowPassPredict_2ndOrderEquivalentU95',...
            'LowPassPredict_samplerateEquivalent','LowPassPredict_samplerateEquivalentL95','LowPassPredict_samplerateEquivalentU95',...
            'HighPassPredict','HighPassPredictL95','HighPassPredictU95',...
            'HighPassPredict_1stOrderEquivalent','HighPassPredict_1stOrderEquivalentL95','HighPassPredict_1stOrderEquivalentU95',...
            'HighPassPredict_2ndOrderEquivalent','HighPassPredict_2ndOrderEquivalentL95','HighPassPredict_2ndOrderEquivalentU95',...
            'FclippedI','FclippedE','FclippedTotal'};

T = [Subject T];

        %% Plot
        
        set(groot,'defaultAxesfontname','arial narrow');
        set(groot,'defaultAxesfontsize',13);
        set(groot,'defaultFigureColor',[1 1 1]);
        set(groot,'defaultAxesBox','off');

        N=1:1:(size(T,1));
        figure(11); clf(11); set(gcf,'color',[1 1 1]);
        subplot(1,3,1);
        h=fill([1 size(T,1) size(T,1) 1],[20 20 50 50],[0.7 1 0.7],'edgealpha',0);
        hold on
        Y = T.LowPassPredict_samplerateEquivalent;
        Y(Y>50)=50;
        h=scatter(N,Y,7,[0 0 0],'filled','markerfacealpha',0.3); 
        ylabel('Sampling Rate Equivalent (Smoothing), Hz');
        %set(gca,'YScale','log')
        box off;
        xlim([min(N) max(N)*1.1]);
        
        plotbarwitherrorsmedian(nanmedian(Y),max(N)*1.05,prctile(Y,75),prctile(Y,25),[0 0 0],NaN,max(N)*0.05,50/200)



        subplot(1,3,2);
        ylow=0.001;
        yhigh=0.4;
        h=fill([1 size(T,1) size(T,1) 1],[ylow ylow 0.0128 0.0128],[0.7 1 0.7],'edgealpha',0);
        hold on
        Y = T.HighPassPredict_1stOrderEquivalent;
        Y(Y>yhigh)=yhigh;
        Y(Y<ylow)=ylow;
        scatter(N,Y,7,[0 0 0],'filled','markerfacealpha',0.3); hold on
       
        ylabel('Baseline distortion (HighPass), Hz');
        set(gca,'YScale','log')
        box off;
        xlim([min(N) max(N)*1.1]);
        ylim([ylow yhigh]);
        
        plotbarwitherrorsmedian(nanmedian(Y),max(N)*1.05,prctile(Y,75),prctile(Y,25),[0 0 0],NaN,max(N)*0.05,yhigh/400)
        set(gca,'ytick',[0.001 0.003 0.01 0.03 0.1 0.3],'yticklabels',{'0.001','0.003','0.01','0.03','0.1','0.3'});
        
        subplot(1,3,3);
        h=fill([1 size(T,1) size(T,1) 1],[0 0 0.5 0.5],[0.7 1 0.7],'edgealpha',0);
        hold on
        Y = T.FclippedTotal*100;
        Y(Y>5)=5;
        scatter(N,Y,7,[0 0 0],'filled','markerfacealpha',0.3); hold on
       
        ylabel('Clipping, %');
        %set(gca,'YScale','log')
        box off;
        xlim([min(N) max(N)*1.1]);
        ylim([0 5]);
        plotbarwitherrorsmedian(nanmedian(Y),max(N)*1.05,prctile(Y,75),prctile(Y,25),[0 0 0],NaN,max(N)*0.05,5/200)
        
        set(gcf,'position',[ 378   295   960   573])

        %% save table and figure
       
        
        
 idcs   = strfind(mainFolder,'\');
 newdir = mainFolder(1:idcs(end)-1);
 if ~(exist([newdir '\QC'], 'dir') == 7)
        mkdir([newdir '\QC']);
 end
    
save([newdir '\QC\' 'FlowFilterDetection.mat'],'T');
saveas(11,[newdir '\QC' '\FlowFilterDetection'],'fig');
saveas(11,[newdir '\QC' '\FlowFilterDetection' ],'png');
        
        
        
        