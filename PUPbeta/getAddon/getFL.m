function [FLFinalT,HistDataT]=getFL(MrangeOverride,skipIEratio)

if ~exist('skipIEratio')
    skipIEratio=1;
end
% RUN StartHere.m first

% this function runs Flow limitation on already analyzed files
% creates Flow limitation metric table for the entire patient cohort.

global settings AMasterSpreadsheet


t_start = clock;

%% Load AMasterSpreadsheet

[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');

[num,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,1),'UniformOutput',false)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,1:2);

NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));

settings = ImportSettings(settings,AMasterSpreadsheet);


MaxM = size(num,1);

try
    M = max(settings.Mrange);
catch
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

if exist('MrangeOverride')
    Mrange=MrangeOverride;
else
    Mrange=settings.Mrange;
end
% intialize variables to zeros/NaN incase FL computation failed or no
% analysis file present
success = zeros(M,1);

%% load flow drive models
load('FinalModel_Table.mat');
load('FinalModel_Table10.mat');
load('FLcertModel.mat');
FlowDriveModels.FinalModelTable = FinalModelTable;
FlowDriveModels.FinalModelTable10 = FinalModelTable10;
FlowDriveModels.FLcert = FLcert;

%% Load analyzed data
for i=Mrange
    try
        filedir = [path{:} settings.savename '_' num2str(i) '.mat'];
        
        if exist(filedir)==2
            disp(['Processing: ' num2str(i) ': ' settings.savename '_' num2str(i) '.mat']);
            load(filedir,'BreathDataTable','BreathFLDataTable');
            
            %% check for flow quality -- calculate SNR and IE ratio
            if size(BreathDataTable,1)==1 && size(BreathDataTable,2)==1
                BreathDataTable2=BreathDataTable{1,1};
            else
                BreathDataTable2=BreathDataTable;
            end
                  
            % calculate SNR
            for jj=1:length(BreathDataTable2)
                clear BreathDataTableWin FlowWin TimeWin
                BreathDataTableWin=BreathDataTable2{1,jj};
                try
                    SNRwindow(jj,1)=nanmedian(BreathDataTableWin.SNRwindow);
                catch me
                    SNRwindow(jj,1)=NaN;
                end
            end
            
            % check if IE ratio is already present or not
            temp2=BreathDataTable2;
            
            clear IEratioW IEratioW2 IEratioEstimatedW IEratioEstimatedW2
            
            IEratiofound=zeros(size(temp2,2),1);
            for ii=1:length(temp2) % remove any empty cells in BreathDataTable
                if isempty(temp2{1,ii})
                    temp2{1,ii}=NaN;
                end
                try
                    IEratiofound(ii)=sum(ismember(temp2{1,ii}.Properties.VariableNames,'IEratio')); %check if IEratio exist.
                end
            end
            BreathDataTable2=temp2;
            
            % calculate IE ratio (untested code)
            if sum(IEratiofound)>0
            for jj=1:length(BreathDataTable2)
                clear BreathDataTableWin
                BreathDataTableWin=BreathDataTable2{1,jj};
                try
                    IEratioW(jj,1)=nanmean(BreathDataTableWin.IEratio);
                    IEratioW2(jj,1)=nanmean(10.^abs(log10(BreathDataTableWin.IEratio)));
                    IEratioEstimatedW(jj,1)=nanmean(BreathDataTableWin.IEratioEstimated);
                    IEratioEstimatedW2(jj,1)=nanmean(10.^abs(log10(BreathDataTableWin.IEratioEstimated)));
                catch me
                    IEratioW(jj,1)=NaN;
                    IEratioW2(jj,1)=NaN;
                    IEratioEstimatedW(jj,1)=NaN;
                    IEratioEstimatedW2(jj,1)=NaN;
                end
            end
            elseif sum(IEratiofound)==0 && ~skipIEratio
                warnText =sprintf('IE ratio not found\n Please run getIERatio hack tool to get IE variables');
                warning(warnText);
                IEratioW(jj,1)=NaN;
                IEratioW2(jj,1)=NaN;
                IEratioEstimatedW(jj,1)=NaN;
                IEratioEstimatedW2(jj,1)=NaN;
                
            end
            
            
            %% Compute Flow Drive
            
            % Load VE data from Analyzed Tables and generate NonOvlapped tables
            clear FlowDrive BreathDataTableFulls FlowLimCertOut FLdata HistData DeltaTime
            
            %legacy fix code:
            if iscell(BreathDataTable) && numel(BreathDataTable)==1
                BreathDataTable = BreathDataTable{1};
            end
            if iscell(BreathFLDataTable) && numel(BreathFLDataTable)==1
                BreathFLDataTable = BreathFLDataTable{1};
            end
            
            [~,~,~,~,BreathDataTableFulls]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
            
            % test BreathDataTable2 and BreathFLDataTable2, if they aren't tables
            % and if the size is [1 1] and the value isnan, then return
            if 0
                if ~isa(BreathDataTable2,'table') && all(size(BreathDataTable2)==[1 1]) && isnan(BreathDataTable2) && ...
                        ~isa(BreathFLDataTable2,'table') && all(size(BreathFLDataTable2)==[1 1]) && isnan(BreathFLDataTable2)
                    disp('Missing BreathData and BreathFLData');
                end
            end
            
            % Flow Drive
            [FlowDrive,~] = computeFlowDrive(FlowDriveModels.FinalModelTable,BreathDataTableFulls);
            [FlowLimCertOut,PrFL] = computeFLcertainty(FlowDriveModels.FLcert,BreathDataTableFulls); %not latest model
            
            FlowDrive(FlowDrive>1)=1;
            FlowDrive(FlowDrive<0.1)=0.1;
            FlowDrive(BreathDataTableFulls.VI<0.1)=0.1;
            FlowDrive(BreathDataTableFulls.Etype==2)=0.1; %apneas are v.severe FL
            FlowDrive(BreathDataTableFulls.Etype==3)=NaN; %central apneas are excluded
            %FlowDrive(BreathDataTableFulls.SNRwindow<10)=NaN; %remove low SNR
            %FlowDrive(abs(BreathDataTableFulls.IEratioActual)>1.9)=NaN; %remove breaths without either insp or exp signals.
            
            FlowDrive2=FlowDrive;
            FlowDrive(BreathDataTableFulls.hypnog_B==4)=NaN;
            FlowDrive2(BreathDataTableFulls.hypnog_B~=4)=NaN;
            
            FlowLimCertOut = double(FlowLimCertOut);
            FlowLimCertOut(BreathDataTableFulls.VI<0.1)=1; %very low flow assumed obstructive and called certain FL
            FlowLimCertOut(BreathDataTableFulls.Etype==2)=1; %apneas are certain FL
            FlowLimCertOut(BreathDataTableFulls.hypnog_B==4) = NaN;
            FlowLimCertOut(BreathDataTableFulls.Etype==3)=NaN; %central apneas are excluded
            
            FLdata.FLCfrequency = sum(FlowLimCertOut==1)./sum(FlowLimCertOut>0)*100;
            FLdata.FLfrequency = sum(FlowLimCertOut<3)./sum(FlowLimCertOut>0)*100;
            FLdata.FlowDriveMedian = nanmedian(FlowDrive*100);
            FLdata.FlowDriveUnder50p = sum(FlowDrive<0.5)./sum(FlowDrive>0)*100;
            FLdata.FlowDriveUnder70p = sum(FlowDrive<0.7)./sum(FlowDrive>0)*100;
            if 1
                FLdata.FlowDriveMedianW = nanmedian(FlowDrive2*100); %low value might mean noisy data
                FLdata.WakeBreathsN=nansum(BreathDataTableFulls.hypnog_B==4);
                DeltaTime=(BreathDataTableFulls.Time_end-BreathDataTableFulls.Time_start);
                FLdata.WakeBreathsTime=nansum(DeltaTime(BreathDataTableFulls.hypnog_B==4));
            end
            
            %% Plot
            figure(1111); clf(1111); set(gcf,'color',[1 1 1]);
            % blue-apnea; orange- controls
            color1f=[.3 .3 .3];
            color2f=[.85 .45 .05];
            color3f=[.85 .75 .05];
            dStep=2;
            Centers=10+dStep/2:dStep:100-dStep/2;
            Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2); Edges(end)=Inf; Edges(1)=-Inf;
            PlotData = FlowDrive*100;
            clear HistDataA
            PlotCrit = (FlowLimCertOut)==1; %certain FL
            [HistDataA(1,:),edges] = histcounts(PlotData(PlotCrit),Edges);
            PlotCrit = double(FlowLimCertOut)==2;
            [HistDataA(2,:),edges] = histcounts(PlotData(PlotCrit),Edges);
            PlotCrit = double(FlowLimCertOut)==3; %normal
            [HistDataA(3,:),edges] = histcounts(PlotData(PlotCrit),Edges);
            
            C = HistDataA(1,:);
            P = HistDataA(2,:);
            N = HistDataA(3,:);
            HistData = [array2table(C) array2table(P) array2table(N)];
            
            %HistDataA1 = reshape(HistDataT{:,:},[50 3])';  %code to get back data for hist plots for a single subject
            
            NbreathsUsed = sum(sum(HistDataA))/100;
            h1 = HistDataA(1,:)/NbreathsUsed;
            h2 = HistDataA(2,:)/NbreathsUsed;
            h3 = HistDataA(3,:)/NbreathsUsed;
            
            bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',1,'BarWidth',1,'FaceColor',color1f,'facealpha',0.7);
            hold on
            bar(Centers,h2,'EdgeAlpha',0,'FaceAlpha',1,'BarWidth',1,'FaceColor',color2f,'facealpha',0.7);
            bar(Centers,h3,'EdgeAlpha',0,'FaceAlpha',1,'BarWidth',1,'FaceColor',color3f,'facealpha',0.7);
            hold off
            
            set(gca,'box','off','tickdir','out','fontname','arial narrow');
            xlim([10 100]);
            
            %% get summary values for each subject
            try
                FLdata.IEratio = nanmedian(IEratioW);
                FLdata.IEratio2 = nanmedian(IEratioW2);
                FLdata.IEratioEstimated = nanmedian(IEratioEstimatedW);
                FLdata.IEratioEstimated2 = nanmedian(IEratioEstimatedW2);
            end
            try
                FLdata.SNRwindow = nanmedian(SNRwindow);
                FLdata.SNRwindowNoiseMean = 10*log10(1./nanmean(1./(10.^([SNRwindow]/10)))); %odd numbers
                FLdata.FSNRover15 = nansum(SNRwindow>15)/nansum(SNRwindow>-Inf);
                FLdata.FSNRover20 = nansum(SNRwindow>20)/nansum(SNRwindow>-Inf);
                FLdata.FSNRover25 = nansum(SNRwindow>25)/nansum(SNRwindow>-Inf);
                FLdata.FSNRover30 = nansum(SNRwindow>30)/nansum(SNRwindow>-Inf);
            end
            
            
            
            %% save Tables
            FLdataT = struct2table(FLdata);
            Varnames1=FLdataT.Properties.VariableNames;
            if ~exist('FLFinalT')
                FLFinalT = nan(max(Mrange),size(FLdataT,2)); %
                FLFinalT(i,:)=table2array(FLdataT);
            else
                FLFinalT(i,:)=table2array(FLdataT);
            end
            
            
            
            Varnames2=HistData.Properties.VariableNames;
            
            if ~exist('HistDataT')
                HistDataT= nan(max(Mrange),135);
                HistDataT(i,:)=table2array(HistData);
            else
                HistDataT(i,:)=table2array(HistData);
            end
           
            disp(['Completed getFL for study ' settings.savename '_' num2str(i) '.mat']);
            
            success(i)=1;
            
        else
            disp(['No Analyzed file ' settings.savename '_' num2str(i) '.mat']);
        end
    catch
        disp(['Failed FL Calculation for study ' settings.savename '_' num2str(i) '.mat']);
        
    end
end
FLFinalT=array2table(FLFinalT);
FLFinalT.Properties.VariableNames=Varnames1;
HistDataT=array2table(HistDataT);
HistDataT.Properties.VariableNames=Varnames2;

savefilename = [settings.workdir '\Summary\getFL' '_' datestr(now,'yyyymmdd HHMMSS')];
savefilename(end-6)='T';
save(savefilename,'FLFinalT','HistDataT','-v7.3')

savefilename = [settings.workdir '\Summary\getFL'];
save(savefilename,'FLFinalT','HistDataT','-v7.3')


delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['FL Calculation Complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext);
