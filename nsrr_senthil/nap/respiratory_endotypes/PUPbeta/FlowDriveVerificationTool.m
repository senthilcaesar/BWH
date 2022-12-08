
disp('Requires Completed Event Analysis. Runs Studies with Analyze=1');

%% Which Data
dirc=[settings.workdir 'Analyzed\EventAnalyzed\'];
settings = ImportSettings(settings,AMasterSpreadsheet);
%[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C58');
%settings.savename = char(raw{1});
%[num,~,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
%analyzelist = logical(num(:,2));

%% Loop

M=1:length(settings.patients);

if exist('MrangeOverride')==1
    M=MrangeOverride;
end

clear Tout;
minflowdrive = 0.1;
%M=15;
for i=M %1:11 13:19 fails: 20 23 24 33
    i
    %if analyzelist(i)
    try
    success(i)=1;
   filen = [dirc settings.savename '_' num2str(i)];
   load(filen);
   
   T = BreathDataTableFulls;
   
   %QC on FlowDrive
   T.FlowDrive(T.FlowDrive<minflowdrive)=minflowdrive;
   T.FlowDrive(T.FlowDrive>1.2)=1.2;
   T.FlowDrive(T.ApneaB==1)=minflowdrive;
   T.FlowDrive(T.VI<minflowdrive & T.Ecentralapnea==0)=minflowdrive;
   
   Nrecoverybreaths = 6;
   Events = 1*(T.Etype>0);
   EventLastBreath = find([(diff(Events)==-1)]);
   EventRecovery = 0*Events;
   for j=1:length(EventLastBreath)
       temp = EventLastBreath(j);
       temp2 = temp + [1:Nrecoverybreaths];
       temp2(temp2>height(T))=[];
       wintemp = T.Time0(temp);
       wintemp2 = T.Time0(temp2);
       samewin = wintemp2==wintemp;
       temp2(samewin==0)=[];
       EventRecovery(temp2)=1;
   end
   EventRecoveryHighV = EventRecovery;
        EventRecoveryHighV(T.VI<1.5)=0;
        %EventRecoveryHighV(T.VI>3)=0;
   EventsHypopneas = 1*(T.Etype==4);
   
   if isfield (settings,'AutoScoredEventsAllCentral')& settings.AutoScoredEventsAllCentral==1
       % for use in central hypopneas- added 3/14/2022- RMA
        EventsHypopneas = 1*(T.Etype==6);
   end
        EventsHypopneas(T.VI>0.5)=0;
        EventsHypopneas(T.VI<=0.1)=0;
   dataN = [sum(EventsHypopneas) sum(EventRecoveryHighV)]
   dataFD = [nanmean(T.FlowDrive(EventsHypopneas==1)) nanmean(T.FlowDrive(EventRecoveryHighV==1))]
   
   T1 = array2table(T.FlowDrive(EventsHypopneas==1));
   T1.Properties.VariableNames = {'FlowDrive'};
   T1.Hypopnea = ones(height(T1),1);
   
   T2 = array2table(T.FlowDrive(EventRecoveryHighV==1));
   T2.Properties.VariableNames = {'FlowDrive'};
   T2.Hypopnea = zeros(height(T2),1);
   
   T3 = [T1;T2];
   T3.Subj = ones(height(T3),1) + i;
   
   if exist('Tout')   %I = GiantTable;
       Tout = [Tout;T3];
   else
       Tout = T3;
   end
    catch me
       success(i)=0;
    end
 %   end
end

%% Plot

N=height(Tout);

figure(1)
bw=0.05/(2^log10(N/100));
trans = 0.7 - 0.2* log10(N/100);
    trans(trans>1)=1;
    trans(trans<0.1)=0.1;
mn = 1-(1./(N.^0.5));
data = Tout.FlowDrive(Tout.Hypopnea==1);
catIdx = ones(length(data),1);
[~,dataOut,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
data = Tout.FlowDrive(Tout.Hypopnea==0);
catIdx = ones(length(data),1);
[~,dataOut2,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
close(1)

figure(1); clf(1); set(gcf,'color',[1 1 1]);
scatter(dataOut(:,1),dataOut(:,2),10,[0.9 0.2 0.05],'filled','markerfacealpha',trans);
hold on
scatter(dataOut2(:,1)+1,dataOut2(:,2),10,fliplr([0.9 0.2 0.05]),'filled','markerfacealpha',trans);
set(gcf,'position',[ 991   545   249   433])
[splithres,AUC,~,~,~,~]=ROCAUCSEM(Tout.Hypopnea,Tout.FlowDrive)

ylabel('Flow:Drive')
set(gca,'xtick',[1 2],'xticklabels',{'Hypopnea','Recovery'},'fontsize',12);
xtickangle(45);

title(['AUC: ' num2str(AUC,2)])

%% Save Figure
if ~(exist([settings.workdir,'QC'], 'dir') == 7)
        mkdir([settings.workdir,'QC']);
end
saveas(1,[settings.workdir 'QC\FlowDriveVerification'],'tiff')
save([settings.workdir 'QC\FlowDriveVerification'],'Tout');


