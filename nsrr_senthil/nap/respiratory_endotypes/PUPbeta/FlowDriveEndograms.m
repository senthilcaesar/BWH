% function [BreathDataTable3 ,Data,DataCI,DataN,varlist,AHItotal,Fstates,...
%     Veupnea,medianV,VEdeciles,Vdrivedeciles,AHIdata,Fsupine,x_,VEdecilesUpper,VEdecilesLower,GGdata,GGpdeciles,GGtdeciles,EdecilesMean,TdecilesMean] = SummaryAnalysisOnePes(subj,settings)
global SAfontsize settings

%% data for endogram selection
   settings.selectstate=1
    WSUthreshold=0.0123  %
    WSLthreshold=0
    settings.constantEupnea=1
    settings.WakeSleepcriteria=1
state=settings.selectstate; %1=nrem1, 4 nrem, 5 rem, 8 ignore state


%% State, 1=nrem1, 2=nrem2, 3=nrem3, 4 = all nrem, 5 rem, 8 ignore state
%defaults
minFnrem1=-Inf;
minFnrem2=-Inf;
minFnrem3=-Inf;
switch settings.selectstate
    case 1
        maxFREM=0;
        minFrem=-Inf;
        minFnrem1=0.5;%settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7);
        PlotString = 'N1';
        hypok = [2];
    case 2
        maxFREM=0;
        minFrem=-Inf;
        minFnrem2=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7);
        PlotString = 'N2';
        hypok = [1];
    case 3
        maxFREM=0;
        minFrem=-Inf;
        minFnrem3=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7);
        PlotString = 'N3';
        hypok = [0];
    case 4 % at the moment this restricts to any window which does not have any REM.
        maxFREM=0; %1-settings.statetolerance;
        minFrem=-Inf;
        AHIstate=AHIdata(88);
        DurationInState=AHIdata(88-7);
        PlotString = 'NREM';
        hypok = [0 1 2];
    case 5 % at the moment this restricts to a window with >50% of REM, and the rest can be anything.
        maxFREM=Inf;
        minFrem=0.5; %settings.statetolerance;
        AHIstate=AHIdata(96);
        DurationInState=AHIdata(96-7);
        PlotString = 'REM';
        hypok = [3];
    case 8 % at the moment this does not restrict at all, and could include periods which are predominantly wake unless a wake threshold is applied elsewhere.
        maxFREM=Inf;
        minFrem=-Inf;
        AHIstate=AHIdata(80);
        DurationInState=AHIdata(80-7);
        PlotString = 'ALL';
        hypok = [0 1 2 3];
end

%% Find first and last windows containing sleep
TimeOfNightCriteria = ones(height(WinT),1);
if settings.selecttimeofnight
    XTiles=settings.selecttimeofnight_XTiles;
    NthXTile=settings.selecttimeofnight_NthXTile;
    
    FwakePerWindow=WinT.FWake;
    SleepWin1=find(FwakePerWindow<1,1,'first');
    SleepWinN=find(FwakePerWindow<1,1,'last');
    
    rangeWin=round((SleepWinN-SleepWin1+1)/XTiles);
    lowerWinNs=round(SleepWin1+((1:XTiles)-1)*rangeWin);
    upperWinNs=round(lowerWinNs+rangeWin-1);
    
    TimeOfNightCriteria=0*TimeOfNightCriteria;
    TimeOfNightCriteria(lowerWinNs(NthXTile):upperWinNs(NthXTile))=1;
end

%%
usemediannotmeanLG1=1;
% %

%% Position

% converting a single cell to table, but if cell array of tables,
% then handle differently, and don't {} convert  (DLM)
if isa(BreathDataTable,'cell') && all(size(BreathDataTable)==([1,1]))
    BreathDataTable=BreathDataTable{1};
end

% old position coding
PosWin = zeros(size(BreathDataTable,2),1);
Fsupinetemp = zeros(size(BreathDataTable,2),1);
% new position coding
try
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{subj});
catch % RMA added--for running single subject A to Z mode
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(subj,:));
end

for xx=1:size(BreathDataTable,2)
    if size(BreathDataTable{xx},1)<=1
        PosWin(xx)=NaN;
        Fsupinetemp(xx)=NaN;
        continue
    end
    temp = BreathDataTable{xx}.pos_B;
    temp2 = PositionTranslator(positioncodes,settings.positioncodesout,temp);
    PosWin(xx) = mode(temp2);
    Fsupinetemp(xx) = nanmean(PositionSelector(temp2,'Supine'));
end
Poscriteria = PositionSelector(PosWin,settings.selectposition);
str=[num2str(nnz(Poscriteria)), ' of ', num2str(length(Poscriteria)) ,' windows in selected position'];
disp(str);


%% fractions of each sleep state
%SleepData(winNum+1,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
if 1
    Funknownsleep=round(1-sum([WinT.FWake WinT.FNREM WinT.FREM],2),2);
else
    Funknownsleep=(1-sum([WinT.FWake WinT.FNREM WinT.FREM],2)); %old code, removed 5/7/2021. delete option in future
end
containsunknownsleep=1*(Funknownsleep>0);


%% This criteria defines a minimum proportion of particular non-REM sleep states.
nremXcriteria= WinT.FNREM1>=minFnrem1 &...
    WinT.FNREM2>=minFnrem2 & ...
    WinT.FNREM3>=minFnrem3;


%%

%% criteria for LG and UA phenotype
minNeventsUA = settings.minNeventsUA; %has been default 1, was default 0 for a long time previously
maxwakethresLG = settings.maxwakethresLG; % 30
maxwakethresUA = settings.maxwakethresUA; % 300; more tolerant to more wake for UA measures

try %added SS 11/2/2021
    CPAPoffCriteria = WinT.CPAPoff==1;
catch
    CPAPoffCriteria = zeros(1,height(WinT))+1;
    disp('warning: CPAPoff detection failed, assuming CPAP is off always');
end


criteriaAll = criteriacomb & ...
    Poscriteria==1 & ...
    CPAPoffCriteria==1 & ...
    nremXcriteria==1 & ...
    WinT.FREM<=maxFREM &...
    WinT.FREM>=minFrem & ...
    containsunknownsleep==0 & ...  %recently added contains unknownsleep
    TimeOfNightCriteria==1;
if settings.verbose
    disp(['criteriaAll = ', num2str(nnz(criteriaAll)), ' / ', num2str(length(criteriaAll))]);
end

criteriaUA = criteriaAll & ...
    WinT.Nevents>=minNeventsUA & ... % nnz(isfinite(N_events))
    WinT.LongestWake<=maxwakethresUA;
if settings.verbose
    disp(['criteriaUA = ', num2str(nnz(criteriaUA)), ' / ', num2str(length(criteriaUA))]);
end



%% Bring in LG/VRA/ArThres window data

WinT.VRA(WinT.Narousals<1)=NaN;
WinT.VRA2(WinT.Narousals<1)=NaN;
WinT.Poscriteria = Poscriteria;
WinT.nremXcriteria = nremXcriteria;
WinT.TimeOfNightCriteria = TimeOfNightCriteria;
WinT.criteriacomb = criteriacomb;
WinT.Fsupine = Fsupinetemp;
WinT.EndogramCriteria = criteriaUA;



%%

%%
% Fstates3=[mean(WinT.FNREM1(criteriaUA)) mean(WinT.FNREM2(criteriaUA)) mean(WinT.FNREM3(criteriaUA)) mean(WinT.FREM(criteriaUA))]';
% Fsupine3 = nanmean(Fsupinetemp(criteriaUA));

NwindowsForUA=sum(criteriaUA);
str=['Using ', num2str(NwindowsForUA), ' of ', num2str(length(criteriaUA)) ,' windows for Endograms'];
disp(str);

for ii=(Nperwindowvars+1):length(varlist)
    eval([varlist{ii} 'N=NwindowsForUA;']);
end
fontsize_= SAfontsize;

%% UA phenotype using model drive
if 0
    subplot(1,3,1);
    %make subset table
    criteriaUAI = find(criteriaUA==1);
    criteriaRow = sum((BreathDataTable2.Win == criteriaUAI'),2)>0;
    BreathDataTable3 = BreathDataTable2(criteriaRow==1,:);
    
    [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,Vcomp,VcompCI]=SummaryAnalysisOne_UAmodel(length(criteriaUAI),BreathDataTable3,ArThres,varlist,settings.plotfigs,settings);
end

%% UA Phenotype, gold standard drive

try
    
    ArThresActive=NaN;
 
    
    if NwindowsForUA>0
        
      
        
        
        
        
        
        if settings.plotfigs
            plot(1,1,'marker','none'); hold('on');
        end
        hold('on');
        %%[16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak 24:FlowEdi, 23:FlowPes]
        dividedrivebyttot=0;
             
        %% Wake data, used to normalize drive to eupneic flow (later)
        %calculate "g" = nanmedian(y(a1)./x(a1)); y = flow, x=drive, a1=wake breaths
        BreathDataTable2;
        
        x=BreathDataTable2.(settings.DriveSignal);
        y = BreathDataTable2.VE;
        a = BreathDataTable2.AR3;
        hyp = BreathDataTable2.hypnog_B;
        
        win = floor(BreathDataTable2.UniqueID/1000);
        
        %need y, x, a1, win
        if settings.scoredarousalsinwake==0
            a=(a==1)|(hyp==4);
        end
        [ad]=howfarawayfromsleep(BreathDataTable2.AR3,win);
        %find breaths during wakefulness and arousals
        minNwakebreaths = 50;
        hh=hist(ad,[1:11]); hh(end)=NaN; th=find(hh>minNwakebreaths,1,'last');
        threshold = min([4 th]);
        a1 = ad>threshold;
        g = nanmedian(y(a1)./x(a1));
        x_wake=nanmedian(x(a1));
        y_wake=nanmedian(y(a1));
        G_N = length(y(a1)./x(a1));
        
        %% GET BREATH-BY-BREATH SLEEP DATA
        
        
        win1 = floor(BreathDataTable2.UniqueID/1000);
        criteriaRow = sum((win1 == (find(criteriaUA==1))'),2)>0;
        BreathDataTable2.Include=criteriaRow;
     %   BreathDataTable3=BreathDataTable2(BreathDataTable2.Include==1,:);
        
        t=BreathDataTable2.Time_start(criteriaRow);
        x=BreathDataTable2.(settings.DriveSignal)(criteriaRow);        
        y = BreathDataTable2.VE(criteriaRow);
        a = BreathDataTable2.AR3(criteriaRow);
        veup = BreathDataTable2.Veup(criteriaRow);
        win = floor(BreathDataTable2.UniqueID(criteriaRow)/1000);
        t2=BreathDataTable2.Time_end(criteriaRow) - BreathDataTable2.Time_start(criteriaRow);
        e=BreathDataTable2.E1(criteriaRow);
        nota = BreathDataTable2.notAR3(criteriaRow);
        GGp = BreathDataTable2.GGpeak(criteriaRow);
        GGt = BreathDataTable2.GGtonic(criteriaRow);
                
       
        if settings.constantEupnea %the single-value for Veupnea will stay constant regardless of selected states
            criteriaEupnea = WinT.LongestWake<=settings.maxwakeUA & WinT.FREM<=0; %use all data except wake and REM; should adapt with window length
            win1 = floor(BreathDataTable2.UniqueID/1000);
            criteriaRowEupnea = sum((win1 == (find(criteriaEupnea==1))'),2)>0;
            veup1 = BreathDataTable2{criteriaRowEupnea==1,'Veup'};
        end
        
        if dividedrivebyttot
            x=x./t2;
        end
        
        %normalize data
        if ~settings.normalizeusingconstantEupnea %eupnea stays constant (does not adapt to signal size), window by window
            y = y./veup; %%window-based normalization. note Y is unnormalized VE from LGfromflow
            y_ = nanmedian(veup);
        else
            if settings.constantEupnea
                y = y./nanmedian(veup1); %study-based normalization, but removes windows predominately wake/REM. note Y is unnormalized VE from LGfromflow
                y_ = nanmedian(veup1);
            else
                y = y./nanmedian(veup); %never actually used this, study-based normalization of all windows Veupnea
                y_ = nanmedian(veup);
            end
        end
        Veupnea = y_*60; %
          y_n = y_;
            x_ = y_/g; %xeup: drive eupneic value, veup divided by wake scale factor
       
        
        if ~exist('EupneicDriveIsUnity')
            EupneicDriveIsUnity=0;
        end
        if settings.DriveSignal=="VdriveModel"
            EupneicDriveIsUnity=1;
        end
        if EupneicDriveIsUnity==1
            x_=1;
        end
        
        x_n=x_
        
        %Arousal threshold
        
        %% Remove wakefulness breaths for further analysis
        if scoredarousalsinwake
            criteriabreath = BreathDataTable2.Include & BreathDataTable2.notAR3==1;
        else
            criteriabreath = BreathDataTable2.Include & BreathDataTable2.notAR3==1 & BreathDataTable2.hypnog_B~=4;
        end
        %               %
    %    criteriabreath = BreathDataTable2.Include2 %ones(length(x),1)==1; %default
        
        if settings.breathlevelinclusion %excludes individual breaths if not in appropriate state as selected; new default
            %hypok = [0 1 2 3];
            criteriabreath = criteriabreath & sum(BreathDataTable2.hypnog_B==hypok,2)>0; % criteria for sleep stages in breath analysis.
        end
        
        %settings.TonicREMonly=1;
        %settings.PhasicREMonly=1;
        disp(['N breaths: ' num2str(sum(criteriabreath))]);
        if state==5
            if isfield(settings,'PhasicREMonly') && settings.PhasicREMonly==1 %Phasic only
                NbreathsPhasic = sum(BreathDataTable2.REMphasic==1&criteriabreath)
                NbreathsTotal = sum(criteriabreath)
                criteriabreath = criteriabreath & (BreathDataTable2.REMphasic==1);
            end
            
            if isfield(settings,'TonicREMonly') && settings.TonicREMonly==1 %Tonic REM only
                NbreathsTonic = sum(REMtonic==1&criteriabreath)
                NbreathsTotal = sum(criteriabreath==1)
                criteriabreath = criteriabreath & (BreathDataTable2.REMtonic==1);
            end
        end
        
        % fun new options :)
        if isfield(settings,'EndogramDriveslopecriteria')
            if settings.EndogramDriveslopecriteria>0
                criteriabreath = criteriabreath & (BreathDataTable2.DriveSlope>0);
            elseif settings.EndogramDriveslopecriteria<0
                criteriabreath = criteriabreath & (BreathDataTable2.DriveSlope<0);
            end
        end
        
        if isfield(settings,'EndogramEventcriteria')
            if settings.EndogramEventcriteria==1 %event only
                criteriabreath = criteriabreath & (BreathDataTable2.Etype>0);
            elseif settings.EndogramEventcriteria==-1 %non event only
                criteriabreath = criteriabreath & (BreathDataTable2.Etype==0);
            elseif settings.EndogramEventcriteria==2 %apnea
                criteriabreath = criteriabreath & (BreathDataTable2.Etype==2);
            elseif settings.EndogramEventcriteria==3 %hypopnea
                criteriabreath = criteriabreath & (BreathDataTable2.Etype==4);
            end
        end
        
        if isfield(settings,'StableBreathingcriteria')
            if settings.StableBreathingcriteria==0
                criteriabreath = criteriabreath & (BreathDataTable2.StableBreathing==settings.StableBreathingcriteria);
            elseif settings.StableBreathingcriteria>0
                criteriabreath = criteriabreath & (BreathDataTable2.StableBreathing>=settings.StableBreathingcriteria);
            end
        end
        
          if isfield(settings,'WakeSleepcriteria')
         
         criteriabreath = criteriabreath & (BreathDataTable2.WakeSleep<=WSUthreshold)& BreathDataTable2.WakeSleep>=WSLthreshold;
       
          end
        
        
        BreathDataTable2.Include2 = criteriabreath; 
        idx=BreathDataTable2.Include2( BreathDataTable2.Include==1);
        disp(['N breaths: ' num2str(sum(criteriabreath))]);
        
        x = x(criteriabreath(BreathDataTable2.Include==1));
        y = y(criteriabreath(BreathDataTable2.Include==1));
        GGp = GGp(criteriabreath(BreathDataTable2.Include==1));
        GGt = GGt(criteriabreath(BreathDataTable2.Include==1));
        veup = veup(criteriabreath(BreathDataTable2.Include==1));
        
        
        medianV = 100*median(y);
        
        x=x/x_;
        
        
        if length(x)<20
            length(x)
            error('not enough breaths left to analyze')
        end
        %NcileLimits=[5 95];
        %settings.Nciles=131;
        
        %NcileLimits=[2 80];
        %settings.Nciles=119;
        %
        %%   Endogram Analysis and Plot
        
        if ~isfield(settings,'SummaryEndogramFilt121')
            filt121 = 1;  %default on, since 20201027
        else
            filt121 = settings.SummaryEndogramFilt121;
        end
        
        ploton=settings.plotfigs;
        if ploton==1
            figure(1);
           subplot(1,3,1);
            hold on
        end
        
        if ~isfield(settings,'NcileLimits')
            NcileLimits=[];
        else
            NcileLimits=settings.NcileLimits;
        end
        
        if ~DriveincmH2OnotEupnea
            if settings.plotfigs
                plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
                hold('on');
                plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
            end
            [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,VEdecilesUpper,VEdecilesLower]=VEVdriveAnalysis(100*x,100*y,100*ArThresActive,ploton,100,settings.Nciles,Nbootstrap,NcileLimits,filt121);
        else
            if settings.plotfigs
                plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
                hold('on');
                plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
            end
            [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,VEdecilesUpper,VEdecilesLower]=VEVdriveAnalysis(x*x_,100*y,ArThresActive*x_,ploton,x_,settings.Nciles,Nbootstrap,NcileLimits,filt121);
        end
        VpassiveCI=VpassiveCI(:);
        VactiveCI=VactiveCI(:);
        
        
        if settings.plotfigs
            set(gca,'fontsize',fontsize_)
            ylim([0 max([prctile(100*y,95) 100])]);
            ylabel('Ventilation, %eupnea');
            if ~DriveincmH2OnotEupnea
                xlim([0 max([prctile(100*x,95) 300])]);
                xlabel('Drive, %eupnea');
            else
                xlim([0 max([prctile(x*x_,95) 8*x_])]);
                xlabel('Drive, absolute units');
            end
            %hold('off');
            %try title(['AHI=',num2str(round(AHIstate),'%u')]); end
            
            subplot(2,6,9);
            xedges3=0.7:0.1:3;
            
            dx=0.1;
            h=histogram(100*ArThresvalues/x_,100*xedges3,'normalization','probability','EdgeAlpha',0);
            
            box('off');
            xlim(100*[min(xedges3) max(xedges3)]);
            %ylim([0 max([0.2;h.Values(:)])]);
            xlabel('ArThres');
            hold('on');
            
            
            subplot(2,6,10);
            xedges=0:0.1:1.4;
            xedges2=0:0.1:4.0;
            dx=0.1;
            h=histogram(100*y(x>(1-dx)&x<(1+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
            box('off');
            xlim(100*[min(xedges) max(xedges)]);
            ylim([0 max([0.2;h.Values(:)])]);
            xlabel('Vpassive');
            hold('on');
            
            
            subplot(2,6,11);
            h=histogram(100*y(x>(ArThresActive-dx)&x<(ArThresActive+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
            box('off');
            xlim(100*[min(xedges) max(xedges)]);
            ylim([0 max([0.2;h.Values(:)])]);
            xlabel('Vactive');
            hold('on');
            
            
        end
        Vcomp = Vactive-Vpassive;
        VcompCI = [VactiveCI(1)-VpassiveCI(2);VactiveCI(2)-VpassiveCI(1)];
        VcompCI = sort(VcompCI);
        
        figure(2); clf(2);
        varlist2 = {'veonvdrive','veonvwake','driveondrivewake'};
        xlimupper = [140 140 300];
        
        for i=1:length(varlist2)
            if 0
                h=histogram(eval(varlist2{i}),0:10:xlimupper(i),'normalization','probability','EdgeAlpha',0);
            else
                Data = eval(varlist2{i});
                dStep=10;
                Centers=5:dStep:(xlimupper(i)-dStep/2);
                Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2); Edges(end)=Inf;
                [h1,edges] = histcounts(Data,Edges);
                bar(Centers,h1/sum(h1),'EdgeAlpha',0,'FaceAlpha',0.6,'BarWidth',1);
                set(gca,'xtick',[0:50:300],'tickdir','out');
                
            end
            hold('on');
            eval([varlist2{i} '(isnan(' varlist2{i} '))=[];']);
            if Nbootstrap>0
                eval([varlist2{i} 'CI=bootci(Nbootstrap,@median,' varlist2{i} ');']);
            end
        end
        box('off');
        xlim([0 max(xlimupper)]);
        
        xlabel('Value, %');
        ylabel('Proportion of Sleep');
        
        ylims=get(gca,'YLim');
        i=1;
        plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0 0.45 0.74]);
        i=2;
        plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0.85 0.33 0.1]);
        i=3;
        plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0.93 0.69 0.13]);
        
        h=legend('Actual Airflow / Intended Airflow','Actual Airflow / Wake Baseline','Intended Airflow / Wake Baseline');
        
        %% GG
        try
            
            ploton=settings.plotfigs;
            figure(4);
            
            
            if ~DriveincmH2OnotEupnea
                if settings.plotfigs
                    plot([100 100],[0 prctile(GGp,95)],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    %plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGppassive,GGpactive,GGpdeciles,GGpdeciles_drive,GGppassiveCI,GGpactiveCI,GGpdecilesUpper,GGpdecilesLower,GGpdecilesMean]=VEVdriveAnalysis(100*x,GGp,100*ArThresActive,ploton,100,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            else
                if settings.plotfigs
                    plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGppassive,GGpactive,GGpdeciles,GGpdeciles_drive,GGppassiveCI,GGpactiveCI,GGpdecilesUpper,GGpdecilesLower,GGpdecilesMean]=VEVdriveAnalysis(x*x_,GGp,ArThresActive*x_,ploton,x_,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            end
            ylim([0 max([GGpdecilesUpper;5])]);
            set(gca,'fontsize',fontsize_);
            
            %ylim([0 max([prctile(100*y,95) 100])]);
            ylabel('GG peak, %max');
            if ~DriveincmH2OnotEupnea
                xlim([0 max([prctile(100*x,95) 300])]);
                xlabel('Drive, %eupnea');
            else
                xlim([0 max([prctile(x*x_,95) 8*x_])]);
                xlabel('Drive, absolute units');
            end
            
            hold('on');
            if ~DriveincmH2OnotEupnea
                if settings.plotfigs
                    plot([100 100],[0 prctile(GGp,95)],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    %plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGtpassive,GGtactive,GGtdeciles,GGtdeciles_drive,GGtpassiveCI,GGtactiveCI,GGtdecilesUpper,GGtdecilesLower]=VEVdriveAnalysis(100*x,GGt,100*ArThresActive,ploton,100,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            else
                if settings.plotfigs
                    plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGtpassive,GGtactive,GGtdeciles,GGtdeciles_drive,GGtpassiveCI,GGtactiveCI,GGtdecilesUpper,GGtdecilesLower]=VEVdriveAnalysis(x*x_,GGt,ArThresActive*x_,ploton,x_,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            end
            
            ylim([0 max([GGpdecilesUpper;5])]);
            
            GGppassiveCI=GGppassiveCI(:);
            GGtactiveCI=GGtactiveCI(:);
            
            set(gca,'fontsize',fontsize_);
            
            %ylim([0 max([prctile(100*y,95) 100])]);
            ylabel('GG, %max');
            if ~DriveincmH2OnotEupnea
                xlim([0 max([prctile(100*x,95) 300])]);
                xlabel('Drive/Effort, %eupnea');
            else
                xlim([0 max([prctile(x*x_,95) 8*x_])]);
                xlabel('Drive/Effort, %absolute');
            end
            GGdata = [GGppassive;GGpactive;GGtpassive;GGtactive];
            
        catch me
            %disp(me.message);
            disp('GG analysis failed. Possibly no GG data available')
            
            GGdata = [NaN;NaN;NaN;NaN];
            GGpdeciles = NaN*ones(settings.Nciles,1);
            GGtdeciles = NaN*ones(settings.Nciles,1);
            close(4);
        end
        
        
        
        %% Event probability
                
        e=e(Ikeep);
        e=e(criteriabreath);
        [Epassive,Eactive,Edeciles,Edeciles_drive,EpassiveCI,EactiveCI,EdecilesUpper,EdecilesLower,EdecilesMean]=VEVdriveAnalysis(x*x_,e,ArThresActive*x_,0,x_,settings.Nciles,Nbootstrap,NcileLimits);
        
        if 1 %add colorbar to individual plot
            figure(1)
            subplot(1,3,1);
            DispEventLikelihoodOnEndogram(Vdrivedeciles,VEdeciles,1-EdecilesMean,0)
        end
        
        decileTime=t(Ikeep);
        decileTime=decileTime(criteriabreath);
        [Tpassive,Tactive,Tdeciles,Tdeciles_drive,TpassiveCI,TactiveCI,TdecilesUpper,TdecilesLower,TdecilesMean]=VEVdriveAnalysis(x*x_,decileTime,ArThresActive*x_,0,x_,settings.Nciles,Nbootstrap,NcileLimits);
        
        
    else
        Vpassive=NaN;
        Vactive=NaN;
        Vcomp=NaN;
        VpassiveCI=[-Inf;Inf];
        VactiveCI=[-Inf;Inf];
        VcompCI=[-Inf;Inf];
        GGdata = [NaN;NaN;NaN;NaN];
        VEdeciles = NaN*ones(settings.Nciles,1);
        Vdrivedeciles = NaN*ones(settings.Nciles,1);
        GGpdeciles = NaN*ones(settings.Nciles,1);
        GGtdeciles = NaN*ones(settings.Nciles,1);
        VEdecilesUpper= NaN*ones(settings.Nciles,1);
        VEdecilesLower= NaN*ones(settings.Nciles,1);
        x_ = NaN;
        y_ = NaN;
        x_wake = NaN;
        y_wake = NaN;
        Veupnea = NaN;
        medianV = NaN;
        EdecilesMean = NaN*ones(settings.Nciles,1);
        TdecilesMean = NaN*ones(settings.Nciles,1);
    end
    
catch me
    disp('failed UA analysis');
    disp(me.message);
    Vpassive=NaN;
    Vactive=NaN;
    Vcomp=NaN;
    VpassiveCI=[-Inf;Inf];
    VactiveCI=[-Inf;Inf];
    VcompCI=[-Inf;Inf];
    GGdata = [NaN;NaN;NaN;NaN];
    VEdeciles = NaN*ones(settings.Nciles,1);
    Vdrivedeciles = NaN*ones(settings.Nciles,1);
    GGpdeciles = NaN*ones(settings.Nciles,1);
    GGtdeciles = NaN*ones(settings.Nciles,1);
    VEdecilesUpper= NaN*ones(settings.Nciles,1);
    VEdecilesLower= NaN*ones(settings.Nciles,1);
    x_ = NaN;
    y_ = NaN;
    x_wake = NaN;
    y_wake = NaN;
    Veupnea = NaN;
    medianV = NaN;
    EdecilesMean = NaN*ones(settings.Nciles,1);
    TdecilesMean = NaN*ones(settings.Nciles,1);
    
end
%%
disp(['Eupneic VE, Vdrive sleep: ' num2str([y_ x_ ]) '; wake: ' num2str([y_wake x_wake]) ' (raw)']);

%% Save
criteriaset=[1 1 1 2 3 3 3 3 2];
clear Data DataCI Fstates Fsupine
for i=1:length(varlist)
    Data(i) = eval(varlist{i});
    if settings.getCIs
        DataCI(:,i) = eval([varlist{i} 'CI']);
    else
        DataCI(:,i) = NaN;
    end
    DataN(i) = eval([varlist{i} 'N']);
    Fstates(:,i) = eval(['Fstates' num2str(criteriaset(i))]);
    Fsupine(:,i) = eval(['Fsupine' num2str(criteriaset(i))]);
end
%pause