%% startup
clear
close all

%% specify a user and individual settings
user = 'Scotty'; % 'Scotty'; % 
switch user
    case 'Scotty'
        % path to code
        CodeDirectory = cd;
        addpath(CodeDirectory);
        cd(CodeDirectory);

        % path to additional code
        addpath('J:\PEOPLE\FACULTY\SANDS\MatlabFunctionDatabase');
        addpath('J:\PEOPLE\FACULTY\SANDS\PUPFlowandPes\PUPbeta 20161221 FlowAndPes\old matlab functions');

        % select data set.
        DataSet = 'FlowDrive'; % 'OralAppliance';  % 
        
        LocalDataDir = 'E:\Work\Projects SS\QAO';
        LocalData = [LocalDataDir,'\FS_FD_1_noSignals.mat'];
        savestr_SVMPrep = [LocalDataDir, '\SA_FD_Testing_Regression_Prep.mat']; 
        savestr_SVMRun = [LocalDataDir, '\SA_FD_Testing_Regression_Run.mat'];

    case 'dwayne'
        % path to code
        CodeDirectory = 'C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO';
        addpath(CodeDirectory);
        cd(CodeDirectory);

        % select data set.
        DataSet = 'FlowDrive';%'OralAppliance';%
        
        LocalDataDir = 'C:\PSG_Data\QAO';
        LocalData = [LocalDataDir,'\FS_FD_1_noSignals.mat'];
        savestr_SVMPrep = [LocalDataDir, '\SA_FD_Testing_Regression_Prep.mat']; 
        savestr_SVMRun = [LocalDataDir, '\SA_FD_Testing_Regression_Run.mat'];
end

switch DataSet
    case 'FlowDrive'
        DriveAnalysis = 1; % set to one for DriveAnalysis
    case 'OralAppliance'
        DriveAnalysis = 0; % set to 0 for non-drive data
    otherwise
        disp('Please select a data set to work with');
        return;
end

%options
ShowFigures = 0;
ShowExtraFigures = 0;

% selectively load variables of interest
load(LocalData,'SleepData','LG_QualityInfo','DataOut',...
    'BreathDataTable','BreathFLDataTable','LocalSignals');

n_pts=length(SleepData); % how many pts are we processing

% % set defaults for figures
% set(groot,'defaultAxesfontname','arial narrow');
% set(groot,'defaultFigureColor',[1 1 1]); % set background colour to white
% set(groot,'defaultAxesBox','off');
% set(groot,'defaultAxesTickDir','out');
% set(groot,'Units','inches');
% r=groot;
% % then later, use
% fig = r.Chilren;

%% Description of loaded variables
% each variable is a 1 x p cell for each of the p patients
% in any given patient, there may be w windows analyzed
% each window has b breaths, and e events, and t time
% 'SleepData' {1,p}(w,7)
% 'LG_QualityInfo' {1,p}(w,13)
% 'DataOut' {1,p}{1,w}(b,24)
% 'BreathDataTable' {1,p}{1,w}(b,22)
% 'BreathFLDataTable' {1,p}{1,w}(b,100)
% 'LocalSignals' {1,p}{1,w}(t,2) col 1 is Flow, col 2 is Edi
% Note that the BreathDataTable and BreathFLDataTable may be
% shorter than the others if no data in the end windows

% Variables available, but currently not loaded
% 'MiscEpochData' {1,p}{1,w}, each being a single value
% 'AnalysisIndex' {1,p} each being (w,2) pair of start stop values
% 'LGplusinfo' {1,p}(w,17)
% 'EventsInfo' {1,p}(w,10)
% 'ArousalDat' {1,p}{1,w}(e,2) where e is the number of arousals?
% 'AHIdata' {1,p}(1,128)
% 'CPAPData' {1,p}(w,4)

%% Processing steps
% for each pt
% get the breath data, VE, VDrive_Edi, VDrive_Pes,
% and all the other breath data that is later matched up
%  (BB_time,BB_Ttot,Ar,win,notAr,Veup,hypnog,pcw,Evt,FL)
% then, determine the best VE/Vdrive ratio for Edi and Pes
%  this is done with moving median value over wake breaths
% all breaths are then divided by these ratios,
% values >1 should occur infrequently, and could be excluded
% 1 would indicate good breathing, and <1 indicates FL.
% remove wakeful breaths with FL (<0.75) (only those >4 breaths from sleep)
%
%
% things to check:
%  breaths are not normalised to any eupnea
%  drive values (Edi nor Pes) are not divided by Ttot
%  what is the threshold for upper cut-off of normalised VE_Vdrive
%    given that >1 means VE was greater than Vdrive
%  what is the lower cut-off for normalised VE_Vdrive, in wake breaths
%    given we don't want wake FL breaths to skew results
%   and.. should these be excluded prior to or after determining the ratio
%

%% declare variables
CombinedPtData = [];
Nwindows=NaN(n_pts,1);
threshold=NaN(n_pts,1);

%% patient processing
for pt=1:n_pts
    str = ['Processing patient ', num2str(pt)];
    disp(str);
    %     if Exclude(n)==1||LG1_N(n)<=minNwindows
    %         continue
    %     end
    try
        supinecodes = [0 2 -5]; %Profusion:[1],Spike:[0 2]
        maxwakethres = 360;  %minwakearthres = 60;
        minNevents = 0;
        maxFREM=0;
        N_events = LG_QualityInfo{pt}(:,2);
        Pos = LG_QualityInfo{pt}(:,5);
        Fwake = SleepData{pt}(:,1);
        FREM = SleepData{pt}(:,6);
        minwake = SleepData{pt}(:,7);
        criteria_type = 3;
        switch criteria_type
            case 1
                % apply criteria to select windows for analysis
                criteria = N_events>=minNevents & minwake<maxwakethres & FREM<=maxFREM ; %& (~isnan(Pos)&supine==1)
            case 2
                criteria = N_events>=minNevents & (Pos==0|Pos==2) & minwake<=maxwakethres;
                supine = 0*Pos;
                for j=1:length(supinecodes)
                    supine = supine + 1*(Pos==supinecodes(j));
                end
            case 3
                % get all data, (all wake, part wake part sleep, and all sleep)
                criteria = minwake>=0; %(~isnan(Pos)&supine==1)
            otherwise
                criteria = [];
        end
        Nwindows(pt)=sum(criteria);
        if sum(criteria)==0
            continue
        end
        % get the breath data, VE, VDrive_Edi, VDrive_Pes,
        [BB_time,VdriveEdi,VdrivePes,VE,Ar,win,notAr,Veup,hypnog,pcw,BB_Ttot,FL,EType] = ...
            VE_VdriveFLArray(pt,criteria,BreathFLDataTable, BreathDataTable);
        % and all the other breath data that is later matched up
        %  (BB_time,BB_Ttot,Ar,win,notAr,Veup,hypnog,pcw,Evt,FL)

        [ad]=howfarawayfromsleep(Ar,win); %replace a with (a==1)&(hyp==4) if arousals not scored in wake...
        if 1
            minNwakebreaths = 50;
            hh=hist(ad,[1:11]); hh(end)=NaN; th=find(hh>minNwakebreaths,1,'last');
            threshold(pt) = min([4 th]);
        else
          %find breaths during wakefulness and arousals
            threshold(n) = 2;
        end
        a1 = ad>threshold(pt); % threshold # breaths away from sleep
        a2 = ad>=2; % at least two breaths away from sleep

        if DriveAnalysis

        % Figure, VE/Vdrive ratio with Edi and Pes
        % and an approximation of drive ratio
        if ShowFigures
            figure(100); clf(figure(100));
            ax100(1) = subplot(2,1,1);
            plot(BB_time, VE./VdriveEdi,'kd');hold('on'); % plot everything
            plot(BB_time(a2),VE(a2)./VdriveEdi(a2),'rd'); % at least two breaths away from sleep
            plot(BB_time(a1),VE(a1)./VdriveEdi(a1),'gd');  % more than ~4 breaths away from sleep
            refline(0,nanmedian(VE(a2)./VdriveEdi(a2)));
            ylabel('VE / VdriveEdi');

            ax100(2) = subplot(2,1,2); linkaxes(ax100, 'x')
            plot(BB_time, VE./VdrivePes,'kd');hold('on'); % plot everything
            plot(BB_time(a2),VE(a2)./VdrivePes(a2),'rd');
            plot(BB_time(a1),VE(a1)./VdrivePes(a1),'gd');
            refline(0,nanmedian(VE(a2)./VdrivePes(a2)));
            ylabel('VE / VdrivePes');
        end

        % then, determine the best VE/Vdrive ratio for Edi and Pes
        %  this is done with moving median value over the wake breaths
        %  that are more than 2 breaths away from sleep
        % edi first
        temp = BB_time(a2);[temp,ia,ic] = unique(temp);
        data_edi = VE(a2)./VdriveEdi(a2);
        data_edi = data_edi(ia);
        data_edi(data_edi>prctile(data_edi,99))=NaN;
        data_edi(data_edi<prctile(data_edi,1))=NaN;
        temp(isnan(data_edi))=[]; data_edi(isnan(data_edi))=[];
        maxt = 1800;
        Gmedian_edi=NaN*temp;
        for j=1:length(temp)
            temp2 = abs(temp(j)-temp);
            weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
            Gmedian_edi(j)=nansum(data_edi(weights>0).*weights(weights>0));
        end
        G_edi=sortrows([temp,Gmedian_edi]);
        % then pes
        temp = BB_time(a2);[temp,ia,ic] = unique(temp);
        data_pes = VE(a2)./VdrivePes(a2);
        data_pes = data_pes(ia);
        data_pes(data_pes>prctile(data_pes,99))=NaN;
        data_pes(data_pes<prctile(data_pes,1))=NaN;
        temp(isnan(data_pes))=[]; data_pes(isnan(data_pes))=[];
        maxt = 1800;
        Gmedian_pes=NaN*temp;
        for j=1:length(temp)
            temp2 = abs(temp(j)-temp);
            weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
            Gmedian_pes(j)=nansum(data_pes(weights>0).*weights(weights>0));
        end
        G_pes=sortrows([temp,Gmedian_pes]);

        if ShowFigures
            figure(100);
            ax(1) = subplot(2,1,1); hold on;
            plot(G_edi(:,1),G_edi(:,2),'g-');
            ax(2) = subplot(2,1,2); hold on;
            plot(G_pes(:,1),G_pes(:,2),'g-');
        end

        if ShowExtraFigures
            figure(101); clf(figure(101));
            ax101(1) = subplot(2,1,1);
            plot(BB_time, VE./VdriveEdi,'kd');hold('on'); % plot everything
            plot(BB_time(a2),VE(a2)./VdriveEdi(a2),'rd'); % at least two breaths away from sleep
            plot(BB_time(a1),VE(a1)./VdriveEdi(a1),'gd'); % more than 4 breaths away from sleep
            plot(BB_time(Evt==1), VE(Evt==1)./VdriveEdi(Evt==1),'b.');
            refline(0,nanmedian(VE(a2)./VdriveEdi(a2)));
            ylabel('VE / VdriveEdi');

            win_nums = unique(win);
            win_nums(isnan(win_nums))=[];
            ax101(2) = subplot(2,1,2); linkaxes(ax101, 'x')
            for w = 1:length(win_nums)
                % use pt number n, and win number to get corresponding Flow and Edi
                [Flow, Edi] = getFlowEdiSignals(LocalSignals, pt, win_nums(w));
                [StTime, BB_Times] = getBBTimes(BreathDataTable, pt, win_nums(w));
                Samples = length(Flow);
                breathsinwin = sum(win==win_nums(w));
                % the signals are for the number of breathsinwin
                loc = find(win==win_nums(w),1,'first');
                EndTime = BB_time(loc+breathsinwin-1);
                deltaT = EndTime-StTime;
                %dt = deltaT/Samples;
                dt = 0.008; % force 1/125
                winlength = (length(Flow)*dt);
                %FlowTime = StTime:dt:StTime+((length(Flow)-1)*dt);
                FlowTime = StTime:dt:StTime+winlength;
                FlowTime(end)=[];

                plot(FlowTime, Flow./max(Flow), 'b'); hold on;
                plot(FlowTime, Edi./max(Edi), 'r');
                plot(StTime+BB_Times(:,1),0,'go');
                plot(StTime+BB_Times(:,2),0,'md');
                plot(StTime+BB_Times(:,3),0,'r.');
                refline(0,0);
                if mod(w,2)==0
                    plot(FlowTime, -1+FlowTime*0, 'k'); text(mean(FlowTime), -0.75, num2str(win_nums(w)));
                else
                    plot(FlowTime, -0.75+FlowTime*0, 'k'); text(mean(FlowTime), -1, num2str(win_nums(w)));
                end
            end
            ylabel('Flow and Edi');
        end

        % all breaths are then normalised, i.e. divided by their corresponding VE/Vdrive ratio
        % Edi first
        g_Edi_Adj = NaN(length(VE),1);
        Gt=interp1(G_edi(:,1),G_edi(:,2),BB_time);
            Gt(BB_time<G_edi(1,1))=G_edi(1,2);
            Gt(BB_time>G_edi(end,1))=G_edi(end,2);
        g_Edi_Adj = (VE./VdriveEdi)./Gt;

        % then Pes
        g_Pes_Adj = NaN(length(VE),1);
        clear Gt
        Gt=interp1(G_pes(:,1),G_pes(:,2),BB_time);
            Gt(BB_time<G_pes(1,1))=G_pes(1,2);
            Gt(BB_time>G_pes(end,1))=G_pes(end,2);
        g_Pes_Adj = (VE./VdrivePes)./Gt;

        % potentially remove major outliers, i.e. only use 5th to 95th centile
        %g_Edi_Adj(g_Edi_Adj>prctile(g_Edi_Adj,95))=NaN;
        %g_Pes_Adj(g_Pes_Adj>prctile(g_Pes_Adj,95))=NaN;

        % dividing by the 'best' ratio, we have normalised to 1.
        % values much >1 should occur infrequently, as this would suggest
        % VE exceeded Vdrive, as such, they should be capped (x>=1.5 == 1.5 )
        g_Edi_Adj(g_Edi_Adj>1.5)=1.5;
        g_Pes_Adj(g_Pes_Adj>1.5)=1.5;

        % a value of 1 would indicate good breathing, and <1 indicates FL.
        % to help model fitting, and not keeping misleading data, breaths
        % in wake that are clearly FL (x<0.7) are removed, but..
        % only if they are >threshold # of breaths away from sleep
        g_Edi_Adj(g_Edi_Adj<0.7 & (a1))=NaN; % far away from sleep, must be good
        g_Pes_Adj(g_Pes_Adj<0.7 & (a1))=NaN;
        %g_Edi_Adj(g_Edi_Adj<0.5 & (a2))=NaN; % more tolerant when closer to sleep
        %g_Pes_Adj(g_Pes_Adj<0.5 & (a2))=NaN;

        % now redo the figure for VE/Vdrive ratio with adjusted Edi and Pes
        % the wake data should be approximately 1
        if ShowFigures
            figure(100); clf(figure(100));
            ax100(1) = subplot(2,1,1);
            plot(BB_time, g_Edi_Adj,'kd');hold('on'); % plot everything
            plot(BB_time(a2),g_Edi_Adj(a2),'rd'); % at least two breaths away from sleep
            plot(BB_time(a1),g_Edi_Adj(a1),'gd'); % more breaths away from sleep
            refline(0,nanmedian(g_Edi_Adj(a2)));
            ylabel('VE / VdriveEdi');

            ax100(2) = subplot(2,1,2); linkaxes(ax100, 'x')
            plot(BB_time, g_Pes_Adj,'kd');hold('on'); % plot everything
            plot(BB_time(a2),g_Pes_Adj(a2),'rd'); % at least two breaths away from sleep
            plot(BB_time(a1),g_Pes_Adj(a1),'gd'); % more than 4 breaths away from sleep
            refline(0,nanmedian(g_Pes_Adj(a2)));
            ylabel('VE / VdrivePes');
        end

        else % not doing drive analysis, i.e. just using this for testing and/or training purposes
            %manually set these as NaN's
            g_Edi_Adj = NaN(length(BB_time),1);
            g_Pes_Adj = NaN(length(BB_time),1);

        end

        % Make a table of all breaths from this window
        VarNames = {'PT','Window','Ar','NotAr','A1','A2','Hypnog','Etype','Veup','VE','DriveEdi','DrivePes',...
            'g_Edi_Adj','g_Pes_Adj','BB_time','BB_Ttot'};
        PT = pt*ones(length(VE),1);
        PtDataWin = table(PT, win, Ar, notAr, a1, a2, hypnog, EType, Veup, VE, VdriveEdi, VdrivePes,...
            g_Edi_Adj, g_Pes_Adj, BB_time, BB_Ttot, 'VariableNames', VarNames);

        % add the Flow limited table to the current data table
        PtDataWinFL = [PtDataWin,FL];

        % add this window data with the FL data, to a growing combined table
        CombinedPtData = cat(1,CombinedPtData,PtDataWinFL);

    catch me
    end
end

%% ensure unique breaths in the combined table
[~,ia,~] = unique(CombinedPtData.BB_time, 'rows','stable');
dups=length(CombinedPtData.BB_time)-length(ia);
str = ['Removing ', num2str(dups), ' duplicates and sorting data']; disp(str);
PtData = CombinedPtData(ia,:);

%% tidy up the workspace, i.e. only need PtData from here onwards
clearvars -except PtData n_pts DataSet savestr_SVMPrep savestr_SVMRun
%clear LocalSignals

%% remove breaths with < threshold% flow below eupnoea
% where VE < 10% of Veup
% ToDo: change to 10% when testing is complete
LowFlow = PtData.VE < (PtData.Veup*0.1);
str = ['Removing ', num2str(nnz(LowFlow)), ' apneic breaths']; disp(str);
PtData(LowFlow,:)=[];

%% Remove extreme outliers: Part 1
% for e.g. VTi_VTe should be around 1, so values of 100 are clearly abnormal. 
% Indeed, values of 10 are also probably not normal. 
% Ideally, reasonable limits for each feature should be set.  

%mu = repmat(iqr(PtData{:,16:end}),size(PtData,1),1);
if 0 
%original simplistic method
centileupper = prctile(PtData{:,16:end},99.9);
centilelower = prctile(PtData{:,16:end},0.1);
%
%DM comment: in the matlab documentation, outlier detection using 
% quantiles is F1=Q1-1.5*IQR, and F2=Q3+1.5*IQR, with F1 > outliers > F2
IQR_ = iqr(PtData{:,16:end});
Q1_ = quantile(PtData{:,16:end}, 0.25);
Q3_ = quantile(PtData{:,16:end}, 0.75);
prctilefactor=1.5; %e.g. set to about 5 to be equivalent to Scotty's method
centileupper_ = Q3_ + prctilefactor*IQR_;
centilelower_ = Q1_ - prctilefactor*IQR_;
else %SS comment: SS revised version labels outliers e.g. if 2 fold above the 95-50 difference, e.g. median = 5, 95th centile = 10; outliers >15.
prctilethres=99; %e.g. 95
prctilefactor=3; %e.g. 2
centileupper = (prctile(PtData{:,16:end},prctilethres)-prctile(PtData{:,16:end},50))*prctilefactor + prctile(PtData{:,16:end},50);
centilelower = prctile(PtData{:,16:end},50) - (prctile(PtData{:,16:end},50)-prctile(PtData{:,16:end},100-prctilethres))*prctilefactor;
end
outliers = (PtData{:,16:end} > centileupper | ...
            PtData{:,16:end} < centilelower);
filler = false(size(PtData,1),15); % exclude the first 15 columns of PtData 
outliers = ([filler outliers]);
perBB = any(outliers,2);
str = [ num2str(nnz(outliers)),' feature outliers were found in ', num2str(nnz(perBB)), ' breaths']; disp(str);
%SS comment: I think better to simply set the values to be equal to the upper/lower boundary

%% Remove extreme outliers: Part 2
% this is a horribly inneficient method of covering the entire mat.
% logical indexing would be much better, but I can't see an easy way
% to do this in the table structure
tic
if 0
    for m=1:size(PtData,1)
        for n=1:size(PtData,2)
            if outliers(m,n)
                PtData{m,n}=NaN;
            end
        end
    end
else %SS version, replaces with 
    for n=1:length(centileupper)
        PtData{:,n+16-1}(PtData{:,n+16-1} > centileupper(n)) = centileupper(n);
        PtData{:,n+16-1}(PtData{:,n+16-1} < centilelower(n)) = centilelower(n);
    end     
end
toc
%histogram(PtData{:,21})
%probplot(PtData{:,19})
%normplot(PtData{:,21})

%% how many hypops and arousals do we have
% 5=Central 6=Obstructive 7=OHypopnea 15=Mixed 16=CentralHypopnea
Hypop = (PtData.Etype==6 | PtData.Etype==7);
Arous = (PtData.Etype==0 & PtData.Ar==1);
table(nnz(Hypop), nnz(Arous),'VariableNames',{'Hypop','Arousal'})

%% set up categorical thresholds for Mild, Moderate and Severe FL
% ToDo: correct for Drive data
% inclusive categorical thresholds, i.e.
Mild_FL = PtData.VE < (PtData.Veup*0.7);
Mod_FL = PtData.VE < (PtData.Veup*0.5);
Sev_FL = PtData.VE < (PtData.Veup*0.3);

% exclusive categorical thresholds, i.e.
% 50% < Mild_FL_e < 70% and  30% < Mod_FL_e < 50%  and  10% < Sev_FL_e < 30%
Mild_FL_e = (PtData.VE < (PtData.Veup*0.7) & PtData.VE > (PtData.Veup*0.5));
Mod_FL_e = (PtData.VE < (PtData.Veup*0.5) & PtData.VE > (PtData.Veup*0.3));
Sev_FL_e = (PtData.VE < (PtData.Veup*0.3) & PtData.VE > (PtData.Veup*0.1));
ALL_FL_e = zeros(length(Mild_FL_e),1);
ALL_FL_e(Mild_FL_e)=1;
ALL_FL_e(Mod_FL_e)=2;
ALL_FL_e(Sev_FL_e)=3;

table([nnz(Mild_FL);nnz(Mild_FL_e)],[nnz(Mod_FL);nnz(Mod_FL_e)],[nnz(Sev_FL);nnz(Sev_FL_e)],...
    'VariableNames',{'Mild','Mod','Sev'},'RowNames',{'Inclusive','Exclusive'})

table(nnz(Mild_FL), nnz(Arous),'VariableNames',{'Mild_FL','Arousal'})
table(nnz(Mild_FL), nnz(Hypop),'VariableNames',{'Mild_FL','Hypop'})
table(nnz(Mild_FL | Hypop), nnz(Mild_FL & Hypop), 'VariableNames',{'Mild_FL_OR_Hypop','Mild_FL_AND_Hypop'})

%% how many breaths of each type/category do we have per patient
BBperPt = NaN(n_pts,11);
for pt = 1:n_pts
    BBperPt(pt,1) = pt;
    BBperPt(pt,2) = size(PtData(PtData.PT==pt,1),1); %total
    BBperPt(pt,3) = size(PtData(PtData.PT==pt & PtData.Ar==1,1),1); % arousal
    BBperPt(pt,4) = size(PtData(PtData.PT==pt & PtData.A1==1,1),1); % A1
    BBperPt(pt,5) = size(PtData(PtData.PT==pt & PtData.A2==1,1),1); % A2
    BBperPt(pt,6) = size(PtData(PtData.PT==pt & PtData.Etype~=0,1),1); % event
    BBperPt(pt,7) = size(PtData(PtData.PT==pt & Sev_FL_e==1,1),1); % severe only breaths
    BBperPt(pt,8) = size(PtData(PtData.PT==pt & Mod_FL_e==1,1),1); % moderate only breaths
    BBperPt(pt,9) = size(PtData(PtData.PT==pt & Mild_FL_e==1,1),1); % mild only breaths
    BBperPt(pt,10) = size(PtData(PtData.PT==pt & Hypop,1),1); % hypops
    BBperPt(pt,11) = size(PtData(PtData.PT==pt & Arous,1),1); % arousals
    BBperPt(pt,12) = size(PtData(PtData.PT==pt & (Hypop | Mild_FL),1),1); % hypop or MildFL
end
varnames = {'PT','Total','Ar','A1','A2','Evts','SevereFL','ModFL','MildFL','Hypop','Arousal','HypMildFL'};
BBperPt_table = table(BBperPt(:,1),BBperPt(:,2),BBperPt(:,3),BBperPt(:,4),...
    BBperPt(:,5),BBperPt(:,6),BBperPt(:,7),BBperPt(:,8),BBperPt(:,9),BBperPt(:,10),BBperPt(:,11),BBperPt(:,12),...
    'variablenames',varnames);
clear varnames

table(min(BBperPt(:,2)), max(BBperPt(:,2)), 'VariableNames', {'Min','Max'})

figure(110); clf(figure(110));
ax110(1)=subplot(1,4,1:3);
plot(BBperPt(:,1),BBperPt(:,2), 'ks', 'markersize', 8); hold on;
refline(0,median(BBperPt(:,2)));
xlabel('PT #'); ylabel('# of Breaths');
ax110(2)=subplot(1,4,4);
boxplot(BBperPt(:,2),ones(BBperPt(end,1),1));
linkaxes(ax110, 'y')

figure(111); clf(figure(111));
ax111(1)=subplot(3,1,1);
y=[BBperPt(:,5) BBperPt(:,3)-BBperPt(:,5) BBperPt(:,2)-BBperPt(:,3);];
bar(y,'stacked');
legend('A2','Ar','Total');
ylabel('# of Ar and Total Breaths');
ax111(2)=subplot(3,1,2);
y=[BBperPt(:,7) BBperPt(:,8) BBperPt(:,9);];
bar(y,'stacked');
legend('Severe','Moderate','Mild');
ylabel('# of FL Breaths');
ax111(3)=subplot(3,1,3);
y=[BBperPt(:,12),BBperPt(:,11), BBperPt(:,2)];
bar(y);
legend('HypMildFL','Arousal','Total');
xlabel('PT #');
ylabel('# of Breaths');

%% make a backup of the data
PtDataBackup=PtData;

%% re-instate the original data from backup
PtData=PtDataBackup;


%% set up continuous input
switch DataSet
    case 'OralAppliance'
        Gtest = PtData{:,10}./PtData{:,9};
    case 'FlowDrive'
        Gtest = PtData{:,13};
end

criteriaR=1*Gtest;
criteriaR(isnan(Gtest))=NaN;
Yvariable=single(Gtest);


%% Weight the input breaths - for regression
% if a patient only has few breaths of a category being tested, then that
% will heavily bias weighted data, as such they must be controlled. The
% approach here is to set these to have no weight, and then remove them.
% In LDA, unequal sample sizes are acceptable. However,
% the sample size of the smallest group needs to exceed the number of
% predictor variables.  As a rule of thumb, the smallest sample size
% should be at least 20 for a few (4 or 5) predictors.
min_breaths = 20;
N=size(PtData,1);
weights=zeros(N,1);
pt_out = [];
for pt_in=1:n_pts
    I1=find(pt_in==PtData{:,1});
    if length(I1)<=min_breaths % this pt doesn't have enough breaths
        weights(I1)=0; % so set the weights to zero, and
        pt_out = cat(1,pt_out,pt_in); % keep a list of pts to be excluded
    else
        weights(I1)=(1/length(I1)); % otherwise, set the weights
    end
end
str = ['After weighting the data, ', num2str(sum(weights)), ' patients remain']; disp(str)
% remove the data
exclude=(weights==0);
str = ['As such, ', num2str(sum(exclude)), ' breaths were removed during weighting']; disp(str)
PtData(exclude,:)=[];
% Amatrix(exclude,:)=[]; % used for svm
criteriaR(exclude)=[];
Yvariable(exclude)=[];
Gtest(exclude)=[];
weights(exclude)=[];
% make a new list of pts to process, used for leave one out processing
pts = [1:1:n_pts];
pts(pt_out)=[];
n_pts=round(sum(weights));
% end of dwayne's weighting method

if 0
% give all breaths a weight of one
pts = [1:1:n_pts];
weights = ones(size(PtData,1),1);
end

%% report on how many NonFL and FL breaths from each pt 
% by Yvariable >1 or <1 definition

SumData = NaN(length(pts), 4);
for n=1:length(pts)
    SumData(n,1) = pts(n);
    SumData(n,2) = nnz(table2array(PtData(:,1))==pts(n) & criteriaR>=1);
    SumData(n,3) = nnz(table2array(PtData(:,1))==pts(n) & criteriaR<1);
    SumData(n,4) = nnz(table2array(PtData(:,1))==pts(n));
end
varnames = {'PT','NonFL','FL','Total'};
FL_BBperPt = table(SumData(:,1),SumData(:,2),SumData(:,3),SumData(:,4),...
    'variablenames',varnames)
clear varnames

%% set up the features list
ftrset = 8;
switch ftrset
    case 0
        % for the initial testing, we used a skip subset of all data
        rangei=18:5:115;
    case 1 % run1, then we manaully selected some features
        rangei = [18, 22, 30, 43, 46, 49, 50, 66, 96];
        % we forced the svm to select three ftrs
        % it found index 9, 4, 6, ie. 96, 43, and 49
        % MostPromPeak, SSArea, SkewDistInsp
        % the performance was saved as run1_performance
    case 2 % run2, again we manually set features
        rangei = [18,30,96];
        % we forced the svm to select three ftrs
        % as such features out was = features in (but the order may have changed)
        % MostPromPeak, InspFlat8020, TiTtot
        % the performance was saved as run2_performance
    case 3 % run3
        % we note that feature 49 (in run1) was pretty rubbish independently,
        % so this time we try again, but with those that are good independently
        rangei = [96, 40, 107, 22, 37, 75, 43, 48, 75];
        % we forced the svm to select three ftrs
        % it found index 2, 4, 7,which was
        % MIF50, MIF_PIF, and SSArea
        % the performance was saved as run3_performance
    case 4 % run4, now with all the features
        rangei = [16:1:114];
        Tbl_rangei = PtData(:,rangei);
        % we forced the svm to select three ftrs
        % it found MIF50, SkewDataInsp, InspFlutPow_Vpeak2
        % the performance was saved as run4_performance
    case 5 % run5, we manually set five features
        %rangei = [40,96,22,43,53];
        rangei = [96,22,43,53,40]; % in reverse order, just for testing
        % we forced the svm to select three ftrs
        % it found 1, 5, and 3, (on 2 occasions, it found 1,5, and 4)
        % MIF50, SkewDataInsp, MIF_PIF, (occasionally SS_area)
        % the performance was not saved

        % this was a test run, to check output was collected from each pt
        % so that we can do histogram/distributions on # of times each feature
        % is selected, and also look at per pt sens, spec, acc, ppv, npv, etc...
    case 6
        % the best 20 from previous svm runs, in no particular order,
        % this is actually 25 ftrs... so were did this list come from?
        %          1, 2, 3, 4, 5, 6, 7 ,8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22, 23, 24, 25
        rangei = [18,23,25,26,27,35,38,40,51,53,55,56,71,74,76,77,78,79,82,86,91,95,108,114,115];
        % eqvlnt [1, 6, 8, 9, 10,18,21,23,34,36,38,39,54,57,59,60,61,62,65,69,74,78, 91, 97, 98]
        % should [1, 4, 5, 6, 7,  8, 9,10,16,23,36,38,39,59,60,61,62,63,65,69,74,91, 92, 97, 98]

    case 7
        % used in regression testing
        % the "best" features from the binary classifier
        rangei = [];
    case 8 %SS test
        % used in regression testing
        rangei = 16-1+[23 66 37 93 80 7 1 31 61 75 69 74 91 24];    
    otherwise
        rangei = [];
        disp('Please set features list');
end


%% compile training data based on features list
N=size(PtData,1); % do again, as weighting above may have changed the number
Amatrix = PtData{:,rangei};
% set up the ignored and forced variables
tempignorevars=[];
forcedvars=[];

%% set up prior learning and prior performance (if exists)
% this is if it has been run through all patients, but we are now extending
% the maxNfeaturesJ variable.
% note: for consistency, must keep the input range (features) the same

rerun = 0; % set to zero to strat new run, or one if rerunning (then manually add prior results)
switch rerun
    case 0
        prior_Ilist = [];
        prior_perf1 = [];
        prior_perf2 = [];
    case 1
        prior_Ilist = Ilist_;
        prior_perf1 = PerPTPerformance1;
        prior_perf2 = PerPTPerformance2;
    otherwise
        prior_Ilist = [];
        prior_perf1=[];
        prior_perf2=[];
end

if 0 % re-instate these variables, useful during debugging
    Ilist_ = prior_learning;
    PerPTPerformance1 = prior_perf1;
    PerPTPerformance2 = prior_perf2;
end
% Other conditions that may arise, and some thoughts on solutions
% (1) resuming processing, after not completing all pts
% (perhaps manually stopped or interrupted). In this condition, both the
% input range (features) and output features list length must be unchanged.
%  solution:
%   load prior learning, then change the main pt for loop so that
%   "pt_in=pts" is the remaining pts.
%   also need to change "count" to reflect actual number through processng
% (2) resuming processing, after (a) not completing all pts, and (b) now
% extending the maxNfeaturesJ variable. In this condition, the input range
% (features) must be the same
%  solution:
%   ? - a combination of both, but I suspoect that the time taken to
%   write the code to handle this scenario will be longer than the time to
%   run the analysis (unless we start doing lots of big n analysis)


%% >>>>>>>>>>>>>>>>>>>   SAVE/LOAD   <<<<<<<<<<<<<<<<<<<<<<<<<
save(savestr_SVMPrep);
load(savestr_SVMPrep);
if 0
    close all
    clear
    LocalDataDir = 'C:\PSG_Data\QAO';
    LocalData = [LocalDataDir,'\FS_2_noSignals.mat'];
    savestr_SVMPrep = [LocalDataDir, '\SA_2_Reg_Prep.mat'];
    savestr_SVMRun = [LocalDataDir, '\SA_2_Reg_Run.mat']; 
end

%% SVM forward select run - regression
savestr_SVMRun = [LocalDataDir, '\SA_8R_Run.mat'];


rerun = 0;

maxNfeaturesJ=20;
Ilist_= NaN(length(pts),maxNfeaturesJ); % the list of features for each pt
%PredTAll_=NaN(N,1); % this is the predicted output, made by concatenating the output from each pt
PerPTPerformance_Train{length(pts)}=[]; % this is the performance in training data
PerPTPerformance_Test{length(pts)}=[]; % this is the performance in the left out pt
PerPTModel_Train{length(pts)}=[]; % this is the training model output for each pt
PerPTModel_Test{length(pts)}=[]; % this is the test model output for each pt
count=0;
% pt_in = 1;
for pt_in=pts
    count=count+1; disp(' ');
    str=(['Processing patient ID ', num2str(pt_in), '; number ', num2str(count),' of ', num2str(length(pts))]); disp(str);
    % set up the cross validation range
    rangevalidate=pt_in==PtData{:,1};
    
    % run the forward selection process, and return a list of features and
    % performance as each feature is added

    if rerun % this option includes prior learning
    [Ilist,PtPerf1,PtPerf2,SVMModel1,SVMModel2]=svmforwardselect_perpatient_regression(Amatrix,criteriaR,...
        tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,forcedvars,weights,...
        prior_Ilist(count,:), prior_perf1{count}, prior_perf2{count});
    else % no prior learning
        [Ilist,PtPerf1,PtPerf2,SVMModel1,SVMModel2]=svmforwardselect_perpatient_regression(Amatrix,criteriaR,...
        tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,forcedvars,weights,...
        [],[],[]);
    end
    % keep a record of the output i.e.
    % the feature list
    Ilist_(count,:)=Ilist;
    % the predT for each pt, using their model, built from their Ilist
    %PredTAll_(rangevalidate,1)=PredTAll(rangevalidate);
    % the history of performance in the training data
    PerPTPerformance_Train{count}=PtPerf1;
    % the history of performance in the left out data
    PerPTPerformance_Test{count}=PtPerf2;
    % the model
    PerPTModel_Train{count} = SVMModel1;
    PerPTModel_Test{count} = SVMModel2;
    % and save it
%     if exist(savestr_SVMRun,'file')==2
%             save(savestr_SVMRun,'PtData','PtDataBackup',...
%                 'Ilist_','PerPTPerformance_Train','PerPTPerformance_Test',...
%                 'PerPTModel_Train', 'PerPTModel_Test',...
%                 'count','pt_in','pts','-append');
%         else
%             save(savestr_SVMRun,'PtData','PtDataBackup',...
%                 'Ilist_','PerPTPerformance_Train','PerPTPerformance_Test',...
%                 'PerPTModel_Train', 'PerPTModel_Test',...
%                 'count','pt_in','pts');
%     end
end

%% save and/or load
save(savestr_SVMRun, '-v7');

if 0
    clear
    close all
    LocalDataDir = 'C:\PSG_Data\QAO';
    savestr_SVMRun = [LocalDataDir, '\FS_4_noSignals.mat'];
    load(savestr_SVMRun);
end


%% what features appear most frequently
% note: this doesn't sort by order of appearance in the list
centers = 1:1:100;
Ilist_linear = Ilist_( : );
counts = hist(Ilist_linear,centers);
%pcts = 100 * counts / sum(counts);
figure(209); clf(figure(209));
bar(centers,counts);
fig = gcf;
fig.Color = [1 1 1]; % set background colour to white
fig.Units = 'inches';
fig.Position = [1 1 12 6];
ax = gca;
ax.TickDir = 'out';
ax.FontSize=7;
ax.XTick=[0:2:length(centers)];
ax.YTick=[0:1:30];
ylabel('Frequency');
xlabel('Feature #');
title('Histogram of features selected (98 c 20)');
box off

%savefig(fig, 'c:\Temp\HistogramFtrs_98c20.fig');

%% Table of features selected
for j=1:6%size(Ilist_,2)
    disp(' ');
    str_title=['Fwd Select #',num2str(j)]; disp(str_title);
    ftrs_found = unique(Ilist_(:,j));
    for k=1:size(ftrs_found,1)
        str_data=['Ftr:', num2str(ftrs_found(k)), ' QTY: ',num2str(nnz(Ilist_(:,j)==ftrs_found(k)))];
        disp(str_data);
    end
end

%% Run on all data, with Ilist, and create one model for entire dataset
% get the top 20 ftrs, and use only these
[i,a]=sort(counts,'descend'); [i(1:20)',a(1:20)']

Amatrix_reduced = Tbl_rangei{:,a(1:20)};    % mat format
Tbl_rangei_reduced =  Tbl_rangei(:,a(1:20));% table format

[Ilist,PtPerf_train,PtPerf_test,SVMModel_train,SVMModel_test]=svmforwardselect_perpatient_regression...
    (Amatrix_reduced,criteriaR,tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,...
    forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights,[],[],[]);

% what does the AWP look like as we add features
PtPerf = PtPerf_train;
awp_per_ft = (PtPerf(:,3)+PtPerf(:,4)) / 2;
figure(200); clf(figure(200));
stairs(awp_per_ft);

% then make a prediction using the current model, for all the data
PredTsvmr = predict(SVMModel_test,Amatrix_reduced);

figure()
scatter(PredTsvmr, criteriaR);

performance = PredictiveValue(criteriaR,PredTAll_,Yvariable);
run6_performance = performance;


%% Performance analysis
%use the training or testing data, for perfromance analysis below
TrainTestData = 1; % zero for Training data, one for left out patient data
switch TrainTestData
    case 0
        PerPTPerformance = PerPTPerformance_Train;
    case 1
        PerPTPerformance = PerPTPerformance_Test;
        %PerPTPerformance = PerPTPerformance2;
end
% remove any empty cells
PerPTPerformance(cellfun('isempty',PerPTPerformance)) = [];

%% per pt performance (with either training or testing data from above)
% this makes a three dimensional matrix, with each column being the
% performance metric, each row being a patient, and each sheet being the
% results at each feature added
PTPerfData=NaN(length(PerPTPerformance),size(PerPTPerformance{1},2));
for ftr=1:size(PerPTPerformance{1},2)
    for pt=1:length(PerPTPerformance)
        if isempty(PerPTPerformance{1,pt})
            continue
        end
        PTPerfData(pt,ftr) = PerPTPerformance{1,pt}{1,ftr}; %PPV
    end
end

figure(201); clf(figure(201));
stairs(mean(PTPerfData)); hold on;
%refline(0,median(table2array(PTPerfTable(:,5))));
suptitle('Average weighted performance vs number of features');
title('AWP = mean (each pt ((sens+spec)/2) for ftr n)');
ylabel('R');
xlabel('Ftrs');
TidyPlot();
if 0
x = 4; y = 0.875;[m,i] = max(mean(awp_each));
txt = (['Max AWP = ', num2str(m), ', at Feature ', num2str(i)]);
text(x,y,txt);

x = 11; y = 0.84;
txt = ([{'Ftr #    AWP'},...
    {['1       ',num2str(mean(awp_each(:,1)))]},...
    {['2       ',num2str(mean(awp_each(:,2)))]},...
    {['3       ',num2str(mean(awp_each(:,3)))]},...
    {['4       ',num2str(mean(awp_each(:,4)))]},...
    {['5       ',num2str(mean(awp_each(:,5)))]},...
    {['6       ',num2str(mean(awp_each(:,6)))]},...
    {['7       ',num2str(mean(awp_each(:,7)))]},...
    {['8       ',num2str(mean(awp_each(:,8)))]},...
    {['9       ',num2str(mean(awp_each(:,9)))]},...
    {['10     ',num2str(mean(awp_each(:,10)))]},...
    {['11     ',num2str(mean(awp_each(:,11)))]},...
    {['12     ',num2str(mean(awp_each(:,12)))]},...
    {['13     ',num2str(mean(awp_each(:,13)))]},...
    {['14     ',num2str(mean(awp_each(:,14)))]},...
    {['15     ',num2str(mean(awp_each(:,15)))]},...
    {['16     ',num2str(mean(awp_each(:,16)))]},...
    {['17     ',num2str(mean(awp_each(:,17)))]},...
    {['18     ',num2str(mean(awp_each(:,18)))]},...
    {['19     ',num2str(mean(awp_each(:,19)))]},...
    {['20     ',num2str(mean(awp_each(:,20)))]}]);
text(x,y,txt);
end

if 0
figure(202); clf(figure(202));
%ax201(2)=subplot(3,1,2);
stairs(table2array(PTPerfTable(:,[6,8])));
legend('Pos','Neg');
ylabel('True Neg/Pos');
xlabel('Pt #');
figure(203); clf(figure(203));
%ax201(3)=subplot(3,1,3);
stairs(table2array(PTPerfTable(:,[7,9])));
legend('Pos','Neg');
ylabel('False Neg/Pos');
xlabel('Pt #');
%linkaxes(ax201, 'x');
end

%% how does performance change with number of features (ftrs out)
PerfPPV=NaN(length(PerPTPerformance),length(Ilist));
PerfNPV=NaN(length(PerPTPerformance),length(Ilist));
PerfSens=NaN(length(PerPTPerformance),length(Ilist));
PerfSpec=NaN(length(PerPTPerformance),length(Ilist));
PerfAcc=NaN(length(PerPTPerformance),length(Ilist));
for pt=1:length(PerPTPerformance)
    if isempty(PerPTPerformance{1,pt})
        continue
    end
    st = ((pt*length(Ilist))-(length(Ilist)))+1;
	en = (st+length(Ilist))-1;
    for n_=1:length(Ilist)
    PerfPPV(st:en,n_) = PerPTPerformance{1,pt}(n_,1); %PPV
    end
    for n_=1:length(Ilist)
    PerfNPV(st:en,n_) = PerPTPerformance{1,pt}(n_,2); %NPV
    end
    for n_=1:length(Ilist)
    PerfSens(st:en,n_) = PerPTPerformance{1,pt}(n_,3); %Sens
    end
    for n_=1:length(Ilist)
    PerfSpec(st:en,n_) = PerPTPerformance{1,pt}(n_,4); %Spec
    end
    for n_=1:length(Ilist)
    PerfAcc(st:en,n_) = PerPTPerformance{1,pt}(n_,5); %Acc
    end
end

figure(204); clf(figure(204));
%ax202(1)=subplot(5,1,1);
boxplot(PerfPPV, 'Whisker',inf);
%title('Performance per feature');
ylabel('PPV');
xlabel('# of Features');
figure(205); clf(figure(205));
%ax202(2)=subplot(5,1,2);
boxplot(PerfNPV, 'Whisker',inf);
ylabel('NPV');
xlabel('# of Features');
figure(206); clf(figure(206));
%ax202(3)=subplot(5,1,3);
boxplot(PerfSens, 'Whisker',inf);
ylabel('Sens');
xlabel('# of Features');
figure(207); clf(figure(207));
%ax202(4)=subplot(5,1,4);
boxplot(PerfSpec, 'Whisker',inf);
ylabel('Spec');
xlabel('# of Features');
figure(208); clf(figure(208));
%ax202(5)=subplot(5,1,5);
boxplot(PerfAcc, 'Whisker',inf);
ylabel('Acc');
xlabel('# of Features');
%linkaxes(ax202, 'y');



%% what does the data look like?
% let's visualise the data, at least the first three variables
%

%% try the svm 3d plot function

% start with the most frequently found features
%Ilist_mode = mode(Ilist_);
%Ilist_mode = Ilist_mode(1,[1:3]);
Ilist_mode = [23, 61, 66]; % force
Ilist_mode = [80, 1, 7]; % force
SVMModel3 = fitcsvm(Amatrix(:,Ilist_mode),criteriaR,'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelfunction,'KernelScale',sigma, 'Weights', weights);
PredT3 = predict(SVMModel3,Amatrix(:,Ilist_mode));


SVMModel4 = svmtrain(Amatrix(:,Ilist_mode),criteriaR,'method',svpmethodstr,...
        'kernel_function',kernelfunction,'rbf_sigma',sigma);
PredT4 = svmclassify(SVMModel4,Amatrix(:,Ilist_mode));

svm_3d_matlab_vis_mod(SVMModel3,Amatrix(:,Ilist_mode),criteriaR,0);
svm_3d_matlab_vis_mod(SVMModel4,Amatrix(:,Ilist_mode),criteriaR,0);

if 0 % original limits
    xlim([0 15]);
    ylim([-1 1]);
    zlim([0 0.06]);
end
% modified limits
if 0
xlim([0 1.2]);
ylim([-1 1]);
zlim([0 0.01]);

xlim([0 1]);
ylim([0 2]);
zlim([0 0.01]);


end
xlabel(['Ftr #',num2str(Ilist_mode(1))]);
ylabel(['Ftr #',num2str(Ilist_mode(2))]);
zlabel(['Ftr #',num2str(Ilist_mode(3))]);


%fig=gcf; savefig(fig, 'c:\Temp\3dscatter_withisosurface.fig');


%% try a standard 3d scatter plot
PredTAlldata = predict(SVMModel3,Amatrix(:,Ilist_mode));

Ilist_mode_perft = mode(Ilist_);
Ilist_mode = Ilist_mode_perft(1,[1:3]);
C = NaN(length(criteriaR),3);
err = NaN(length(criteriaR),1);
% this is a terribly inefficient method to set colours for criteria
for k = 1:length(criteriaR)
    if criteriaR(k)==1
        C(k,1:3)=[1,0,0];
    else
        C(k,:)=[0,0,1];
    end
    % also set a flag if we failed to predict this one
    err(k) = (criteriaR(k)~=PredTAlldata(k));
end

% restrict the size down a little, a bit more manageable to view
subset = 1:5:size(Amatrix,1);
C = C(subset,:); err = err(subset,:);
S = ones(length(C),1)*4;
X = Amatrix(subset,Ilist_mode(1));
Y = Amatrix(subset,Ilist_mode(2));
Z = Amatrix(subset,Ilist_mode(3));

figure(210); clf(figure(210));
%h=scatter3sph(X,Y,Z,'color',C,'size',S); % this method is very slow
h = scatter3(X, Y, Z, S, C, 'filled'); hold on;
xlabel(['Ftr #',num2str(Ilist_mode(1))]);
ylabel(['Ftr #',num2str(Ilist_mode(2))]);
zlabel(['Ftr #',num2str(Ilist_mode(3))]);
% these limits will almost certainly need changing for different Ilist
%xlim([0 1.2]);
%zlim([0 0.0025]);

% mark the misclassifications
scatter3(X(logical(err)), Y(logical(err)), Z(logical(err)));
box on

%figure(210); fig=gcf;
%savefig(fig, 'c:\Temp\3dscatter_withmisclass.fig');

%%
ctrs = 1:1:20;
IList_counts_overall = hist(Ilist_( : ), ctrs);

AllFoundFtrs = Ilist_( : );
[counts,vals]= hist(AllFoundFtrs);

%%
if 0
rangevalidate=ones(N,1); %all data, all patients
[IlistAll,PredTAlldata,tempXX2]=svmforwardselect_perpatient_regression(Amatrix,criteriaR,...
    tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights);
performanceAll = PredictiveValue(criteriaR,PredTAlldata,Yvariable);
end


%%  multiclass classification
% regression model
SVMModel_reg = fitrsvm(Amatrix(:,[23,36,65]),criteriaR,'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction','polynomial','KernelScale',sigma, 'Weights', weights);%'OptimizeHyperparameters','auto'

PredTsvm_reg = predict(SVMModel_reg,Amatrix);
figure(211); clf(figure(211)); plot(PredTsvm_reg);
[thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(criteriaR,PredTsvm_reg,1); %

%% continuous prediction
% set Yvariable as the continuous variable

%% experimenting...
% naive bayes
NB_data = PtDataBackup(:,[96, 40, 107, 22, 37, 75, 43, 48, 75]);
Mdl_nb = fitcnb(NB_data,ALL_FL_e, 'ClassNames',{'0','1','2','3'});
PredT_nb = predict(Mdl_nb,NB_data);
PredT_nb_mat = str2num(cell2mat(PredT_nb));
PredT_nb_mat_resp = [PredT_nb_mat,ALL_FL_e];
c = confusionmat(PredT_nb_mat_resp(:,1),PredT_nb_mat_resp(:,2));


%% Distibution within features

% some more boxplots
figure(212); clf(figure(212));
group = [ones(nnz(Yvariable==0),1)*1;ones(nnz(Yvariable==1),1)*2];
subplot(1,3,1);
data = [Amatrix(Yvariable==0,Ilist(1)); Amatrix(Yvariable==1,Ilist(1))];
boxplot(data,group);
xlabel('PkProm');

subplot(1,3,2);
data = [Amatrix(Yvariable==0,Ilist(2)); Amatrix(Yvariable==1,Ilist(2))];
boxplot(data,group);
xlabel('SSArea');

subplot(1,3,3);
data = [Amatrix(Yvariable==0,Ilist(3)); Amatrix(Yvariable==1,Ilist(3))];
boxplot(data,group);
xlabel('Skew');

% some more scatterhists
figure(213); clf(figure(213)); % 36, 74
scatterhist(Amatrix(:,79), Amatrix(:,74), 'Color','br', 'Group', Yvariable, 'Marker','..','MarkerSize',[4,4]); hold on;
xlabel('PkProm'); ylabel('Tpeak1Ti');
%axis square; r(1) = refline(1,0); r(1).Color='y'; r(1).LineWidth = 1.5;
legend off;
fig = gcf;
fig.Color = [1 1 1]; % set background colour to white
fig.Units = 'inches';
fig.Position = [1 1 9 8];
ax = gca;
ax.TickDir = 'out';
title('Scatter plot and Histogram of Ali''s ftrs');
box off
legend('NonFL','FL');

%savefig(fig, 'c:\Temp\ScatterHistFtrs_PkPromTpeak1Ti.fig');

%% SS test
if 1
    %loaded FS_FD_1_noSignals, then ran to line 658. Data also in 'E:\Work\Projects SS\QAO\SSworkspace.mat'
    
    %% SVM regress
    
    %weights and performance are not based on individual patients, but rather individual breaths
   
    %setup
    Gtest = 
    criteriaR=1*Gtest; %VE/Vdrive
    %criteriaR=PtData.VE./PtData.Veup; %VEnormalized
                 %[24 66 1 37 80 7 23 31 61 75 69 74 91 93 25 40 58 26 62 8 13]
    rangei = 16-1+[24 66 1 37 80 7 23 31 61 75 69 74 91 93 25 40 58 26 62 8 13]; %some features that are usually reliable  1 37 80 7 23 31 61 75 69 
    Amatrix = PtData{:,rangei};
    
    %remove NaN (AIF)
    I = isnan(criteriaR);
    Amatrix(I,:)=[];
    criteriaR(I,:)=[];
    weights_=weights;
    weights_(I,:)=[];
    
    %remove features with lots of NaNs or Infs
    temp=sum(isnan(Amatrix)|(Amatrix==Inf)|(Amatrix==-Inf))/size(Amatrix,1);
    Amatrix(:,temp>0.01)=[];
    
    %Downsample for speedup
    downsamplefactor=3; %downsample 10 takes just 2s [similar for 1-20 features],5->11s, 4->15s, 3->25s, none->3-5 min
    I=ceil(rand(1,1)*downsamplefactor):downsamplefactor:length(criteriaR); %downsample 5 takes
%     if 1 %too few FL breaths.
%     thres=0.7;
%     I2=find(criteriaR<0.7);
%     end
    %Run SVM regress
    tic
    svpmethodstr = 'SMO';
    kernelf = 'rbf';%'polynomial'
    sigma=2;
    criteriaR_ = criteriaR;
    SVMModel_reg = fitrsvm(Amatrix(I,:),criteriaR_(I),'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelf,'KernelScale',sigma,'Weights',weights_(I));%'OptimizeHyperparameters','auto' ,'Weights',weights
    toc
    
    tic
    PredTsvm_reg = predict(SVMModel_reg,Amatrix(:,:));
    toc
    
    %figure(); 
    %plot(criteriaR_(I),PredTsvm_reg(I)-criteriaR_(I),'.','markersize',1);
    
    if 0
    %run equivalent SVM classifier
    w=1+0*criteriaR_(I);
    ratio=sum((criteriaR_(I)<thres))/sum((criteriaR_(I)>=thres));
    w(criteriaR_(I)<thres)=w(criteriaR_(I)<thres)/ratio;
    %ratio=sum(w(criteriaR_(I)<thres))/sum(w(criteriaR_(I)>=thres)); %test ratio=1
    SVMModel_c = fitcsvm(Amatrix(I,:),(criteriaR_(I))<thres,'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelf,'KernelScale',sigma,'Weights',w);%'OptimizeHyperparameters','auto' ,'Weights',weights
    PredTsvm_c = predict(SVMModel_c,Amatrix(:,:));
    performance_c = PredictiveValue(criteriaR_<thres,PredTsvm_c,PredTsvm_c);
    end
    
    %plot results: good
    figure(211); clf(figure(211));
    subplot(1,3,1);
    %plot(criteriaR_,PredTsvm_reg,'.','markersize',1);
    % [Rvalue,Pvalue,Slope,Int]=plotregressionwithSEM(PredTsvm_reg,criteriaR_);
    % % the Fn above may be a newer version, with "Int" returned
    [Rvalue,Pvalue,Slope]=plotregressionwithSEM(PredTsvm_reg,criteriaR_);
    h=get(gca, 'Children'); set(h(1),'markersize',0.5,'color',[0.2 0.6 0.9])
    
    subplot(1,3,2); plot(criteriaR_,PredTsvm_reg-criteriaR_,'.','markersize',1); hold('on'); box('off');
    plot([0 max(criteriaR_)],nanmean(PredTsvm_reg-criteriaR_)*[1 1],'k-'); 
    plot([0 max(criteriaR_)],nanmean(PredTsvm_reg-criteriaR_)+1.96*nanstd(PredTsvm_reg-criteriaR_)*[1 1],'k--'); 
    plot([0 max(criteriaR_)],nanmean(PredTsvm_reg-criteriaR_)-1.96*nanstd(PredTsvm_reg-criteriaR_)*[1 1],'k--'); 
    
    subplot(1,3,3);
    [thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(criteriaR_<thres,PredTsvm_reg,1);
    performance = PredictiveValue(criteriaR_<thres,PredTsvm_reg<thresX,PredTsvm_reg);
    title(['AUC: ' num2str(AUC,2) ' Acc*: ' num2str(mean(sensspec),2)])
    
    subplot(1,3,1);
    dx=0.1;
    xbins=[0 0.15:dx:1.05 1.5];
    title(['Rvalue: ' num2str(Rvalue,2)])
    %plot binned data
    for i=1:length(xbins)-1
        I=PredTsvm_reg>xbins(i)&PredTsvm_reg<xbins(i+1);
        medianX(i)=prctile(PredTsvm_reg(I),50);
        medianY(i)=prctile(criteriaR_(I),50);
        upperIQRY(i)=prctile(criteriaR_(I),75);
        lowerIQRY(i)=prctile(criteriaR_(I),25);
        upperIQRY2(i)=prctile(criteriaR_(I),90);
        lowerIQRY2(i)=prctile(criteriaR_(I),10);
    end
    plotbarwitherrorsmedian(medianY,medianX,upperIQRY2,lowerIQRY2,zeros(length(medianY),3),NaN,0.005,0.01)
    plotbarwitherrorsmedian(medianY,medianX,upperIQRY,lowerIQRY,zeros(length(medianY),3),NaN,dx/2,0.01)

    %create classes from regression SVM results
    
    %classcutoffs = [0.7 0.4];
    classcutoffs = [0.9 0.7 0.5]; %mild moderate severe
    Nclasses=length(classcutoffs)+1; 
    
    g = NaN*criteriaR_;
    ghat = NaN*PredTsvm_reg;
    g(criteriaR_>classcutoffs(1))=1;
    ghat(PredTsvm_reg>classcutoffs(1))=1;
    for i=2:Nclasses
        g(criteriaR_<classcutoffs(i-1))=i;
        ghat(PredTsvm_reg<classcutoffs(i-1))=i;
    end 
    
    C = confusionmat(g,ghat)
    %sumactual=sum(C')';
    sumactual=sum(C,2);
    sumestimated=sum(C);
    %C_Factual=C./sum(C')'
    C_Factual=C./sum(C,2)*100
    C_Festimated=C./sum(C)*100;
    %plotconfusion(g,ghat) %needs neuralnet tools, and a lot of memory too.
    % I get an "Out of memory" error when attempting plotconfusion with the
    % current input arguments
    
    if 0
    %try multiclass--not working seems biased to scoring all as the
    %majority class despite weights
    for i=1:Nclasses
        tic
        indx=(g(I)==i);
        w=1+0*indx;
        currentratio=sum(indx)/length(indx); %fudge
        currentratio=currentratio^2;
        w(indx)=w(indx)/currentratio;
        [SVMModel_reg_{i}] = fitcsvm(Amatrix(I,:),g(I)==i,'Standardize',true,'Solver',svpmethodstr,...
            'KernelFunction',kernelf,'KernelScale',sigma,'Weights',w);%'OptimizeHyperparameters','auto' ,'Weights',weights
        toc
        [PredTsvm_reg_(:,i),score] = predict(SVMModel_reg_{i},Amatrix(:,:));
        Scores(:,i) = score(:,2);
        toc
    end
    [~,ghat_class] = max(Scores,[],2);
    C_class = confusionmat(g,ghat_class)
    end
    
%%
end



