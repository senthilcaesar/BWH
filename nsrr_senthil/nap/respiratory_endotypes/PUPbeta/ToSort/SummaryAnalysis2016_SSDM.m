%% startup
close all
clear global
clear
clc

%% specify a user and individual settings
user = 'dwayne'; %'Scotty'; %
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
        DataSet = 'OralAppliance';  %'FlowDrive'
        
        LocalDataDir = '';  %set
        LocalData = [LocalDataDir,'\FS_OA_2_noSignals.mat']; %set
        savestr_SVM = [LocalDataDir, '\SVMRun_.mat']; %set
        
    case 'dwayne'
        % path to code
        CodeDirectory = 'C:\Users\uqdmann\Dropbox\QAO\Code';
        addpath(CodeDirectory);
        cd(CodeDirectory);
        
        % path to data
        LocalDataDir = 'C:\PSG_Data\QAO';
        
        if 0 % combining two feature sets
            % load one set of data
            %LocalData = [LocalDataDir,'\FS_FD_6.mat']; % ttran timing
            LocalData = [LocalDataDir,'\SA_FD_10.mat']; % ttran timing
            load(LocalData); %
            PtData_Ttran = PtDataBackup(:,1:17);
            FtrsData_Ttran = PtDataBackup(:,18:end);
            FtrsData_Ttran(:,ftrsToExclude) = [];
            
            clearvars -except LocalDataDir PtData_Ttran FtrsData_Ttran
            
            % load the ftrs list
            load('C:\PSG_Data\QAO\Ftrs_OA.mat');
            %load('C:\PSG_Data\QAO\Ftrs_FD.mat');
            
            FtrsData_Ttran = FtrsData_Ttran(:,Ttran_Ftrs_OA);
            
            % load the second set of data
            LocalData = [LocalDataDir,'\FS_FD_5.mat']; % original timing
            %LocalData = [LocalDataDir,'\SA_FD_9.mat']; % original timing
            
            load(LocalData);
            
            PtData_Orig = PtDataBackup(:,1:17);
            FtrsData_Orig = PtDataBackup(:,18:end);
            FtrsData_Orig(:,ftrsToExclude) = [];
            FtrsData_Orig = FtrsData_Orig(:,Original_Ftrs_OA);
            
            clearvars -except LocalDataDir ...
                PtData_Ttran FtrsData_Ttran ...
                PtData_Orig FtrsData_Orig
            
            if ~(isequaln(PtData_Orig, PtData_Ttran)) % check for equality, treating nan's as equals
                disp('PtData is not equivalent');
                keyboard;
            else
                PtData = PtData_Orig;
                
                % append ftr number and _O to original features
                FeatureNames_OriginalList = FtrsData_Orig.Properties.VariableNames';
                for n = 1:length(FeatureNames_OriginalList)
                    FeatureNames_OriginalList{n}=[FeatureNames_OriginalList{n,1},'_',num2str(n),'_O'];
                end
                FtrsData_Orig.Properties.VariableNames = FeatureNames_OriginalList';
                
                % append ftr number and _T to ttran features
                FeatureNames_TtranList = FtrsData_Ttran.Properties.VariableNames';
                for n = 1:length(FeatureNames_TtranList)
                    FeatureNames_TtranList{n}=[FeatureNames_TtranList{n,1},'_',num2str(n),'_T'];
                end
                FtrsData_Ttran.Properties.VariableNames = FeatureNames_TtranList';
                
                % concatenate the tables
                FtrsData = [FtrsData_Orig FtrsData_Ttran];
                
                clearvars -except PtData FtrsData LocalDataDir
                
                n_pts=length(unique(PtData.PT));
            end
            
        else % just working with one feature set
           
            LocalData = [LocalDataDir,'\FS_FD_9.mat']; % original timing
            %LocalData = [LocalDataDir,'\FS_FD_10.mat']; % ttran timing
            
            %LocalData = [LocalDataDir,'\FS_FD_wFlow_OriginalT.mat']; % original timing, Flow
            %LocalData = [LocalDataDir,'\FS_FD_wPnasal_OriginalT_noSQRT.mat']; % original timing, Pnasal
            %LocalData = [LocalDataDir,'\FS_FD_wPnasal_OriginalT_wSQRT.mat']; % original timing, Pnasal
            
            % selectively load variables of interest
            load(LocalData,'SleepData','LG_QualityInfo','DataOut',...
                'BreathDataTable','BreathFLDataTable','LocalSignals');
             
            n_pts=length(SleepData); % how many pts are we processing
            
            %savestr = [LocalDataDir, '\FS_FD_9_OriginalT.mat'];
            
            savestr = [LocalDataDir, '\FS_FD_g_wFlow_OriginalT.mat'];
            %savestr = [LocalDataDir, '\FS_FD_g_wPnasal_OriginalT_noSQRT.mat'];
            %savestr = [LocalDataDir, '\FS_FD_g_wPnasal_OriginalT_wSQRT.mat'];
        end
        
        % data set, path to data, and load data
        DataSet = 'FlowDrive'; %'OralAppliance'; %
        switch 1
            case 1
                DriveThreshold_upper = 0.7;
                DriveThreshold_lower = 0.7;
                %savestr_SVM = [LocalDataDir, '\FS_FD_x.mat']; %16.mat']; % save analysis output as...
            case 2
                DriveThreshold_upper = 0.9;
                DriveThreshold_lower = 0.7;
                %savestr_SVM = [LocalDataDir, '\SA_FD_x.mat']; %17.mat']; % save analysis output as...
        end
end

switch DataSet
    case 'FlowDrive'
        DriveAnalysis = 1; % set to one for DriveAnalysis
    case 'OralAppliance'
        DriveAnalysis = 0; % set to zero for non-drive data
    otherwise
        disp('Please select a data set to work with');
        return;
end

%options
ShowFigures = 0;
ShowFlowFigures = 0; %show the flow and edi signals during processing

%%
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
% 'DataOut' {1,p}{1,w}(b,24) % not used anymore
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
for pt=9%1:30%25%[3, 5, 6, 10, 17, 22, 23, 25, 26]% n_pts %22 has high outliers  
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
        % and all the other breath data that is later matched up
        %  (BB_time,BB_Ttot,Ar,win,notAr,Veup,hypnog,pcw,Evt,FL)
        [BB_time,VdriveEdi,VdrivePes,VE,Ar,win,notAr,Veup,hypnog,pcw,BB_Ttot,FL,EType,ApneaB] = ...
            VE_VdriveFLArray(pt,criteria,BreathFLDataTable, BreathDataTable);
        
        if 0 % debugging inverted pnasal and bad breath timing matching
        figure(1); 
        plot(time, Pnasal); hold on;
        plot(BB_time, 0, 'go');
            plot(time, Flow);
        end
        
        %% away from sleep
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
        
        
         %% the hypnog on a per breath basis
         if 0 % not req'd, already output above from VE_VdriveFLArray...
         numWindows = length(BreathDataTable{pt}); HypnogBB=[];
         for w=1:numWindows
             if (size(BreathDataTable{pt}{w},1)==1&&isnan(BreathDataTable{pt}{w}))||...
                     isempty(BreathDataTable{pt}{w})
                 continue
             else
                 HypnogBB = [HypnogBB; [BreathDataTable{pt}{w}.Time_start, BreathDataTable{pt}{w}.hypnog_B]];
             end
         end
         end
         
         %%
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
            
            
            %% then, determine the best VE/Vdrive ratio for Edi and Pes
            
           
            
            
            %% this is the old method for determining the adjusted drive
            %  this is done with moving median value over the wake breaths
            %  that are more than 2 breaths away from sleep
            
            if 1
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
                G_edi_w=sortrows([temp,Gmedian_edi]);
                G_edi_w_indiv=sortrows([temp,data_edi]);
                
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
                G_pes_w=sortrows([temp,Gmedian_pes]);
                G_pes_w_indiv=sortrows([temp,data_pes]);
                
            end
            
            %% this is the new method for determining adjusted drive
            % this is done with the manual reference breath scoring
            if 1                        
                % load reference scoring for this patient
                refscoring = ['..\..\DOS_Scoring\DriveScoring\DriveScoring_Pt_',num2str(pt),'_Ref_TG.mat'];
                load(refscoring);
                
                if exist('BB_DOS','var') % check if the file opened
                    %% confirm times match
                    % BB_Times_All are the breaths in the reference scoring
                    % handle any dropped breaths through indexing by 'ib'
                    [matchedBB, ia_, ib] = intersect(BB_Times_All, BB_time);
                    hypnogBB = hypnog(ib);
                    if size(matchedBB,1) ~= size(BB_time,1)               
                        figure(10); clf(figure(10));
                        subplot(2,1,1);
                        stairs(BB_time); hold on;
                        stairs(BB_Times_All);
                        subplot(2,1,2);
                        stairs(BB_time(ib)); hold on;
                    end  
                     
                    %% find VE/Vdrive for all WAKE ONLY reference breaths
                    refBB = find((BB_DOS(:,2)==4) & (hypnogBB==4));
                    times = BB_time(refBB);
                    [times_,ia,ic] = unique(times);

                    %% edi first
                    data_edi = VE(refBB)./VdriveEdi(refBB);
                    data_edi = data_edi(ia);
                    data_edi(data_edi>prctile(data_edi,99))=NaN;
                    data_edi(data_edi<prctile(data_edi,1))=NaN;
                    times = times_;
                    times(isnan(data_edi))=[]; data_edi(isnan(data_edi))=[];
                    maxt = 1800;
                    Gmedian_edi=NaN*times;
                    for j=1:length(times)
                        temp2 = abs(times(j)-times);
                        weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
                        Gmedian_edi(j)=nansum(data_edi(weights>0).*weights(weights>0));
                    end
                    G_edi_r=sortrows([times,Gmedian_edi]);
                    G_edi_r_indiv=sortrows([times,data_edi]);

                    %% then pes
                    data_pes = VE(refBB)./VdrivePes(refBB);
                    data_pes = data_pes(ia);
                    data_pes(data_pes>prctile(data_pes,99))=NaN;
                    data_pes(data_pes<prctile(data_pes,1))=NaN;
                    times = times_;
                    times(isnan(data_pes))=[]; data_pes(isnan(data_pes))=[];
                    maxt = 1800;
                    Gmedian_pes=NaN*times;
                    for j=1:length(times)
                        temp2 = abs(times(j)-times);
                        weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
                        Gmedian_pes(j)=nansum(data_pes(weights>0).*weights(weights>0));
                    end
                    G_pes_r=sortrows([times,Gmedian_pes]);
                    G_pes_r_indiv=sortrows([times,data_pes]);
            
                    %% find VE/Vdrive for all reference breaths
                    refBB = find(BB_DOS(:,2)==4);
                    times = BB_time(refBB);
                    [times_,ia,ic] = unique(times);

                     %% edi first
                    data_edi = VE(refBB)./VdriveEdi(refBB);
                    data_edi = data_edi(ia);
                    data_edi(data_edi>prctile(data_edi,99))=NaN;
                    data_edi(data_edi<prctile(data_edi,1))=NaN;
                    times = times_;
                    times(isnan(data_edi))=[]; data_edi(isnan(data_edi))=[];
                    maxt = 1800;
                    Gmedian_edi=NaN*times;
                    for j=1:length(times)
                        temp2 = abs(times(j)-times);
                        weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
                        Gmedian_edi(j)=nansum(data_edi(weights>0).*weights(weights>0));
                    end
                    G_edi_r_all=sortrows([times,Gmedian_edi]);
                    G_edi_r_all_indiv=sortrows([times,data_edi]);

                    %% then pes
                    data_pes = VE(refBB)./VdrivePes(refBB);
                    data_pes = data_pes(ia);
                    data_pes(data_pes>prctile(data_pes,99))=NaN;
                    data_pes(data_pes<prctile(data_pes,1))=NaN;
                    times = times_;
                    times(isnan(data_pes))=[]; data_pes(isnan(data_pes))=[];
                    maxt = 1800;
                    Gmedian_pes=NaN*times;
                    for j=1:length(times)
                        temp2 = abs(times(j)-times);
                        weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
                        Gmedian_pes(j)=nansum(data_pes(weights>0).*weights(weights>0));
                    end
                    G_pes_r_all=sortrows([times,Gmedian_pes]);
                    G_pes_r_all_indiv=sortrows([times,data_pes]);
      
                end
            end
            
            %% Pick the G_edi and G_pes to use hereafter
            if 1 % using automated process based on wake 
                G_edi = G_edi_w;
                G_pes = G_pes_w;
            else % using manual reference scoring
                G_edi = G_edi_r;
                G_pes = G_pes_r;
            end
            
            %%
            if 1 %ShowFigures
                
                % show the windows that we analysed
                numWindows = length(BreathDataTable{pt});
                DataWindows = NaN(numWindows,1);
                for w=1:numWindows
                    if (size(BreathDataTable{pt}{w},1)==1&&isnan(BreathDataTable{pt}{w}))||...
                            isempty(BreathDataTable{pt}{w})
                        continue
                    else
                        DataWindows(w) = BreathDataTable{pt}{w}.Time0(1);
                    end
                end
                 

                % do figure
                
                figure(100); clf(figure(100));
                
                ax(1) = subplot(2,1,1); hold on;
                str=['Patient ', num2str(pt)]; title(str);
                plot(G_edi_w(:,1),G_edi_w(:,2),'g-');
                plot(G_edi_r(:,1),G_edi_r(:,2),'r-');
                plot(G_edi_r_all(:,1),G_edi_r_all(:,2),'b-');
                plot(G_edi_w_indiv(:,1),G_edi_w_indiv(:,2),'g.');
                plot(G_edi_r_all_indiv(:,1),G_edi_r_all_indiv(:,2),'b.');
                plot(G_edi_r_indiv(:,1),G_edi_r_indiv(:,2),'r.');
                plot([DataWindows, DataWindows+180], [0 0], 'k', 'LineWidth',10);
                ylabel('Edi');
                
                ax(2) = subplot(2,1,2); hold on;
                plot(G_pes_w(:,1),G_pes_w(:,2),'g-');
                plot(G_pes_r(:,1),G_pes_r(:,2),'r-');
                plot(G_pes_r_all(:,1),G_pes_r_all(:,2),'b-');
                plot(G_pes_w_indiv(:,1),G_pes_w_indiv(:,2),'g.');
                plot(G_pes_r_all_indiv(:,1),G_pes_r_all_indiv(:,2),'b.');
                plot(G_pes_r_indiv(:,1),G_pes_r_indiv(:,2),'r.');
                plot([DataWindows, DataWindows+180], [0 0], 'k', 'LineWidth',10);
                ylabel('Pes');
                legend('Automated', 'Manual (wake)', 'Manual (all)', 'location','southeast');
                linkaxes(ax, 'x');
                
                fig = gcf;
                fig.Color = [1 1 1]; % set background colour to white
                fig.Units = 'inches';
                fig.Position = [20.5 2.5 12 7];
                
                %str = ['..\Figures\', str];
                %savefig(str);
                
            end
            
            if ShowFlowFigures
                figure(101); clf(figure(101));
                ax101(1) = subplot(2,1,1);
                plot(BB_time, VE./VdriveEdi,'kd');hold('on'); % plot everything
                plot(BB_time(a2),VE(a2)./VdriveEdi(a2),'rd'); % at least two breaths away from sleep
                plot(BB_time(a1),VE(a1)./VdriveEdi(a1),'gd'); % more than 4 breaths away from sleep
                plot(BB_time(EType==1), VE(EType==1)./VdriveEdi(EType==1),'b.');
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
            
            %% all breaths are then normalised, i.e. divided by their corresponding VE/Vdrive ratio
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
            
            %% potentially remove major outliers, i.e. only use 5th to 95th centile
            %g_Edi_Adj(g_Edi_Adj>prctile(g_Edi_Adj,95))=NaN;
            %g_Pes_Adj(g_Pes_Adj>prctile(g_Pes_Adj,95))=NaN;
            
            %% dividing by the 'best' ratio, we have normalised to 1.
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
            % the reference data should be exactly one.
            
            if ShowFigures
                figure(100); clf(figure(100));
                ax100(1) = subplot(2,1,1);
                plot(BB_time, g_Edi_Adj,'kd');hold('on'); % plot everything
                plot(BB_time(refBB), g_Edi_Adj(refBB), 'rd'); 
               
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
        %% not doing drive analysis, i.e. just using this for testing and/or training purposes 
        else 
            %manually set these as NaN's
            g_Edi_Adj = NaN(length(BB_time),1);
            g_Pes_Adj = NaN(length(BB_time),1);
            
        end
        
        % Make a table of all breaths from this window
        VarNames = {'PT','Window','Ar','NotAr','A1','A2','Hypnog','Etype','Veup','VE','DriveEdi','DrivePes',...
            'g_Edi_Adj','g_Pes_Adj','ApneaB','BB_time','BB_Ttot'};
        PT = pt*ones(length(VE),1);
        PtDataWin = table(PT, win, Ar, notAr, a1, a2, hypnog, EType, Veup, VE, VdriveEdi, VdrivePes,...
            g_Edi_Adj, g_Pes_Adj, ApneaB, BB_time, BB_Ttot, 'VariableNames', VarNames);
        
        % add the Flow limited table to the current data table
        PtDataWinFL = [PtDataWin,FL];
        
        % add this window data with the FL data, to a growing combined table
        CombinedPtData = cat(1,CombinedPtData,PtDataWinFL);
        
    catch me
    end
end

%% ensure unique breaths in the combined table
% if zeroed, because it is finding the same BB_time in different patients
if 0
    [~,ia,ic] = unique(CombinedPtData.BB_time, 'rows','stable');
    dups=length(CombinedPtData.BB_time)-length(ia);
    str = ['Removing ', num2str(dups), ' duplicates and sorting data']; disp(str);
    PtData = CombinedPtData(ia,:);
else
    % just copy the table
    PtData = CombinedPtData;
end

%% tidy up the workspace, i.e. only need PtData from here onwards
clearvars -except PtData n_pts DataSet savestr_SVM savestr ...
    DriveThreshold_lower DriveThreshold_upper DriveAnalysis PtData_Orig
%clear LocalSignals

%% save PtData
if 1
   save(savestr, 'PtData', '-v7'); 
end

%% remove labelled apnea breaths (as labelled by VEfromFlow)
ApneaFlow = (PtData.ApneaB==1); %nnz(PtData.ApneaB)
str = ['Removing ', num2str(nnz(ApneaFlow)), ' apneic breaths']; disp(str);
PtData(ApneaFlow,:)=[];

%% remove breaths with < threshold% flow below eupnoea
% where VE < 10% of Veup
LowFlow = PtData.VE < (PtData.Veup*0.1); %nnz(PtData.VE < (PtData.Veup*0.1))
str = ['Removing an additional ', num2str(nnz(LowFlow)), ' low flow breaths (<10% of eupnea)']; disp(str);
PtData(LowFlow,:)=[];

%% Find extreme outliers: Part 1
% for e.g. VTi_VTe should be around 1, so values of 100 are clearly abnormal.
firstdatacol = 18;  % confirm this during procesing
if 0
    %original simplistic method, using basic prctile
    centileupper = prctile(PtData{:,firstdatacol:end},99.9);
    centilelower = prctile(PtData{:,firstdatacol:end},0.1);
    
    %outlier detection using quantiles, where
    % F1=Q1-1.5*IQR, and F2=Q3+1.5*IQR, with F1 > outliers > F2
    IQR_ = iqr(PtData{:,firstdatacol:end});
    Q1_ = quantile(PtData{:,firstdatacol:end}, 0.25);
    Q3_ = quantile(PtData{:,firstdatacol:end}, 0.75);
    prctilefactor=1.5; %e.g. set to about 5 to be equivalent to Scotty's method
    centileupper = Q3_ + prctilefactor*IQR_;
    centilelower = Q1_ - prctilefactor*IQR_;
else %SS comment: SS revised version labels outliers e.g. if 2 fold above the 95-50 difference, e.g. median = 5, 95th centile = 10; outliers >15.
    prctilethres=99; %e.g. 95
    prctilefactor=3; %e.g. 2
    centileupper = prctile(PtData{:,firstdatacol:end},50) + (prctile(PtData{:,firstdatacol:end},prctilethres)-prctile(PtData{:,firstdatacol:end},50))*prctilefactor;
    centilelower = prctile(PtData{:,firstdatacol:end},50) - (prctile(PtData{:,firstdatacol:end},50)-prctile(PtData{:,firstdatacol:end},100-prctilethres))*prctilefactor;
end
outliers = (PtData{:,firstdatacol:end} > centileupper | ...
    PtData{:,firstdatacol:end} < centilelower);
%upper_outliers = nnz(PtData{:,firstdatacol:end} > centileupper);
%lower_outliers = nnz(PtData{:,firstdatacol:end} < centilelower);
filler = false(size(PtData,1),firstdatacol-1); % exclude the first n columns of PtData, which are not feature data
outliers = ([filler outliers]);
perBB = any(outliers,2);
str = [ num2str(nnz(outliers)),' feature outliers were found in ', num2str(nnz(perBB)), ' breaths']; disp(str);
%SS comment: I think better to simply set the values to be equal to the upper/lower boundary

%% Remove extreme outliers: Part 2
if 0
    % this is a horribly inneficient method of covering the entire mat.
    % logical indexing would be much better, but I can't see an easy way
    % to do this in the table structure
    for m=1:size(PtData,1)
        for n=1:size(PtData,2)
            if outliers(m,n)
                PtData{m,n}=NaN;
            end
        end
    end
else %SS version, replaces outliers with limit value
    for n=1:length(centileupper)
        PtData{:,n+firstdatacol-1}(PtData{:,n+firstdatacol-1} > centileupper(n)) = centileupper(n);
        PtData{:,n+firstdatacol-1}(PtData{:,n+firstdatacol-1} < centilelower(n)) = centilelower(n);
    end
end

%% set up categorical thresholds for Mild, Moderate and Severe FL
% and also count number of hypops and arousals
%
% three category model
%   NFL >90%
%   intermediate, 30<x<90%
%   FL <30%
%
% five category model
%   Normal >90
%   Mild FL 70-90
%   Mod FL 50-70
%   Sev FL 30-50
%   V.Sev <30

% inclusive categorical thresholds, i.e.
Mild_FL = PtData.VE < (PtData.Veup*0.7);
Mod_FL = PtData.VE < (PtData.Veup*0.5);
Sev_FL = PtData.VE < (PtData.Veup*0.3);

if DriveAnalysis
    NonFL_Drive = PtData.g_Edi_Adj > (0.9);
    Mild_FL_Drive = PtData.g_Edi_Adj < (0.9);
    Mod_FL_Drive = PtData.g_Edi_Adj < (0.7);
    Sev_FL_Drive = PtData.g_Edi_Adj < (0.5);
    VSev_FL_Drive = PtData.g_Edi_Adj < (0.3);
end

% exclusive categorical thresholds, i.e.
% 50% < Mild_FL_e < 70% and  30% < Mod_FL_e < 50%  and  10% < Sev_FL_e < 30%
Mild_FL_e = (PtData.VE < (PtData.Veup*0.7) & PtData.VE > (PtData.Veup*0.5));
Mod_FL_e = (PtData.VE < (PtData.Veup*0.5) & PtData.VE > (PtData.Veup*0.3));
Sev_FL_e = (PtData.VE < (PtData.Veup*0.3) & PtData.VE > (PtData.Veup*0.1));
ALL_FL_e = zeros(length(Mild_FL_e),1); ALL_FL_e(Mild_FL_e)=1; ALL_FL_e(Mod_FL_e)=2; ALL_FL_e(Sev_FL_e)=3;

if DriveAnalysis
    NonFL_e_Drive = PtData.g_Edi_Adj > (0.9);
    Mild_FL_e_Drive = (PtData.g_Edi_Adj < 0.7) & (PtData.g_Edi_Adj > 0.5);
    Mod_FL_e_Drive = (PtData.g_Edi_Adj < 0.5) & (PtData.g_Edi_Adj > 0.3);
    Sev_FL_e_Drive = (PtData.g_Edi_Adj < 0.3) & (PtData.g_Edi_Adj > 0.1);
    ALL_FL_e_Drive = zeros(length(Mild_FL_e_Drive),1); ALL_FL_e_Drive(Mild_FL_e_Drive)=1;
    ALL_FL_e_Drive(Mod_FL_e_Drive)=2; ALL_FL_e_Drive(Sev_FL_e_Drive)=3;
end

% 5=Central 6=Obstructive 7=OHypopnea 15=Mixed 16=CentralHypopnea
Hypop = (PtData.Etype==6 | PtData.Etype==7);
Arous = (PtData.Etype==0 & PtData.Ar==1);

FL_byVE = table([nnz(Mild_FL);nnz(Mild_FL_e)],[nnz(Mod_FL);nnz(Mod_FL_e)],[nnz(Sev_FL);nnz(Sev_FL_e)],...
    'VariableNames',{'Mild','Mod','Sev'},'RowNames',{'Inclusive','Exclusive'})

if DriveAnalysis
    FL_byDrive = table([nnz(Mild_FL_Drive);nnz(Mild_FL_e_Drive)],[nnz(Mod_FL_Drive);nnz(Mod_FL_e_Drive)],[nnz(Sev_FL_Drive);nnz(Sev_FL_e_Drive)],...
        'VariableNames',{'Mild','Mod','Sev'},'RowNames',{'Inclusive','Exclusive'})
end

table(nnz(Mild_FL), nnz(Hypop), nnz(Arous),'VariableNames',{'Mild_FL','Hypop','Arousal'})

table(nnz(Mild_FL | Hypop), nnz(Mild_FL & Hypop), 'VariableNames',{'Mild_FL_OR_Hypop','Mild_FL_AND_Hypop'})

%% how many breaths of each type/category do we have per patient
n_pts=54;
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

%table(min(BBperPt(:,2)), max(BBperPt(:,2)), 'VariableNames', {'Min','Max'})

% figures for number of breaths per patient, and breath categories
if 0
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
end

%% BB per pt, based on VE:VDrive
if DriveAnalysis
    BBperPt_Drive = NaN(n_pts,12);
    for pt = 1:n_pts
        BBperPt_Drive(pt,1) = pt;
        BBperPt_Drive(pt,2) = size(PtData(PtData.PT==pt,1),1); %total
        BBperPt_Drive(pt,3) = size(PtData(PtData.PT==pt & PtData.Ar==1,1),1); % arousal
        BBperPt_Drive(pt,4) = size(PtData(PtData.PT==pt & PtData.A1==1,1),1); % A1
        BBperPt_Drive(pt,5) = size(PtData(PtData.PT==pt & PtData.A2==1,1),1); % A2
        BBperPt_Drive(pt,6) = size(PtData(PtData.PT==pt & PtData.Etype~=0,1),1); % event
        BBperPt_Drive(pt,7) = size(PtData(PtData.PT==pt & Sev_FL_e_Drive==1,1),1); % severe only breaths
        BBperPt_Drive(pt,8) = size(PtData(PtData.PT==pt & Mod_FL_e_Drive==1,1),1); % moderate only breaths
        BBperPt_Drive(pt,9) = size(PtData(PtData.PT==pt & Mild_FL_e_Drive==1,1),1); % mild only breaths
        BBperPt_Drive(pt,10) = size(PtData(PtData.PT==pt & Hypop,1),1); % hypops
        BBperPt_Drive(pt,11) = size(PtData(PtData.PT==pt & Arous,1),1); % arousals
        BBperPt_Drive(pt,12) = size(PtData(PtData.PT==pt & (Hypop | Mild_FL_Drive),1),1); % hypop or MildFL
    end
    varnames = {'PT','Total','Ar','A1','A2','Evts','SevereFL','ModFL','MildFL','Hypop','Arousal','HypMildFL'};
    BBperPt_Drive = table(BBperPt_Drive(:,1),BBperPt_Drive(:,2),BBperPt_Drive(:,3),BBperPt_Drive(:,4),...
        BBperPt_Drive(:,5),BBperPt_Drive(:,6),BBperPt_Drive(:,7),BBperPt_Drive(:,8),BBperPt_Drive(:,9),...
        BBperPt_Drive(:,10),BBperPt_Drive(:,11),BBperPt_Drive(:,12),...
        'variablenames',varnames)
    clear varnames
end

%% before we go any further, make a backup of the data
PtDataBackup=PtData;

% save(savestr_SVM, '-v7')

%% re-instate the original data from backup
PtData=PtDataBackup;

%% have a look at some of the data
if 0
    % check agreement between drive measures,
    % i.e. how do the two gold standards compare
    x = PtData.Edi_Adj; x(x<0.01)=NaN;
    y = PtData.Pes_Adj; y(y<0.01)=NaN;
    figure(103); clf(figure(103));
    scatterhist(x, y, 'Group', PtData.Ar==1, 'Marker','..','MarkerSize',[3,3]); hold on;
    xlabel('Edi'); ylabel('Pes');
    axis square; r(1) = refline(1,0); r(1).Color='y'; r(1).LineWidth = 1.5;
    lsline; legend off; legend('..','Ar');
    
    % how does drive look when compared with ventilation
    % grouped by FL magnitude
    x = PtData.Edi_Adj; x(x<0.01)=NaN;
    y = PtData.VE; y(y==0)=NaN;
    figure(106); clf(figure(106));
    %scatterhist(x, y, 'Group', (PtData.Evt~=1 & PtData.Ar==0), 'Marker','..','MarkerSize',[3,3]); hold on;
    scatterhist(x, y, 'Group', ALL_FL_e, 'Marker','..','MarkerSize',[3,3]); hold on;
    xlabel('EdiAdj'); ylabel('VE');
    axis square; r(1) = refline(1,0); r(1).Color='y'; r(1).LineWidth = 1.5;
    lsline;legend off;
    
    % how do the kurtosis and skewness methods compare
    x = PtData.SkewDistInsp;
    y = PtData.SkewDataInsp;
    figure(104); clf(figure(104));
    scatter(x, y, 'k.'); hold on;
    xlabel('Dist'); ylabel('Data');
    scatter(x(PtData.Ar==1), y(PtData.Ar==1), 'r.');
    axis square; r(1) = refline(1,0); r(1).Color='y'; r(1).LineWidth = 1.5;
    title('Skew');
    x = PtData.KurtDistInsp;
    y = PtData.KurtDataInsp;
    figure(105); clf(figure(105));
    scatter(x, y, 'k.'); hold on;
    xlabel('Dist'); ylabel('Data');
    scatter(x(PtData.Ar==1), y(PtData.Ar==1), 'r.');
    axis square; r(1) = refline(1,0); r(1).Color='y'; r(1).LineWidth = 1.5;
    title('Kurt');
    
    % what about other measures - see ScatMe function
    for m=18 % for random columns
        x = PtData.Edi_Adj; x(x<0.01)=NaN;
        y = PtData{:,m}; %y(y==0)=NaN;
        figure(100+m); clf(figure(100+m));
        scatter(x(PtData.NotAr==1), y(PtData.NotAr==1), 'k.'); hold on;
        scatter(x(PtData.Ar==1), y(PtData.Ar==1), 10, 'r.'); hold on;
        xlabel('Edi'); ylabel('...');
        axis square; r(1) = refline(1,0); r(1).Color='y'; r(1).LineWidth = 1.5;
    end
    
end

%% setup the binary response variable for classification
% set up grouping and criteria data.
% could remove pts who are not contributing. See BBperPt_table for guidance
switch DataSet
    case 'FlowDrive'
        Gtest_continuous = PtData.g_Edi_Adj;
        
        Gtest = ones(size(PtData,1),1)*2; % set a vector of two
        Gtest(Gtest_continuous>=DriveThreshold_upper) = 1; % NFL
        Gtest(Gtest_continuous<DriveThreshold_lower) = 0; % FL
        
        exclude = Gtest==2;
        Gtest_continuous(exclude)=[];
        Gtest(exclude)=[];
        PtData(exclude,:)=[];
        if exist('FtrsData','var')
            FtrsData(exclude,:)=[];
        end
        nnz(exclude)
        str = [num2str(nnz(~exclude)), ' breaths were labelled as either FL or NonFL']; disp(str)
        str = ['Breakdown: ', num2str(nnz(Gtest==1)), ' FL, and ', num2str(nnz(Gtest==0)), ' NonFL']; disp(str)
        str = [num2str(nnz(exclude)), ' breaths were unlabelled and removed']; disp(str)
        
        if 0 % figure, comparison of continuous measures
            a = 0.7; % value for line color
            VEVEup = PtData.VE./PtData.Veup;
            figure(99); clf(figure(99));
            scatter(Gtest_continuous, VEVEup, 1, [0.5 0.05 0.05]);
            hold on; lsline;
            currentAxes = gca;
            ylimits = currentAxes.YLim; % call this explicitly
            xlimits = currentAxes.XLim; % call this explicitly
            line([0.9 0.9],ylimits,'LineWidth',1,'Color',[a a a],'Parent',currentAxes);
            line([0.7 0.7],ylimits,'LineWidth',1,'Color',[a a a],'Parent',currentAxes);
            line([0.5 0.5],ylimits,'LineWidth',1,'Color',[a a a],'Parent',currentAxes);
            line(xlimits,[0.9 0.9],'LineWidth',1,'Color',[a a a],'Parent',currentAxes);
            line(xlimits,[0.7 0.7],'LineWidth',1,'Color',[a a a],'Parent',currentAxes);
            line(xlimits,[0.5 0.5],'LineWidth',1,'Color',[a a a],'Parent',currentAxes);
            xlabel('VEVDrive'); ylabel('VEVEup');
            r = corr(Gtest_continuous, VEVEup, 'type', 'Pearson', 'rows', 'pairwise');
            str=['Comparison of continuous measures. R^2 = ', num2str(r^2,1)];
            title(str);
        end
        if 0 % numerical comparison of classifications, i.e. confusion mat
            Th_NFL = 0.9; % greater
            %Intermediate, between these margins, equiv to Mild and Mod
            Th_FL = 0.7; % below, equiv to Sev and VSev
            BB_NFL = Gtest_continuous>=Th_NFL;
            BB_FL = Gtest_continuous<=Th_FL;
            BB_Inter = ~BB_NFL & ~BB_FL;
            
            G_FD(BB_NFL) = 0;
            G_FD(BB_Inter) = 1;
            G_FD(BB_FL) = 2;
            
            G_OA = ones(size(PtData,1),1); % set a vector of ones
            G_OA(PtData.A1==1 & PtData.VE>(PtData.Veup*0.9) & PtData.Etype==0)=0; % set NFL as 0
            G_OA(PtData.NotAr==1 & PtData.VE<(PtData.Veup*0.7) & (PtData.Etype==6 | PtData.Etype==7))=2; % set FL as 2
            
            %[c, order] = confusionmat(G_FD,G_OA); c
            [tbl,chi2,p,lbl] = crosstab(G_FD,G_OA); tbl
            
            %nnz(G_FD==0)
        end
    case 'OralAppliance'
        Gtest = ones(size(PtData,1),1)*2; % set a vector of two's,
        % that we will mark NonFL as zero, FL as one, and remove the twos
        
        % NonFL breaths are:
        %  Ar, by at least 2 breaths (so using A1, which is threshold)
        %  VE > Veup*0.9, VE is good
        %  Etype == 0, not a clinically scored event
        Gtest(PtData.A1==1 & PtData.VE>(PtData.Veup*0.9) & PtData.Etype==0)=0;
        
        % FL breaths are:
        %  NotAr
        %  VE < Veup*0.7,
        %  Etype == 6 or 7
        % use the clinically scored events (hypop/obstructive) as true labels
        Gtest(PtData.NotAr==1 & PtData.VE<(PtData.Veup*0.7) & (PtData.Etype==6 | PtData.Etype==7))=1;
        
        %% What are the Gtest2 breaths
        if 0
            VEVEup = PtData.VE./PtData.Veup;
            Ar1 = PtData.A1==0;
            nnz(VEVeup(Gtest==2)>0.9&Ar1(Gtest==2))
            nnz(Gtest==2)
        end
        
        %% remove breaths that are neither FL, nor NonFL
        exclude = Gtest==2;
        Gtest(exclude)=[];
        PtData(exclude,:)=[];
        nnz(exclude)
        str = [num2str(nnz(~exclude)), ' breaths were labelled as either FL or NonFL']; disp(str)
        str = ['Breakdown: ', num2str(nnz(Gtest==1)), ' FL, and ', num2str(nnz(Gtest==0)), ' NonFL']; disp(str)
        str = [num2str(nnz(exclude)), ' breaths were unlabelled and removed']; disp(str)
        
    otherwise
        pt_exclude = [];
        PtData(pt_exclude,:)=[]; %truncate the PtData
        n_pts = size(BBperPt,1);
end

%%  Set up criteriaR and Yvariable
criteriaR=1*Gtest;
criteriaR(isnan(Gtest))=NaN;
Yvariable=single(Gtest);

%% Weight the input breaths
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
    I1=find(pt_in==table2array(PtData(:,1))&criteriaR==1);
    I2=find(pt_in==table2array(PtData(:,1))&criteriaR==0);
    if length(I1)<=min_breaths || length(I2)<=min_breaths
        weights(I1)=0; % this pt doesn't have enough breaths
        weights(I2)=0; % so set the weights to zero
        pt_out = cat(1,pt_out,pt_in); % keep a list of pts to be excluded
    else
        weights(I1)=(1/length(I1))/2;
        weights(I2)=(1/length(I2))/2;
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

str = ['New Breakdown: ', num2str(nnz(Yvariable)), ' FL, and ', num2str(size(Yvariable,1)-nnz(Yvariable)), ' NonFL']; disp(str)

%% report on how many NonFL and FL breaths from each pt
% by Yvariable (or criteriaR) zero or one definition
SumData = NaN(length(pts), 4);
for n=1:length(pts)
    SumData(n,1) = pts(n);
    SumData(n,2) = nnz(table2array(PtData(:,1))==pts(n) & criteriaR==0);
    SumData(n,3) = nnz(table2array(PtData(:,1))==pts(n) & criteriaR==1);
    SumData(n,4) = nnz(table2array(PtData(:,1))==pts(n));
end
varnames = {'PT','NonFL','FL','Total'};
FL_BBperPt = table(SumData(:,1),SumData(:,2),SumData(:,3),SumData(:,4),...
    'variablenames',varnames)
clear varnames

%% >>>>>>>>>>>>>>>>>>>   STOP   <<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    % unless you want to run a LDA, or univariate analysis,
    % you can skip this section, and go to the SVM section below
    
    %% test performance of a few manually selected feature sets
    % note these numbers no longer align, should use variable names
    VEVeup = table(PtData{:,10}./PtData{:,9},'VariableNames',{'VEVeup'});
    clear Tbl
    TblType = 3;
    switch TblType
        case 1
            % Ti_Ttot, MIF_PIF, AA_InspFlat8020, SS_Area, InvParabI, SkewDistInsp, KurtDistInsp, SeriesIEtime, Ali_MostPromPeak
            Tbl = PtData(:,[18, 22, 30, 43, 46, 49, 50, 66, 96]);
        case 2
            % Ti_Ttot, Teschler, SinI, SeriesIEflow
            Tbl = PtData(:,[18,42,44,65]);
        case 3
            % Ti_Ttot, AA_InspFlat8020, Ali_MostPromPeak
            Tbl = PtData(:,[18,30,96]); %
        case 4
            % this is equivalent to Run_1 in svm
            Tbl = PtData(:,[18, 22, 30, 43, 46, 49, 50, 66, 96]);
        case 5
            % this is equivalent to Run_2 in svm
            Tbl = PtData(:,[96, 40, 107, 22, 37, 75, 43, 48, 75]);
        otherwise
            % Just VE ./ Veup
            Tbl = [VEVeup]; % this should result in near perfect prediction
    end
    
    %% LDA run
    Mdl = fitcdiscr(Tbl,Yvariable);%,'Weights',weights);
    LDA_MDL = predict(Mdl, Tbl);
    
    %% LDA results
    x = Yvariable';
    y = LDA_MDL';
    
    table([nnz(x);length(x)], [nnz(y);length(y)], 'VariableNames',{'Label','Predict'}, 'RowNames', {'Pos','Total'} )
    
    % ROC AUC SEM Sens Spec
    ShowFigures = 0;
    [thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(x,y,ShowFigures);
    
    confusionmat(x,y)
    figure(501); clf(figure(501)); plotconfusion(x,y)
    %figure(); plotroc(x,y)
    
    
    %% Boxplots
    % warning, potentially lots and lots of boxplots
    k=18;
    group = [ones(nnz(Yvariable==0),1)*1;ones(nnz(Yvariable==1),1)*2];
    for outer = 1:11
        figure(504+outer); clf(figure(504+outer));
        for inner = 1:9
            if k==116
                break
            end
            subplot(3,3,inner);
            data = [table2array(PtData(Yvariable==0,k)); table2array(PtData(Yvariable==1,k))];
            boxplot(data,group); hold on;
            refline(0,nanmedian(table2array(PtData(Yvariable==0,k))));
            refline(0,nanmedian(table2array(PtData(Yvariable==1,k))));
            str = [num2str(k)];
            xlabel(str);
            k = k+1;
        end
    end
    
    % handy figure closing
    for outer = 1:11
        close(figure(504+outer));
    end
    
    %% >>>>>>>>>>>>>>>>>>> Univariate analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<
    ShowFigures = 0;
    ftrs=18:1:115;
    AUCData=NaN(size(ftrs,2),7);
    for fts=[ftrs]
        x = PtData{:,fts};
        if any(isinf(x))%isnan(x)|
            continue % this includes inf etc
        end
        Mdl_ = fitcdiscr(x,Yvariable,'Weights',weights);
        Mdl_Pred = predict(Mdl_, x);
        if nnz(Mdl_Pred)==0 || nnz(~Mdl_Pred)==0
            continue
        else
            [thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(Mdl_Pred,Yvariable,ShowFigures); %
            AUCData(fts-17,1)=thresX;
            AUCData(fts-17,2)=AUC;
            AUCData(fts-17,3)=SEM;
            AUCData(fts-17,4)=p;
            AUCData(fts-17,5)=posclass;
            AUCData(fts-17,6)=sensspec(1);
            AUCData(fts-17,7)=sensspec(2);
        end
    end
    
    AUCData(:,1)=ftrs;
    AUCDataSorted=sortrows(AUCData, 2, 'descend');
    exclude_=isnan(AUCDataSorted(:,2))|AUCDataSorted(:,6)==0|AUCDataSorted(:,7)==0;
    AUCDataSortedTrimmed=AUCDataSorted;
    AUCDataSortedTrimmed(exclude_,:)=[];
    AUCDataSortedTrimmed(:,[4,5])=[];
    
end

%% set up the features list
if 1 % set to one is using a single feature set, zero for pooled features
    ftrset = 4;
    switch ftrset
        case 0
            % for the initial testing, we used a skip subset of all data
            rangei=18:5:115;
        case 1 % run1, then we manually selected some features
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
            % was 16:1;115 in MIFLv2
            % was 17:1:116 in MIFLv3
            % rangei = [18:1:118]; % in MIFLv4
            rangei = [firstdatacol:1:size(PtData,2)]; % more robust method
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
        case 9
            % testing multinomial regression, using previous top three ftrs
            rangei = 15+[80, 66, 8];    % ttran timing
            rangei = 15+[24, 66, 7];    % original timing
            
        otherwise
            rangei = [];
            disp('Please set features list');
    end
    Amatrix = PtData{:,rangei};
else % using predictors from two feature sets, as such, 
    % Amatrix is determined from the features data, and is set as:
    Amatrix = FtrsData{:,:};
end


%% remove ftrs with more than SetThreshold of nonfinite data
SetThreshold = 0.01; % 0.01, have also run at 0.001
a=[];
for i = 1:size(Amatrix,2)
    a = [a;nnz(~isfinite(Amatrix(:,i)))];
end
Threshold = round(size(Amatrix,1)*SetThreshold);
ftrsToExclude = find(a>Threshold);

% FD, original timing, 0.001, 21 ftrs
% FD, original timing, 0.01, 20 ftrs
% FD, original timing, 0.05, 17 ftrs

Amatrix(:,ftrsToExclude)=[];

BBwithNans = find(any(isnan(Amatrix),2));
Amatrix(BBwithNans,:)=[];
weights(BBwithNans)=[];
criteriaR(BBwithNans)=[];
Gtest(BBwithNans)=[];
Gtest_continuous(BBwithNans)=[];
Yvariable(BBwithNans)=[];
PtData(BBwithNans,:)=[];
%rangei = [18:1:size(PtDataBackup,2)]; 
rangei(ftrsToExclude)=[];
%Tbl_rangei = PtDataBackup(:,rangei); % make a human readable data table equiv to Amatrix
Tbl_rangei = PtData(:,rangei); % make a human readable data table equiv to Amatrix

%% plot NaN's in data
if 0
    DataMap = zeros(size(Amatrix));
    for i = 1:size(Amatrix,2)
        DataMap(:,i) = (~isfinite(Amatrix(:,i))*90);
    end
    figure(1)
    image(DataMap(:,:))
    ylabel('BB''s');
    xlabel('Features');
    title('Nonfinite data (All Ftrs)');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    fig.Position = [1 1 12 6];
    ax = gca;
    ax.TickDir = 'out';
    box('off');
end

%% Center and scale data, normalise standardize etc.
% (1) look at each feature, does it make sense to standardize
% (2) how will i interpret predictions based on a model built using standardized input
if 0
    % zscore (built-in)
    [Amatrix,MatMu,MatStd] = zscore(Amatrix);
    if 0
        x = Amatrix;
        % zscore as an anonymous function
        Zsc = @(x) (x - mean(x))./std(x);
        Zx = Zsc(x);
        
        % zscore as an anonymous function, with 2 x std
        Zsc2 = @(x) (x - mean(x))./(2*std(x));
        Zx2 = Zsc2(x);
    end
end

%% >>>>>>>>>>>>>>>>>>>   LOAD   <<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    close all
    clear
    clc
    
    LocalDataDir = 'C:\PSG_Data\QAO';
    savestr_SVM = [LocalDataDir, '\SA_FD_12.mat'];
    load(savestr_SVM);
end

%% >>>>>>>>>>>>>>>>>>>   SAVE   <<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    save(savestr, '-v7')
end

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


%% SVM classification settings
coststr = '1-SensplusSpec_';
sigma = 2; %1 (default), Scotty set at 2, app training best was 3.2
svpmethodstr = 'SMO'; %'ISDA' | 'L1QP' | 'SMO'
kernelfunction = 'rbf'; %  'linear'; %'polynomial'; %
maxNfeaturesJ=10;
tempignorevars=[];
forcedvars=[];
if 1
    downsamplefactor=5;
    I=ceil(rand(1,1)*downsamplefactor):downsamplefactor:length(criteriaR);
    % if starting with 80k observations, 100 c 20
    % ds   100 c 20        53 c 15        53 c 10
    % 10   90 s / step
    % 7                    100 s / step
    % 5    300 s / step                    200 s / step
    % 4    540 s / step    340 s / step
    % 3
    % none
else
    I = [1:1:length(criteriaR)]; % set to process all observations
end

%% SVM forward select run
t_start_outer = clock;
D = cell(length(pts),1);
Ilist_= NaN(length(pts),maxNfeaturesJ); % the list of features for each pt
%PredTAll_=NaN(N,1); % this is the predicted output, made by concatenating the output from each pt
PerPTPerformance_Train{length(pts)}=[]; % this is the performance in training data
PerPTPerformance_Test{length(pts)}=[]; % this is the performance in the left out pt
PerPTModel_Train{length(pts)}=[]; % this is the training model output for each pt
PerPTModel_Test{length(pts)}=[]; % this is the test model output for each pt
TrainingBreaths_ = NaN(length(pts),4); % the number of breaths used for each patient
TestingBreaths_ = NaN(length(pts),4); % the number of breaths used for each patient
count=0;
for pt_in=pts
    count=count+1; disp(' ');
    str=(['Processing patient ID ', num2str(pt_in), '; number ', num2str(count),' of ', num2str(length(pts))]); disp(str);
    progressbar('Patient','Loop');  % show progressbar
    progressbar(count/length(pts),[]);
    % set up the cross validation range
    rangevalidate=pt_in==PtData{I,1};
    
    % run the forward selection process, and return a list of features and
    % performance as each feature is added
    t_start_inner = clock;
    if rerun % this option includes prior learning
        [Ilist,PtPerf1,PtPerf2,SVMModel1,SVMModel2, TrainingBreaths, TestingBreaths]=svmforwardselect_perpatient(Amatrix(I,:),criteriaR(I),...
            tempignorevars,Yvariable(I),maxNfeaturesJ,rangevalidate,forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights(I),...
            prior_Ilist(count,:), prior_perf1{count}, prior_perf2{count});
    else % no prior learning
        [Ilist,PtPerf1,PtPerf2,SVMModel1,SVMModel2, TrainingBreaths, TestingBreaths]=svmforwardselect_perpatient(Amatrix(I,:),criteriaR(I),...
            tempignorevars,Yvariable(I),maxNfeaturesJ,rangevalidate,forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights(I),...
            [],[],[]);
    end
    t_end_inner = clock;
    delta_t = etime(t_end_inner, t_start_inner); % delta in seconds
    D{count,1} = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
    str = ['Patient processing time: ', char(D{count,1}), ' (hh:mm:ss)']; disp(str);
    
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
    TrainingBreaths_(count,:) = TrainingBreaths;
    TestingBreaths_(count,:)  = TestingBreaths;
    
end

t_end_outer = clock;
delta_t = etime(t_end_outer, t_start_outer); % delta in seconds
D_ = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
str = ['Total processing time: ', char(D_), ' (hh:mm:ss)']; disp(str);

progressbar(1,1); % force close progressbar

save(savestr_SVM, '-v7')

%% =======================================================================
%               End of SVM stepwise forward select process
% ========================================================================

%% Save
if 0
    LocalDataDir = 'C:\PSG_Data\QAO';
    savestr_SVM = [LocalDataDir, '\SA_FD_15'];
    save(savestr_SVM, '-v7'); %
end

%% Load
if 0
    close all
    clear
    clc
    LocalDataDir = 'C:\PSG_Data\QAO';
    savestr_SVM = [LocalDataDir, '\SA_FD_15'];
    load(savestr_SVM);
end

%% how many breaths were used during training and testing phases
BB_duringTraining = nan(length(pts),maxNfeaturesJ);
for pt=1:length(PerPTModel_Train)
    if isempty(PerPTModel_Train{1,pt})
        continue
    end
    for ftr=1:size(PerPTModel_Train{1,pt},2)
        BB_duringTraining(pt,ftr) = PerPTModel_Train{1,pt}{1,ftr}.NumObservations;
    end
end

BB_duringTesting = nan(length(pts),maxNfeaturesJ);
for pt=1:length(PerPTModel_Test)
    if isempty(PerPTModel_Test{1,pt})
        continue
    end
    for ftr=1:size(PerPTModel_Test{1,pt},2)
        BB_duringTesting(pt,ftr) = PerPTModel_Test{1,pt}{1,ftr}.NumObservations;
    end
end

[BB_duringTraining(:,1) ,BB_duringTesting(:,1)]

%% what features appear most frequently
% note: this doesn't sort by order of appearance in the list
numpts = length(pts);
numftrs = size(Amatrix,2);
centers = 1:1:numftrs;
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
ax.YTick=[0:1:numpts];
ylabel('Frequency');
xlabel('Feature #');
title('Histogram of features (79 c 20)');
box off
%savefig(fig, 'c:\Temp\HistogramFtrs_98c20.fig');

%% Try the builtin "Classification Learner App"
BestFtrs = Tbl_rangei{:,[71 61 25 10 1 40 62 68 78 79 2]}; % as matrix
DataIn = [Yvariable BestFtrs];

BestFtrs = Tbl_rangei(:,[71 61 25 10 1 40 62 68 78 79 2]); % as table
DataIn = [RespVar BestFtrs];

RequiredFeatures = table([1:size(DataIn,2)]', DataIn.Properties.VariableNames', ...
    'VariableNames', {'Ftr','Name'});

classificationLearner % works with 2016b onwards, not sure what else.

regressionLearner % needs Matlab 2017a onwards

%% Top 20 - run in another feature set
% get the top 20 ftrs, and use only these
[i,a]=sort(counts,'descend'); [i(1:20)',a(1:20)']

clearvars -except a
LocalDataDir = 'C:\PSG_Data\QAO';
savestr_SVMPrep = [LocalDataDir, '\SA_OA_4_Prep.mat'];
load(savestr_SVMPrep);
savestr_SVM = [LocalDataDir, '\SA_OA_9_Run.mat']; % save the output as this

coststr = '1-SensplusSpec_'; % 'Sens_' 'Spec_' 'Acc_'
sigma = 2; %1 (default) | 'auto' | positive scalar, Scotty set at 2
svpmethodstr = 'SMO'; %'ISDA' | 'L1QP' | 'SMO'
kernelfunction = 'rbf'; %'linear' (default) | 'gaussian' | 'rbf' | 'polynomial' | function name
maxNfeaturesJ=20;
rerun = 0;

N=size(PtData,1); % do again, as weighting above may have changed the number
rangei = [16:1:114];
Tbl_rangei = PtData(:,rangei);
Amatrix = Tbl_rangei{:,a(1:20)};    % mat format
Tbl_rangei =  Tbl_rangei(:,a(1:20));% table format
% set up the ignored and forced variables
tempignorevars=[];
forcedvars=[];

%% Table of features selected
if 0
    for j=1:6%size(Ilist_,2)
        disp(' ');
        str_title=['Fwd Select #',num2str(j)]; disp(str_title);
        ftrs_found = unique(Ilist_(:,j));
        for k=1:size(ftrs_found,1)
            str_data=['Ftr:', num2str(ftrs_found(k)), ' QTY: ',num2str(nnz(Ilist_(:,j)==ftrs_found(k)))];
            disp(str_data);
        end
    end
end

if 1
    FeatureNames = table([1:size(Tbl_rangei,2)]', Tbl_rangei.Properties.VariableNames', ...
        'VariableNames', {'Ftr','Name'});
else
    FeatureNames = table([1:size(Ftrs2020,2)]', Ftrs2020.Properties.VariableNames', ...
        'VariableNames', {'Ftr','Name'});
end

FeatureNames = table([1:size(FtrsData,2)]', FtrsData.Properties.VariableNames', ...
    'VariableNames', {'Ftr','Name'});

for j=1:10%size(Ilist_,2)
    disp(' ');
    str_title=['Fwd Select #',num2str(j)]; disp(str_title);
    ftrs_found = unique(Ilist_(:,j));
    for k=1:size(ftrs_found,1)
        str_data=['x' , num2str(nnz(Ilist_(:,j)==ftrs_found(k))), ...
            ', #', num2str(ftrs_found(k)), ...
            ', ', char(FeatureNames.Name(ftrs_found(k))) ];
        disp(str_data);
    end
end


%% Run on all data, with best Ilist, and create one model for entire dataset
if 0
    % binary version
    
    % get the top 20 ftrs, and use only these
    [i,a]=sort(counts,'descend'); [i(1:20)',a(1:20)']
    
    char(FeatureNames.Name(a(1:20)))
    
    Amatrix_reduced = Tbl_rangei{:,a(1:20)};    % mat format
    Tbl_rangei_reduced =  Tbl_rangei(:,a(1:20));% table format
    [Ilist,PtPerf_train,PtPerf_test,SVMModel_train,SVMModel_test]=svmforwardselect_perpatient...
        (Amatrix_reduced,criteriaR,tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,...
        forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights,[],[],[]);
    
    % what does the AWP look like as we add features
    PtPerf = PtPerf_train;
    awp_per_ft = (PtPerf(:,3)+PtPerf(:,4)) / 2;
    figure(200); clf(figure(200));
    stairs(awp_per_ft);
    
    % then make a prediction using the current model, for all the data
    PredTsvm = predict(SVMModel_test,Amatrix_reduced);
    performance = PredictiveValue(criteriaR,PredTAll_,Yvariable);
    
    figure(); scatter(PredTsvm, criteriaR);
    [C,order] = confusionmat(criteriaR,PredTsvm); C
    
    plotconfusion(criteriaR',PredTsvm')
    % plotconfusion(targets, outputs) %needs neuralnet tools, data must be in columns not rows
    
    [tbl,chi2,p,lbl] = crosstab(criteriaR',PredTsvm');
    tbl % note rows are criteriaR, columns are PredT
    
else
    %% ====================================================================
    %                  Regression Version
    % =====================================================================
    % Regression version
    
    % get the frequent flyer features i.e. those that occur most often in full selection
    centers = 1:1:size(Amatrix,2);
    Ilist_linear = Ilist_( : ); % convert the Ilist matrix to one linear list
    counts = hist(Ilist_linear,centers); % get the the histogram data
    [i,a]=sort(counts,'descend'); % find the most frequently occuring
    FtrsFreq = a(1:10)'; % limit finding to 10 % ToDo: change to threshold
    
    % get the podium finishing features i.e. those in the first steps of fwd selection
    N_ftrs = 3; % how many steps to take, at least three, probably no more than five
    Ilist_1toN=Ilist_(:,[1:N_ftrs]);
    FtrsFirstN = unique(Ilist_1toN(:),'stable');
    
    % compile one list of these features (as the union of the two sets)
    Ftrs_union = union(FtrsFreq,FtrsFirstN); FeatureNames.Name(Ftrs_union)
    
    tic
    SVMModel1 = fitcsvm(Amatrix(:,Ftrs_union),criteriaR,'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelfunction,'KernelScale',sigma, 'Weights', weights);
    toc % about 5 minutes with 10 - 13 features and 80k observations
    
    % then make a prediction using the current model, for all the data
    tic
    PredTsvm = predict(SVMModel1,Amatrix(:,Ftrs_union));
    toc
    
    % Rsquared can be obtained directly from the model  
    SVMModel1.Rsquared.Ordinary
    SVMModel1.Rsquared.Adjusted
    
    tic
    svpmethodstr = 'SMO';
    kernelf = 'rbf';%'polynomial'
    sigma=2;
    criteriaR_ = Gtest_continuous;
    % can use I to reduce observations
    % 10k observations, 10 secs train, 5 seconds predict
    % 76k observations, 540 secs train, 40 seconds predict
    SVMModel_reg = fitrsvm(Amatrix(:,Ftrs_union),criteriaR_(:),'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelf,'KernelScale',sigma,'Weights',weights(:));
    toc
    
    % can get R-squared direclty from model
    SVMModel_reg.Rsquared.Ordinary
    SVMModel_reg.Rsquared.Adjusted
    
    tic
    PredTsvm_reg = predict(SVMModel_reg,Amatrix(:,Ftrs_union));
    toc
    figure(1); clf(figure(1));
    scatter(criteriaR_, PredTsvm_reg, 1, [0.5 0.05 0.05]);
    xlabel('VE:VDrive'); ylabel('Prediction');
    
    % why two lines appearing in plot?
    
    %figure(); scatter(PredTsvm_reg, Amatrix(:,Ftrs_union(2)), 1, [0.5 0.05 0.05]);
    [RHO,Pval] = corr(criteriaR, PredTsvm_reg, 'type', 'Pearson', 'rows', 'complete');
    RHO^2
    
    % check resdiuals (Residual = observed - fitted)
    % plot of Residual (y) vs Fitted value (x)
    
end

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
PTPerfData=NaN(length(PerPTPerformance),9,size(PerPTPerformance{1},1));
for ftr=1:size(PerPTPerformance{1},1)
    for pt=1:length(PerPTPerformance)
        if isempty(PerPTPerformance{1,pt})
            continue
        end
        PTPerfData(pt,1,ftr) = PerPTPerformance{1,pt}(ftr,1); %PPV
        PTPerfData(pt,2,ftr) = PerPTPerformance{1,pt}(ftr,2); %NPV
        PTPerfData(pt,3,ftr) = PerPTPerformance{1,pt}(ftr,3); %Sens
        PTPerfData(pt,4,ftr) = PerPTPerformance{1,pt}(ftr,4); %Spec
        PTPerfData(pt,5,ftr) = PerPTPerformance{1,pt}(ftr,5); %Acc
        PTPerfData(pt,6,ftr) = PerPTPerformance{1,pt}(ftr,6); %TP
        PTPerfData(pt,7,ftr) = PerPTPerformance{1,pt}(ftr,7); %FP
        PTPerfData(pt,8,ftr) = PerPTPerformance{1,pt}(ftr,8); %TN
        PTPerfData(pt,9,ftr) = PerPTPerformance{1,pt}(ftr,9); %FN
    end
end

if 0 % if only getting data at one specific feature set, it can be a table
    VarNames = {'PPV','NPV','Sens','Spec','Acc','TP','FP','TN','FN'};
    PTPerfTable = table(PTPerfData(:,1),PTPerfData(:,2), PTPerfData(:,3),...
        PTPerfData(:,4),PTPerfData(:,5),PTPerfData(:,6),...
        PTPerfData(:,7),PTPerfData(:,8),PTPerfData(:,9),...
        'VariableNames',VarNames);
    
    
    figure(201); clf(figure(201));
    %ax201(1)=subplot(3,1,1);
    %stairs(table2array(PTPerfTable(:,1:5)));
    %legend('PPV','NPV','Sens','Spec','Acc');
    
    stairs(table2array(PTPerfTable(:,5))); hold on;
    refline(0,median(table2array(PTPerfTable(:,5))));
    title('Performance per pt');
    
    ylabel('Acc');
    xlabel('Pt #');
    
end

% this makes a matrix of average weighted performance (as (sens+spec)/2)
% each column is awp as each feature is added, each row is patient result
% sensitivity is column 3 in PTPerfData, specificity is column 4
awp_each = NaN(size(PTPerfData,1),size(PTPerfData,3));
for ftr=1:size(awp_each,2)
    awp_each(:,ftr) = (PTPerfData(:,3,ftr)+PTPerfData(:,4,ftr))/2;
end

figure(201); clf(figure(201));
stairs(mean(awp_each)); hold on;
%refline(0,median(table2array(PTPerfTable(:,5))));
suptitle('Average weighted performance vs number of features');
title('AWP = mean (each pt ((sens+spec)/2) for ftr n)');
ylabel('AWP'); xlabel('Ftrs');
TidyPlot();

[m,i] = max(mean(awp_each));
x = i; y = m; % Max text position
txt = ([{['Max AWP = ', num2str(m)]},{['at ', num2str(i) ' Ftrs']}]);
text(x,y,txt);

x = 5; y = 0.875; % list position
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


%% svm 3d plot function
% start with the most frequently found features
%Ilist_mode = mode(Ilist_);
%Ilist_mode = Ilist_mode(1,[1:3]);
Ilist_mode = [23, 61, 66]; %
Ilist_mode = [80, 1, 7]; %
Ilist_mode = [69, 59, 7]; % set for SA_OA_20
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
    [IlistAll,PredTAlldata,tempXX2]=svmforwardselect_perpatient(Amatrix,criteriaR,...
        tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights);
    performanceAll = PredictiveValue(criteriaR,PredTAlldata,Yvariable);
end


%%  multiclass classification
% regression model
ftr_set = [23,36,65];
ftr_set = [69,59,7];
tic
SVMModel_pol = fitrsvm(Amatrix(:,ftr_set),criteriaR,'Standardize',true,'Solver',svpmethodstr,...
    'KernelFunction','polynomial','KernelScale',sigma, 'Weights', weights);%'OptimizeHyperparameters','auto'
disp('polynomial time'); %180 secs
toc
tic
SVMModel_lin = fitrsvm(Amatrix(:,ftr_set),criteriaR,'Standardize',true,'Solver',svpmethodstr,...
    'KernelFunction','linear','KernelScale',sigma, 'Weights', weights);%'OptimizeHyperparameters','auto'
disp('linear time'); % 15 secs
toc

PredTsvm_pol = predict(SVMModel_pol,Amatrix(:,ftr_set));
[thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(criteriaR,PredTsvm_pol,1);
PredTsvm_lin = predict(SVMModel_lin,Amatrix(:,ftr_set));
[thresX_,AUC_,SEM_,p_,posclass_,sensspec_]=ROCAUCSEM(criteriaR,PredTsvm_lin,1);

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

% multiple linear regression - regress
Ilist_mode = mode(Ilist_);
Ilist_mode = Ilist_mode(1,[1:3]);
%y = single(Gtest);
y = single(Gtest_continuous);
x1 = Amatrix(:,Ilist_mode(1));
x2 = Amatrix(:,Ilist_mode(2));
%x3 = Amatrix(:,Ilist_mode(3));
X = [ones(size(Amatrix(:,Ilist_mode(1)))) x1 x2 x1.*x2]; %  x3]; %
b = regress(y,X);

figure(5);clf(figure(5));
scatter3(x1(Gtest),x2(Gtest),y(Gtest),1,[0.5 0.05 0.05],'filled'); hold on;
scatter3(x1(~Gtest),x2(~Gtest),y(~Gtest),1,[0.05 0.05 0.05],'filled');
x1fit = min(x1):0.01:max(x1);
x2fit = min(x2):0.01:max(x2);
%x3fit = min(x3):0.01:max(x3);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('ftr 69')
ylabel('ftr 57')
zlabel('drive')
%view(50,10)

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

%% =======================================================================
%                  LDA, Xval with leave one out, fixed feature set
% ========================================================================
% Linear Discriminant Analysis, with Gtest_continuous
% with full cross validation, i.e. leave one out, cycle through all pts
% and no feature selection, i.e. just use the best feature set, Ftrs_union
if 0
   PtData = PtDataBackup;
   FtrsData = FtrsDataBackup;
end
%

DriveThreshold_upper = 0.9;
DriveThreshold_lower = 0.7;
Gtest_continuous = PtData.g_Edi_Adj;
Gtest = ones(size(PtData,1),1)*2; % set a vector of two's
Gtest(Gtest_continuous>=DriveThreshold_upper) = 1; % NFL
Gtest(Gtest_continuous<DriveThreshold_lower) = 0; % FL
%nnz(Gtest==2)
exclude = Gtest==2;
Gtest_continuous(exclude)=[];
Gtest(exclude)=[];
PtData(exclude,:)=[];
Amatrix(exclude,:)=[];
nnz(exclude)
str = [num2str(nnz(~exclude)), ' breaths were labelled as either FL or NonFL']; disp(str)
str = ['Breakdown: ', num2str(nnz(Gtest==1)), ' FL, and ', num2str(nnz(Gtest==0)), ' NonFL']; disp(str)
str = [num2str(nnz(exclude)), ' breaths were unlabelled and removed']; disp(str)

Yvariable = Gtest_continuous;
criteriaR = Gtest;

LDA_results = cell(length(pts),1);
LDA_results_2 = cell(length(pts),1);
count=0;
for pt_in=pts
    count=count+1; %disp(' ');
    str=(['Processing patient ID ', num2str(pt_in), '; number ', num2str(count),' of ', num2str(length(pts))]); disp(str);
    
    % set up the cross validation range
    rangevalidate=pt_in==PtData{:,1};
    
    % run a basic cross validation process, keeping performance for each pt
    
    % train on all data, less one patient, being the test data
    Mdl = fitcdiscr(Amatrix(~rangevalidate,Ftrs_union),criteriaR(~rangevalidate));
    
    % make a prediction on the test patient
    LDA_MDL = predict(Mdl, Amatrix(rangevalidate,Ftrs_union));
    
    x = single(criteriaR(rangevalidate)');
    y = single(LDA_MDL');
    
    LDA_results{count,1} = confusionmat(x,y);
    LDA_results_2{count,1} = PredictiveValue(criteriaR(rangevalidate)',LDA_MDL',Yvariable(rangevalidate)');% 
end

LDA_awp = NaN(size(LDA_results,1),10);
for pt=1:size(LDA_awp,1)
    LDA_awp(pt,1) = LDA_results_2{pt,1}.TP_FP_TN_FN(1);
    LDA_awp(pt,2) = LDA_results_2{pt,1}.TP_FP_TN_FN(2);
    LDA_awp(pt,3) = LDA_results_2{pt,1}.TP_FP_TN_FN(3);
    LDA_awp(pt,4) = LDA_results_2{pt,1}.TP_FP_TN_FN(4);
    LDA_awp(pt,5) = LDA_results_2{pt,1}.PPV_sem_chance_p(1);
    LDA_awp(pt,6) = LDA_results_2{pt,1}.NPV_sem_chance_p(1);
    LDA_awp(pt,7) = LDA_results_2{pt,1}.Sens_sem_chance_p(1);
    LDA_awp(pt,8) = LDA_results_2{pt,1}.Spec_sem_chance_p(1);
    LDA_awp(pt,9) = LDA_results_2{pt,1}.Acc_sem_chance_p(1);
    LDA_awp(pt,10) = (LDA_results_2{pt,1}.Sens_sem_chance_p(1) + LDA_results_2{pt,1}.Spec_sem_chance_p(1))/2;
end

figure(1); clf(figure(1));
stairs(LDA_awp(:,10));
xlabel('Patient number');
ylabel('(Sens+Spec)/2'); 

LDA_summary = array2table(LDA_awp, ...
        'VariableNames', {'TP','FP','TN','FN','PPV','NPV','Sens','Spec','Acc','AWP'});

%% =======================================================================
% LDA, Xval with leave one out, test features within set, no feature selection
% ========================================================================
% Linear Discriminant Analysis, with Gtest_continuous
% with full cross validation, i.e. leave one out, cycle through all pts
% and no feature selection, i.e. run on indiv features in set, Ftrs_union
if 0
   PtData = PtDataBackup;
   FtrsData = FtrsDataBackup;
end

% set up the data for this run
PtData = PtData_Ttran;
Amatrix = FtrsData_Ttran{:,:};
Ftrs_union = [1:size(Amatrix,2)];

DriveThreshold_upper = 0.9;
DriveThreshold_lower = 0.7;
Gtest_continuous = PtData.g_Edi_Adj;
Gtest = ones(size(PtData,1),1)*2; % set a vector of two's
Gtest(Gtest_continuous>=DriveThreshold_upper) = 1; % NFL
Gtest(Gtest_continuous<DriveThreshold_lower) = 0; % FL
table(nnz(Gtest==1),nnz(Gtest==0), nnz(Gtest==2),  'VariableNames', {'NFL','FL','Exclude'})
exclude = Gtest==2;
Gtest_continuous(exclude)=[];
Gtest(exclude)=[];
PtData(exclude,:)=[];
Amatrix(exclude,:)=[];

Yvariable = Gtest_continuous;
criteriaR = Gtest;

LDA_results_each = cell(length(pts),length(Ftrs_union));
LDA_results_2_each = cell(length(pts),length(Ftrs_union));
count=0;
for pt_in=pts
    count=count+1; %disp(' ');
    str=(['Processing patient ID ', num2str(pt_in), '; number ', num2str(count),' of ', num2str(length(pts))]); disp(str);
       
    % set up the cross validation range
    rangevalidate=pt_in==PtData{:,1};
    
    % run a basic cross validation process, keeping performance for each pt
    for ftrNum=1:length(Ftrs_union)
        
        % train on all data, less one patient, that being the test data
        Mdl = fitcdiscr(Amatrix(~rangevalidate,Ftrs_union(ftrNum)),criteriaR(~rangevalidate));
        
        % then make a prediction on the test patient
        LDA_MDL = predict(Mdl, Amatrix(rangevalidate,Ftrs_union(ftrNum)));
        
        x = single(criteriaR(rangevalidate)');
        y = single(LDA_MDL');
        
        LDA_results_each{count,ftrNum} = confusionmat(x,y);
        LDA_results_2_each{count,ftrNum} = PredictiveValue(criteriaR(rangevalidate)',LDA_MDL',Yvariable(rangevalidate)');%
    end
end

LDA_awp_each = NaN(size(LDA_results_2_each,1),length(Ftrs_union));
for pt=1:size(LDA_awp_each,1)
    for ftrNum=1:length(Ftrs_union)
        LDA_awp_each(pt,ftrNum) = (LDA_results_2_each{pt,ftrNum}.Sens_sem_chance_p(1) + LDA_results_2_each{pt,ftrNum}.Spec_sem_chance_p(1))/2;
    end
end

LDA_sens_each = NaN(size(LDA_results_2_each,1),length(Ftrs_union));
for pt=1:size(LDA_sens_each,1)
    for ftrNum=1:length(Ftrs_union)
        LDA_sens_each(pt,ftrNum) = LDA_results_2_each{pt,ftrNum}.Sens_sem_chance_p(1);
    end
end

LDA_spec_each = NaN(size(LDA_results_2_each,1),length(Ftrs_union));
for pt=1:size(LDA_spec_each,1)
    for ftrNum=1:length(Ftrs_union)
        LDA_spec_each(pt,ftrNum) = LDA_results_2_each{pt,ftrNum}.Spec_sem_chance_p(1);
    end
end

LDA_PPV_each = NaN(size(LDA_results_2_each,1),length(Ftrs_union));
for pt=1:size(LDA_PPV_each,1)
    for ftrNum=1:length(Ftrs_union)
        LDA_PPV_each(pt,ftrNum) = LDA_results_2_each{pt,ftrNum}.PPV_sem_chance_p(1);
    end
end

LDA_NPV_each = NaN(size(LDA_results_2_each,1),length(Ftrs_union));
for pt=1:size(LDA_NPV_each,1)
    for ftrNum=1:length(Ftrs_union)
        LDA_NPV_each(pt,ftrNum) = LDA_results_2_each{pt,ftrNum}.NPV_sem_chance_p(1);
    end
end

LDA_Acc_each = NaN(size(LDA_results_2_each,1),length(Ftrs_union));
for pt=1:size(LDA_Acc_each,1)
    for ftrNum=1:length(Ftrs_union)
        LDA_Acc_each(pt,ftrNum) = LDA_results_2_each{pt,ftrNum}.Acc_sem_chance_p(1);
    end
end

LDA_awp_each_mean = nanmean(LDA_awp_each); find(LDA_awp_each_mean>0.675)
LDA_sens_each_mean = nanmean(LDA_sens_each);
LDA_spec_each_mean = nanmean(LDA_spec_each);
LDA_PPV_each_mean = nanmean(LDA_PPV_each);
LDA_NPV_each_mean = nanmean(LDA_NPV_each);
LDA_Acc_each_mean = nanmean(LDA_Acc_each);

figure(1); clf(figure(1));
stairs(LDA_awp_each_mean); hold on;
stairs(LDA_sens_each_mean);
stairs(LDA_spec_each_mean);
stairs(LDA_PPV_each_mean);
stairs(LDA_NPV_each_mean);
stairs(LDA_Acc_each_mean);
xlabel('Ftr number');
ylabel('Average Xval perfromance'); 
legend('AWP','Sens','Spec','PPV','NPV','Acc');
   