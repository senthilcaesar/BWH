function MakeFeatureSpace()
%% This function performs additional analysis on the Flow Drive data
% It is designed to run after Analysis.m. The output from Analysis.m could
% be a single file for each pt, or one large file with each pt in a cell
% array. This script uses the large single file. If the individual files
% were produced, a single large file can be made with the function called
% 'CombineIndivPtAnalysisMats.m'
%
% This procduces a flow drive feature space
%
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
% remove wakeful breaths with FL (<0.7) (only those > 4 breaths from sleep)
%

%% Inputs / Outputs
% Inputs
%   - file, one (potentially very large) file with all patient data
%   - file, the analysis spreadsheet
%
% Outputs
%   - file, Flow Drive feature space table

%% startup
close all
clear global AnalyzeDataSpreadsheet
clear
clc

%% Local variables and settings
if verLessThan('matlab', '9.2')   % 2017a is ver 9.2,
    % i.e. matlab 2016b or anything older
    if 1
        AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_MESA.xlsx';
        %AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_OA.xlsx';
        DriveAnalysis = 0;      % set to one to calc drive ratios, must set as zero in OA data, can be either in FD data
    else
        AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_Pnasal.xlsx';
        DriveAnalysis = 1;      % set to one to calc drive ratios, must set as zero in OA data, can be either in FD data
    end
else
    % i.e. matlab 2017a (only)
    AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet_pneumotach.xlsx';
    DriveAnalysis = 1;      % set to one to calc drive ratios, must set as zero in OA data, can be either in FD data
end

manual_ref_breaths = 0; % set to one to use the manually scored reference breaths, zero uses auto found ref breaths
RefWakeOnly = 1;        % in manual_ref_breaths, set to one to use ref breaths from wake only, zero to use all ref breaths
checkVEmismatch = 0;    % in manual_ref_breaths, set to one to check VE mismatch in alignment, zero is req'd when importing pnasal
ShowFigures = 0;
ShowFlowFigures = 0;    % show the flow and edi signals during processing

%% pre-processing
% read spreadsheet (options worksheet)
[~,~,raw] = xlsread(AnalyzeDataSpreadsheet,2,'C3:C32');
settings.savename = char(raw{1});
settings.OutputDataDirectory = char(raw{18});
if settings.OutputDataDirectory(end)~='\'
    settings.OutputDataDirectory=[settings.OutputDataDirectory '\'];
end
%settings.savestr = [settings.OutputDataDirectory, 'TestingOnly_DeleteAtWill'];

% read spreadsheet (files worksheet)
[num,~,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:K300');
settings.pt_indx = num(:,1);
settings.pt_flow  = logical(num(:,3));
settings.pt_pnasal  = logical(num(:,4));
settings.pt_edi  = logical(num(:,5));
settings.pt_pes  = logical(num(:,6));
settings.n_pts = length(settings.pt_indx); % how many pts are we processing

%% declare variables
FeatureSpace = [];
Nwindows=NaN(settings.n_pts,1);
threshold=NaN(settings.n_pts,1);
BreathsPerPatient=NaN(settings.n_pts,5);

%% processing
try
    filename = [settings.OutputDataDirectory, settings.savename, '_All.mat'];
    if exist(filename,'file') == 2
        displaytext=['Loading saved data: ', filename]; disp(displaytext);
        % selectively load variables of interest
        load(filename,'SleepData','LG_QualityInfo','DataOut',...
            'BreathDataTable','BreathFLDataTable','LocalSignals');
        for pt=1:settings.n_pts
            checker = 0;
            try
                checker = isnan(BreathDataTable{pt});
            catch
                %disp('All good'); % if we get here, there is data
            end
            if checker
                disp(' ');
                str = ['Skipping patient ', num2str(pt)]; disp(str);
                continue
            else
                disp(' ');
                str = ['Processing patient ', num2str(pt)]; disp(str);
                supinecodes = [0 2 -5]; %Profusion:[1],Spike:[0 2]
                maxwakethres = 360;  %minwakearthres = 60;
                minNevents = 0;
                maxFREM=0;
                %N_events = LG_QualityInfo{pt}(:,2);
                %Pos = LG_QualityInfo{pt}(:,5);
                %Fwake = SleepData{pt}(:,1);
                %FREM = SleepData{pt}(:,6);
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
                FeatureSize = 170; %was 168; %was 172; % was 160 
                [BB_time,VdriveEdi,VdrivePes,VE,Ar,win,notAr,Veup,hypnog,pcw,BB_Ttot,FL,EType,ApneaB] = ...
                    VE_VdriveFLArray(pt,criteria,BreathFLDataTable, BreathDataTable, FeatureSize);
                
                if ~isempty(BB_time) % only continue processing this pt if data
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
                    a1 = ad>threshold(pt); % threshold # breaths away from sleep, usually 4 can vary, see above
                    a2 = ad>=2; % at least two breaths away from sleep
                    
                    %% aims to determine the best VE/Vdrive ratio for Edi and Pes
                    if DriveAnalysis
                        %% the new method for determining adjusted drive
                        % is done with the manual reference breath scoring
                        if manual_ref_breaths
                            str_automan = ['ManRef2'];
                            clear BB_DOS
                            % load reference scoring for this patient
                            refscoring = ['..\..\DOS_Scoring_2\DriveScoring\DriveScoring_Pt_',num2str(pt),'_Ref_TG.mat'];
                            load(refscoring);                                                  
                            if exist('BB_DOS','var') % check if the file opened
                                %% confirm times match
                                % BB_Times_All are the breaths in the reference scoring
                                % BB_DOS is the 'degree of shitness' score,
                                % with column 2 being the human score, and
                                % a value of 4 being a good breath.
                                % VE_All is the ventilation for each breath
                                % 
                                
                                if 0 % excluding this as it adds nothing
                                % handle dropped breaths through indexing by 'ib'
                                [matchedBB, ia_, ib] = intersect(BB_Times_All, BB_time);
                                hypnogBB = hypnog(ib);
 
                                if size(matchedBB,1) ~= size(BB_time,1)
                                    figure(10); clf(figure(10));
                                    plot(BB_time, 1, 'k.'); hold on;
                                    plot(BB_Times_All, 2,'k.');
                                    plot(BB_time(ib), 3, 'k.'); 
                                    ylim([0.1 3.9]);
                                    % may need to zoom in to spot diffs
                                end
                                end
                                
                                %% find matching times. ver 2
                                % BB_Times_All are the breaths in the DOS data
                                % BB_times are the breaths if feature data
                                tolerance = 0.5; % look for breaths within a given tolerance range
                                BB_tolerance = @(T) [T-tolerance,T+tolerance];
                                
                                % step through feature time, looking for tolerable matches in  DOS time
                                matchedBB_ = NaN(length(BB_time),1);
                                for bb=1:length(BB_time)
                                    [range]=BB_tolerance(BB_time(bb)); % get range
                                    [candidates] = find(BB_Times_All>range(1) & BB_Times_All<range(2));
                                    switch size(candidates,1)
                                        case 0
                                            %str=['No match for breath ', num2str(bb)]; disp(str);
                                            matchedBB_(bb) = NaN;
                                        case 1
                                            matchedBB_(bb) = candidates;
                                        otherwise % more than one
                                            str=['Multiple matches for breath ', num2str(bb)]; disp(str);
                                    end
                                end
                                str = [num2str(size(BB_DOS,1)), ' breaths in reference data']; disp(str);
                                str = [num2str(size(BB_time,1)), ' breaths in feature data']; disp(str);
                                str = [num2str(nnz(~isnan(matchedBB_))), ' breaths time matched']; disp(str);

                                %% set up indices into both data
                                matches=[];
                                for bb=1:length(matchedBB_)
                                    row = matchedBB_(bb);
                                    if isnan(row)
                                        continue
                                    else
                                        matches=[matches;[row,bb]]; %  Flow, Pnasal
                                    end
                                end
                                
                                % matches(:,1) is the index into DOS data
                                % matches(:,2) is the index into feature data
  
                                % test by showing that VE should be VERY similar (when feature space is real flow)
                                % it will be different when the feature space is Pnasal, because DOS data is real flow
                                if checkVEmismatch
                                    VE_All_=VE_All(matches(:,1));
                                    VE_ = VE(matches(:,2));
                                    nearZero = nnz(abs(round(VE_All_,4) - round(VE_,4))>0.01);
                                    if nearZero > 0
                                        disp('WARNING - VE values do not align for DOS and Feature data (you should check matches)');
                                        keyboard
                                        if 0
                                            odds = find(abs(round(VE_All_,4) - round(VE_,4))>0.001==1);% find the offending items
                                            round(VE_All_(odds),4) - round(VE_(odds),4) % what are their values
                                        end
                                    end
                                    
                                    if ShowFigures&&0
                                        figure(10); clf(figure(10));
                                        plot(BB_time(matches(:,2),:), 1, 'k.'); hold on;
                                        plot(BB_Times_All(matches(:,1),:), 2,'k.');
                                        ylim([0.1 2.9]);
                                    end    
                                end
                                
                                % prepare variables for co-locating ref breaths in both DOS and feature data
                                hypnogBB = hypnog(matches(:,2));
                                Ar_BB = Ar(matches(:,2));
                                BB_DOS_matched = BB_DOS(matches(:,1),:); 
                                
                                %% find VE/Vdrive for all WAKE ONLY reference breaths
                                if RefWakeOnly
                                    if 1 % Wake epoch, and scored Ar within sleep epoch
                                        refBB = find((BB_DOS_matched(:,2)==4) & ((hypnogBB==4)|(Ar_BB==1)));
                                    else
                                        refBB = find((BB_DOS_matched(:,2)==4) & (hypnogBB==4));
                                    end
                                else
                                    refBB = find(BB_DOS_matched(:,2)==4);
                                end
                                
                                % what do the refBB's relate to in the feature data
                                refBBintoftr = matches(refBB,2); % [refBB refBBintoftr]
                                reftimes = BB_time(refBBintoftr);
                                [times_,ia,~] = unique(reftimes);  % control for any duplicates
                                str = [num2str(length(times_)), ' manual reference breaths found']; disp(str);

                                %% edi first
                                data_edi = VE(refBBintoftr)./VdriveEdi(refBBintoftr);
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
                                str = [num2str(length(G_edi_r)), ' Edi reference breaths']; disp(str);
                                
                                %% then pes
                                data_pes = VE(refBBintoftr)./VdrivePes(refBBintoftr);
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
                                str = [num2str(length(G_pes_r)), ' Pes reference breaths']; disp(str);
                                
                            else
                                disp('Reference scoring not found!');
                            end
                        else
                            %% the original method for determining the adjusted drive
                            %  this is done with moving median value over the wake breaths
                            %  that are more than 2 breaths away from sleep
                            str_automan = ['AutoRef2'];
                            
                            % edi first
                            if settings.pt_edi(pt)
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
                                str = [num2str(length(G_edi_w)), ' Edi reference breaths']; disp(str);
                            end
                            
                            % then pes
                            if settings.pt_pes(pt)
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
                                str = [num2str(length(G_pes_w)), ' Pes reference breaths']; disp(str);
                            end
                            a2_locs = find(a2==1);
                            refBBintoftr = a2_locs(ia); % set the index of ref breaths, mainly for plotting
                            
                        end % end of import ref breaths
                        
                        %% Pick the G_edi and G_pes to use hereafter
                        if manual_ref_breaths % using manual reference scoring
                            if settings.pt_edi(pt)
                                G_edi = G_edi_r;
                            end
                            if settings.pt_pes(pt)
                                G_pes = G_pes_r;
                            end
                        else % using automated process based on wake
                            if settings.pt_edi(pt)
                                G_edi = G_edi_w;
                            end
                            if settings.pt_pes(pt)
                                G_pes = G_pes_w;
                            end 
                        end
                        BreathsPerPatient(pt,1) = pt;
                        BreathsPerPatient(pt,3) = size(G_pes,1); % Edi ref breaths
                        
                        %% all breaths are then normalised, i.e. divided by their corresponding VE/Vdrive ratio
                        % Edi first
                        g_Edi_Adj = NaN(length(VE),1);
                        if settings.pt_edi(pt)
                            Gt=interp1(G_edi(:,1),G_edi(:,2),BB_time);
                            Gt(BB_time<G_edi(1,1))=G_edi(1,2);
                            Gt(BB_time>G_edi(end,1))=G_edi(end,2);
                            g_Edi_Adj = (VE./VdriveEdi)./Gt;
                            Gt_edi = Gt;
                        end
                        
                        % then Pes
                        g_Pes_Adj = NaN(length(VE),1);
                        if settings.pt_pes(pt)
                            clear Gt
                            Gt=interp1(G_pes(:,1),G_pes(:,2),BB_time);
                            Gt(BB_time<G_pes(1,1))=G_pes(1,2);
                            Gt(BB_time>G_pes(end,1))=G_pes(end,2);
                            g_Pes_Adj = (VE./VdrivePes)./Gt;
                            Gt_pes = Gt;
                        end
                        
                        %% potentially remove major outliers, i.e. only use 5th to 95th centile
                        %g_Edi_Adj(g_Edi_Adj>prctile(g_Edi_Adj,95))=NaN;
                        %g_Pes_Adj(g_Pes_Adj>prctile(g_Pes_Adj,95))=NaN;
                        
                        str = [num2str(nnz(g_Edi_Adj)), ' breaths with g_Edi']; disp(str);
                        BreathsPerPatient(pt,2) = nnz(g_Edi_Adj); % total breaths
                        
                        %% dividing by the 'best' ratio, we have normalised to 1.
                        % values much >1 should occur infrequently, as this would suggest
                        % VE exceeded Vdrive, as such, they should be capped (x>=1.5 == 1.5 )
                        str = [num2str(nnz(g_Edi_Adj>1.5)), ' breaths capped to g = 1.5']; disp(str);
                        BreathsPerPatient(pt,4) = nnz(g_Edi_Adj>1.5); % breaths capped to 1.5
                        g_Edi_Adj(g_Edi_Adj>1.5)=1.5;
                        g_Pes_Adj(g_Pes_Adj>1.5)=1.5;
                        
                        % a value of 1 would indicate good breathing, and <1 indicates FL.
                        % to help model fitting, and not keeping misleading data, breaths
                        % in wake that are clearly FL (x<0.7) are removed, but..
                        % only if they are >threshold # of breaths away from sleep
                        str = [num2str(nnz(g_Edi_Adj<0.5 & (a1))), ' awake breaths with g < 0.7 removed']; disp(str);
                        BreathsPerPatient(pt,5) = nnz(g_Edi_Adj<0.5 & (a1)); % wake breaths removed
                        g_Edi_Adj(g_Edi_Adj<0.7 & (a1))=NaN; % far away from sleep, must be good
                        g_Pes_Adj(g_Pes_Adj<0.7 & (a1))=NaN;
                        %g_Edi_Adj(g_Edi_Adj<0.5 & (a2))=NaN; % more tolerant when closer to sleep
                        %g_Pes_Adj(g_Pes_Adj<0.5 & (a2))=NaN;
                        
                        % redo the figure for VE/Vdrive ratio with adjusted Edi and Pes
                        % the wake data should be approximately 1
                        % the reference data should be exactly one.
                        if ShowFigures
                            VEVdrive_raw = (VE./VdriveEdi);
                            if manual_ref_breaths
                                titlestr = ['Normalising Pnasal - Pt ', num2str(pt), ' Manually scored ref breaths'];
                            else
                                titlestr = ['Normalising Pnasal - Pt ', num2str(pt), ' Automatically selected ref breaths'];
                            end
                            figure(100); clf(figure(100));
                            ax(1) = subplot(2,1,1);
                            plot(BB_time,VEVdrive_raw,'k.'); hold on; % all breaths
                            plot(BB_time(refBBintoftr),VEVdrive_raw(refBBintoftr),'r.');                          
                            plot(G_edi(:,1),G_edi(:,2),'gd');
                            plot(BB_time, Gt_edi, 'b+', 'Markersize', 2); % quicker to plot, but less accurate...
                            legend('All BB','Ref BB','weighted','normal factor');
                            ylabel('non-normalised gEdi (VE/Vdrive)');
                            title(titlestr);
                            
                            ax(2) = subplot(2,1,2);
                            plot(BB_time, g_Edi_Adj,'k.'); hold on;
                            plot(BB_time(refBBintoftr), g_Edi_Adj(refBBintoftr),'r.');
                            refBB_mean = nanmean(g_Edi_Adj(refBBintoftr));
                            refBB_mean_plot = ones(length(BB_time)).*refBB_mean;
                            plot(BB_time, refBB_mean_plot, 'k:');
                            legend('All BB','Ref BB', 'Ref mean'); 
                            ylabel('normalised gEdi (VE/Vdrive)');
                            linkaxes(ax,'x');
                            
                            fig = gcf;
                            fig.Color = [1 1 1]; % set background colour to white
                            fig.Units = 'inches';
                            fig.Position = [-12.333       1.7708       12.125       8.8021];
                            
                            str = ['..\Figures\', titlestr ];
                            
                            saveas(fig, str, 'png'); %savefig(str);
                        end
                        
                    else
                        %% not doing drive analysis, i.e. just using this for testing and/or training purposes
                        %manually set these as NaN's
                        g_Edi_Adj = NaN(length(BB_time),1);
                        g_Pes_Adj = NaN(length(BB_time),1);
                        str_automan = ['NoRef'];
                    end
                    
                    %% Make a table of all breaths from this window
                    VarNames = {'PT','Window','Ar','NotAr','A1','A2','Hypnog','Etype','Veup','VE','DriveEdi','DrivePes',...
                        'g_Edi_Adj','g_Pes_Adj','ApneaB','BB_time','BB_Ttot'};
                    PT = pt*ones(length(VE),1);
                    PtDataWin = table(PT, win, Ar, notAr, a1, a2, hypnog, EType, Veup, VE, VdriveEdi, VdrivePes,...
                        g_Edi_Adj, g_Pes_Adj, ApneaB, BB_time, BB_Ttot, 'VariableNames', VarNames);
                    
                    % add the Flow limited table to the current data table
                    PtDataWinFL = [PtDataWin,FL];
                    
                    %% check out what this data looks like as an 'image'
                    %NaNSpace = isnan(PtDataWinFL{:,:});
                    %NaNSpaceMap = NaNSpace*90;
                    %figure(pt);image(NaNSpaceMap);
                    if ShowFigures&&0
                        figure(pt); clf(figure(pt));
                        subplot(1,2,1);
                        NaNSpace = isnan(PtDataWinFL{:,:}); NaNSpaceMap = NaNSpace*90;
                        image(NaNSpaceMap); xlabel('Table columns'); ylabel('breath number');
                        title('All breaths');
                        subplot(1,2,2);
                        BreathFLDataTable_NoApnea = PtDataWinFL{:,:};
                        BreathFLDataTable_NoApnea(PtDataWin.ApneaB==1,:) = [];
                        NaNSpace = isnan(BreathFLDataTable_NoApnea); NaNSpaceMap = NaNSpace*90;
                        image(NaNSpaceMap); xlabel('Table columns'); ylabel('breath number');
                        title('Excluding apnoea');
                    end
                    % add this window data with the FL data, to a growing combined table
                    FeatureSpace = cat(1,FeatureSpace,PtDataWinFL);
                    
                else
                    disp(['Pt ', num2str(pt), ' had no usable drive data']);
                end
            end
        end
    else
        displaytext=['No saved Analysis data']; disp(displaytext);
    end
catch MainProcessing
    %     disp(MainProcessing.message);
    disp(MainProcessing.getReport);
end

%% Make a table of breath count summary values
ptlist = unique(FeatureSpace.PT);
BreathsPerPatient = BreathsPerPatient(ptlist,:);
VarNames = {'PT','Total','Ref','Capped','Removed'};
PtDataBBSummary = array2table(BreathsPerPatient, 'VariableNames', VarNames);

%
%% save feature space
disp('Saving Feature Space');
save([settings.OutputDataDirectory, settings.savename,'_FeatureSpace_', str_automan],...
    'FeatureSpace','PtDataBBSummary', '-v7');

%% Unused code - just figures
if 0
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
    
    if ShowFigures
        % determine the windows that we analysed
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
    
end







