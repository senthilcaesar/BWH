clear all

for firstrun = 1:-1:0


if firstrun
Timing = 'O';
load('FS_FD_5.mat');
else
Timing = 'T';
load('FS_FD_6.mat');    
end

DriveAnalysis=1;
ShowFigures=0;
ShowFlowFigures=0;

n_pts=length(SleepData); % how many pts are we processing
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
		% and all the other breath data that is later matched up
        %  (BB_time,BB_Ttot,Ar,win,notAr,Veup,hypnog,pcw,Evt,FL)
        [BB_time,VdriveEdi,VdrivePes,VE,Ar,win,notAr,Veup,hypnog,pcw,BB_Ttot,FL,EType,ApneaB] = ...
            VE_VdriveFLArrayQAO(pt,criteria,BreathFLDataTable, BreathDataTable);
        
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

        if ShowFlowFigures
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
            'g_Edi_Adj','g_Pes_Adj','ApneaB','BB_time','BB_Ttot'};
        PT = pt*ones(length(VE),1);
        PtDataWin = table(PT, win, Ar, notAr, a1, a2, hypnog, EType, Veup, VE, VdriveEdi, VdrivePes,...
            g_Edi_Adj, g_Pes_Adj, ApneaB, BB_time, BB_Ttot, 'VariableNames', VarNames);

        % add the Flow limited table to the current data table
        PtDataWinFL = [PtDataWin,FL];

        % add this window data with the FL data, to a growing combined table
        CombinedPtData = cat(1,CombinedPtData,PtDataWinFL);

    catch me
        disp(me.message)
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
clearvars -except PtData n_pts DataSet savestr_SVM DriveAnalysis firstrun
%clear LocalSignals

%% remove apnea breaths
ApneaFlow = (PtData.ApneaB==1);
str = ['Removing ', num2str(nnz(ApneaFlow)), ' apneic breaths']; disp(str);
PtData(ApneaFlow,:)=[];

%% remove breaths with < threshold% flow below eupnoea
% where VE < 10% of Veup
LowFlow = PtData.VE < (PtData.Veup*0.1);
str = ['Removing an additional ', num2str(nnz(LowFlow)), ' low flow breaths (<10% of eupnea)']; disp(str);
PtData(LowFlow,:)=[];

%% Remove extreme outliers: Part 1
% for e.g. VTi_VTe should be around 1, so values of 100 are clearly abnormal. 
firstdatacol = 18; % confirm this during procesing
if 0
    %original, very simplistic method
    centileupper = prctile(PtData{:,firstdatacol:end},99.9);
    centilelower = prctile(PtData{:,firstdatacol:end},0.1);
    
    %DM comment: in the matlab documentation, outlier detection using
    % quantiles is F1=Q1-1.5*IQR, and F2=Q3+1.5*IQR, with F1 > outliers > F2
    IQR_ = iqr(PtData{:,firstdatacol:end});
    Q1_ = quantile(PtData{:,firstdatacol:end}, 0.25);
    Q3_ = quantile(PtData{:,firstdatacol:end}, 0.75);
    prctilefactor=1.5; %e.g. set to about 5 to be equivalent to Scotty's method
    centileupper_ = Q3_ + prctilefactor*IQR_;
    centilelower_ = Q1_ - prctilefactor*IQR_;
else %SS comment: SS revised version labels outliers e.g. if 2 fold above the 95-50 difference, e.g. median = 5, 95th centile = 10; outliers >15.
    prctilethres=99; %e.g. 95
    prctilefactor=3; %e.g. 2
    centileupper = prctile(PtData{:,firstdatacol:end},50) + (prctile(PtData{:,firstdatacol:end},prctilethres)-prctile(PtData{:,firstdatacol:end},50))*prctilefactor;
    centilelower = prctile(PtData{:,firstdatacol:end},50) - (prctile(PtData{:,firstdatacol:end},50)-prctile(PtData{:,firstdatacol:end},100-prctilethres))*prctilefactor;
end
outliers = (PtData{:,firstdatacol:end} > centileupper | ...
    PtData{:,firstdatacol:end} < centilelower);
upper_outliers = nnz(PtData{:,firstdatacol:end} > centileupper)
lower_outliers = nnz(PtData{:,firstdatacol:end} < centilelower)
filler = false(size(PtData,1),firstdatacol-1); % exclude the first 15 columns of PtData
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
else %SS version, replaces with 
    for n=1:length(centileupper)
        PtData{:,n+firstdatacol-1}(PtData{:,n+firstdatacol-1} > centileupper(n)) = centileupper(n);
        PtData{:,n+firstdatacol-1}(PtData{:,n+firstdatacol-1} < centilelower(n)) = centilelower(n);
    end     
end

%% Remove features

PtData.AA_NumOfPeaks=[];
PtData.AA_Pf_1_3=[];
PtData.AA_Pf_2_3=[];
PtData.AA_Pf_3_3=[];
PtData.AA_IsTerminalPeak=[];
PtData.AA_PeaksRatio=[];
PtData.TTran_i=[];
PtData.TTran_e=[];
PtData.MansourRUA=[];
PtData.TpeakE_Te=[]; %same as Morris, exactly

%%
if firstrun

PtDataO=PtData;
save temp PtDataO
else

PtDataT=PtData;
load temp PtDataO
end

end
%% Merge Tables

FeatureNamesO = PtDataO.Properties.VariableNames(firstdatacol:end);
AmatrixO = table2array(PtDataO(:,firstdatacol:end));

for i=1:length(FeatureNamesO)
    FeatureNamesO{i}=[FeatureNamesO{i} '_O'];
end

FeatureNamesT = PtDataT.Properties.VariableNames(firstdatacol:end);
AmatrixT = table2array(PtDataT(:,firstdatacol:end));

for i=1:length(FeatureNamesT)
    FeatureNamesT{i}=[FeatureNamesT{i} '_T'];
end

FeatureNames = [FeatureNamesO FeatureNamesT];
Amatrix = [AmatrixO AmatrixT];
%%
allnanrows = sum(isnan(Amatrix),2)==size(Amatrix,2);
%sum(allnanrows)
Amatrix(allnanrows,:)=[];
PtData(allnanrows,:)=[];

Fnan=sum(isnan(Amatrix))/size(Amatrix,1);
I=Fnan>0.001;
Amatrix(:,I)=[];
FeatureNames(:,I)=[];

%%
Name = FeatureNames';
Ftr = (1:length(FeatureNames))';
FeatureNames=table(Ftr,Name);
%%
clearvars -except PtData Amatrix FeatureNames
save SSblend_FS_FD_56
