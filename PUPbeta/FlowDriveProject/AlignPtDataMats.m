%% This script aligns rows in original and transition time feature spaces
% This can be applied to the large feature spaces ...
%

%% startup
close all
clear global
clear
clc

%% load the two PtData files
LocalDataDir = 'C:\PSG_Data\QAO';% path to data
LocalData = [LocalDataDir,'\FS_FD_9_PtData.mat']; % load original timing
load(LocalData,'PtData'); % selectively load variables of interest
PtData_Orig = PtData;
clearvars -except PtData_Orig LocalDataDir
LocalData = [LocalDataDir,'\FS_FD_10_PtData.mat']; % load transition timing
load(LocalData,'PtData'); % selectively load variables of interest
PtData_Tran = PtData;
clearvars -except PtData_Orig PtData_Tran LocalDataDir

%% first pass, just work on limited columns
% cols 15:17 takes 1hr 5 mins to run on ~83k rows
% cols 16:17 will take about an hour to run on ~83k rows
% col 17 alone takes about 10 minutes, but is far too innacurate
progressbar('Progress');
t_start = clock;
% find each row from shorter table, in the longer table
% Orig is shorter than Ttran - Ttran is more breaths?
AlignData = [];
for m = 1:size(PtData_Orig,1)
    progressbar(m/size(PtData_Orig,1)); 
    Dat1 = PtData_Orig{m,15:17};
    if (m>=51); k=m-50; else; k=1; end
    for n = k:size(PtData_Tran,1)
        Dat2 = PtData_Tran{n,15:17};
        mismatches = nnz(diff([Dat1; Dat2]));
        if mismatches == 0 % perfect match
            AlignData = [AlignData; [m,n]];
            break 
        end
        if mismatches > 4 % way off, just continue
            continue;
        else % pretty close, check if it's nan matching 
            if isequaln(Dat1, Dat2)
                AlignData = [AlignData; [m,n]];
            break
            end
        end
    end
end
progressbar(1); % force close the progress bar
t_end = clock;
delta_t = etime(t_end, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
str = ['Total processing time: ', char(D), ' (hh:mm:ss)']; disp(str);

%% second pass, work on extended range, using AlignData created above
% this will take about 11 minutes to run on ~83k rows.
% using cols 1 to 12, excluding g_adj_ cols (13, 14) as these can change
tic
progressbar('Progress');
ConfirmData = false(size(AlignData,1),1);
for m = 1:size(AlignData,1)
    progressbar(m/size(AlignData,1));
    i1 = AlignData(m,1);
    Dat1 = PtData_Orig{i1,1:12};
    i2 = AlignData(m,2);
    Dat2 = PtData_Tran{i2,1:12};
    
    mismatches = nnz(diff([Dat1; Dat2]));
    if mismatches == 0 % perfect match
        ConfirmData(m) = true;
        continue;
    end
    if mismatches > 4 % way off, just continue
        ConfirmData(m) = false;
        continue;
    else % pretty close, check if it's nan matching
        if isequaln(Dat1, Dat2)
            ConfirmData(m) = true;
            continue;
        end
    end
end
progressbar(1); % force close the progress bar
toc

%% let's find any errors, and then have a quick look at them
errors = find(ConfirmData==false);
m = errors;  
i1 = AlignData(m,1);
Dat1 = PtData_Orig(i1,1:17);
i2 = AlignData(m,2);
Dat2 = PtData_Tran(i2,1:17);
Dat = [Dat1;Dat2];

%% create aligned data
PtData_O = PtData_Orig(AlignData(:,1),:);
PtData_T = PtData_Tran(AlignData(:,2),:);

clearvars -except PtData_O PtData_T LocalDataDir

%% save the aligned data
savestr = [LocalDataDir, '\FS_FD_A_9_10_1.mat'];
save(savestr, '-v7')
