function [Gtest_pt, PredY_pt, ManScore_pt, match_ok] = ImportDOSscoring(pt,PtData,Gtest_All,predy)

%% load the drive scoring for current pt
%pt = 5;
match_ok = 1;
DOSDataDir = 'C:\Users\uqdmann\Dropbox\DOS_Scoring_2\DriveScoring\';
DoS = load([DOSDataDir,'DriveScoring_Pt_', num2str(pt), '_TG.mat']);

%% match up the DOS scoring with current pt data

% PtData is all pts, so just get the pt of interest
Isubj = PtData.PT == pt;      % get the pt BB's
Gtest_pt = Gtest_All(Isubj);    % actual
PredY_pt = predy(Isubj);        % pred

PtData_BB_Time = PtData.BB_time(Isubj); % BB time
PtData_BB_VE = PtData.VE(Isubj); % BB VE

DOSData_BB_Time = DoS.BB_Times_All;
DOSData_BB_VE = DoS.VE_All;

%% find matching times (within given tolerance). ver 2
% BB_Times_All are the breaths in the DOS data
% BB_times are the breaths if feature data
tolerance = 0.5; % look for breaths within a given tolerance range
BB_tolerance = @(T) [T-tolerance,T+tolerance];

% step through feature time, looking for tolerable matches in  DOS time
matchedBB_ = NaN(length(PtData_BB_Time),1);
for bb=1:length(PtData_BB_Time)
    [range]=BB_tolerance(PtData_BB_Time(bb)); % get range
    [candidates] = find(DOSData_BB_Time>range(1) & DOSData_BB_Time<range(2));
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

str = [num2str(size(DOSData_BB_Time,1)), ' breaths in DOS data']; disp(str);
str = [num2str(size(PtData_BB_Time,1)), ' breaths in feature data']; disp(str);
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

% test by showing that VE should be VERY similar (only when both are based on real flow)
checkVEmismatch = 1;

if checkVEmismatch
    VE_All_=DOSData_BB_VE(matches(:,1));
    VE_ = PtData_BB_VE(matches(:,2));
    nearZero = nnz(abs(round(VE_All_,4) - round(VE_,4))>0.01);
    if nearZero > 0
        disp('WARNING - VE values do not align for DOS and Feature data (you should check matches)');
        match_ok = 0;
        if 0
            odds = find(abs(round(VE_All_,4) - round(VE_,4))>0.001==1);% find the offending items
            round(VE_All_(odds),4) - round(VE_(odds),4) % what are their values
        end
    end
    
    if 0
        figure(10); clf(figure(10));
        plot(PtData_BB_Time(matches(:,2),:), 1, 'k.'); hold on;
        plot(DOSData_BB_Time(matches(:,1),:), 2,'k.');
        ylim([0.1 2.9]);
    end
end

%% generate tidy output
Gtest_pt = Gtest_pt(matches(:,2));    % actual
PredY_pt = PredY_pt(matches(:,2));    % predicted

% and now, the most important bit, the human score
% col 1 is estimate, col 2 is human score
ManScore_pt = DoS.BB_DOS(matches(:,1),2); 
% value 0 is FL, 2 is intermediate, 4 is NotFL, 5 is shitty data
ManScore_pt(ManScore_pt==0) = 0.1;
ManScore_pt(ManScore_pt==2) = 0.5;
ManScore_pt(ManScore_pt==4) = 1.0;
ManScore_pt(ManScore_pt==5) = NaN;

end
