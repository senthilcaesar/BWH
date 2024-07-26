function TidyFeatureSpace()
%% This function performs Feature Space cleaning and preparation for run
% This handles NaN features, apnea breaths, low flow breaths, and also does
% some general housekeeping (extreme outliers, etc.)
%
% Inputs
%  a newly created Feature Space, (i.e. the variable called 'FeatureSpace')
%
% Outputs
%  a tidied Feature Space, ready to run in SStestL1Ox
%   PtData - Table of patient and breath data
%   Amatrix - matrix of Features (per col) and breaths (per row)
%   FeatureNames - list of column titles (i.e. Features) in Amatrix
%
% The process is:
%   apnea
%   low flow
%   NaN - features
%   NaN - breaths
%   outliers
%   continuous to discrete transform - i.e. classification
%   formatting output to suit next step

%% start
close all
clear
clc
addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');

%% runs     FlowDriveData   Use25HzDS   EdiDrive
% OA25Hz        0               1           n/a
% OA125Hz       0               0           n/a
% FD25HzEdi     1               1           1
% FD125HzEdi    1               0           1   

%% options
FlowDriveData = 1;      % set as 1 for flow drive data, or 0 to use OA data

RealFlow = 0;           % set as 1 for real flow, 0 for pnasal flow
Use25HzDS = 1;          % set as 1 for 25Hz data, or 0 to use 125Hz data
EdiDrive = 1;           % set as 1 for Edi Drive data, 0 for Pes drive
ShowFigures = 0;        % set as 1 to show figures, or 0 to not show figures
RemoveRedundantFtrs = 1;% set as 1 to find and remove  redundant ftrs, 0 to ignore (and keep in)
firstdatacol = 18;      % confirm this during procesing

%datadir = '..\FeatureSpaces\';  % 
SaveCleanMat = 1; 

%% open file
if FlowDriveData
    datadir = 'C:\PSG_Data\FlowDrive\Analyzed\';
    if Use25HzDS
        if RealFlow
            %filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_AutoRef2.mat'];
            filename = [datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2.mat'];
        else % Pnasal
            %filename = [datadir, 'PnasalDrive_only25Hz_exp067_FeatureSpace_ManRef2.mat']; 
            %filename = [datadir, 'PnasalDrive_only25Hz_exp1_FeatureSpace_AutoRef2.mat'];
            filename = [datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2.mat'];
        end
    else; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace.mat'];
    end
else
    datadir = 'C:\PSG_Data\OralAppliance\Analyzed\';
    if Use25HzDS; filename = [datadir, 'FlowOA_25Hz_ReNormalized_FeatureSpace_NoRef.mat'];                                      
    else; filename = [datadir, 'FlowOA_only125Hz_FeatureSpace.mat'];                                
    end
end
str=['Loading ' filename]; disp(str);
load(filename, 'FeatureSpace');
PtData = FeatureSpace; % keep FeatureSpace as a backup, work with PtData
str = ['Starting with ', num2str(height(PtData)), ' breaths']; disp(str);

%% remove labelled apnea breaths (as labelled by VEfromFlow)
ApneaFlow = (PtData.ApneaB==1);
pts = unique(PtData.PT(ApneaFlow));
for pt = 1:length(pts)
    current_pt = PtData.PT==pts(pt);
    temp_data = PtData(ApneaFlow & current_pt,:);
    RemovedBB_Ap(pt, 1) = pts(pt);
    RemovedBB_Ap(pt, 2) = height(temp_data);
    RemovedBB_Ap(pt, 3) = sum(temp_data.BB_Ttot);
    RemovedBB_Ap(pt, 4) = nnz(temp_data.Etype==2);% how many of these were clinically scored as (2) obstructive ap
    RemovedBB_Ap(pt, 5) = nnz(temp_data.Etype==3);% 3 central ap
    RemovedBB_Ap(pt, 6) = nnz(temp_data.Etype==4);% 4 hypopnoea ob
    RemovedBB_Ap(pt, 7) = nnz(temp_data.Etype==5);% 5 mixed
end
RemovedBB_Apnoea = array2table(RemovedBB_Ap, 'VariableNames', {'Pt','NumBB','TotalTime','ApO','ApC','HypO','Mixed'})
str = ['Removing ', num2str(nnz(ApneaFlow)), ' apneic breaths']; disp(str);
PtData(ApneaFlow,:)=[];

%% remove breaths with < threshold% flow below eupnoea
% where VE < 10% of Veup
LowFlow = PtData.VE < (PtData.Veup*0.1); %nnz(PtData.VE < (PtData.Veup*0.1))
pts = unique(PtData.PT(LowFlow));
for pt = 1:length(pts)
    current_pt = PtData.PT==pts(pt);
    temp_data = PtData(LowFlow & current_pt,:);
    RemovedBB_Low(pt, 1) = pts(pt);
    RemovedBB_Low(pt, 2) = height(temp_data);
    RemovedBB_Low(pt, 3) = sum(temp_data.BB_Ttot);
    RemovedBB_Low(pt, 4) = nnz(temp_data.Etype==2);% 2 obstructive ap
    RemovedBB_Low(pt, 5) = nnz(temp_data.Etype==3);% 3 central ap
    RemovedBB_Low(pt, 6) = nnz(temp_data.Etype==4);% 4 hypopnoea ob
    RemovedBB_Low(pt, 7) = nnz(temp_data.Etype==5);% 5 mixed
end
RemovedBB_LowFlow = array2table(RemovedBB_Low, 'VariableNames', {'Pt','NumBB','TotalTime','ApO','ApC','HypO','Mixed'})
str = ['Removing ', num2str(nnz(LowFlow)), ' low flow breaths (<10% of eupnea)']; disp(str);
PtData(LowFlow,:)=[];

%% remove breaths that are NaN for ALL features
BBwithALLNans = find(all(isnan(PtData{:,firstdatacol:end}),2));
str = ['Removing ', num2str(length(BBwithALLNans)), ' breaths with ALL NaN features']; disp(str);
PtData(BBwithALLNans,:)=[];

BBwithManyNans = find(sum(isnan(PtData{:,firstdatacol:end}),2)>5);
str = ['Removing ', num2str(length(BBwithManyNans)), ' breaths with many NaN features']; disp(str);
PtData(BBwithManyNans,:)=[];

%% remove features with high NaN counts (using ~isfinite instead of isnan)
FeatureNames = PtData.Properties.VariableNames';  
SetThreshold = 0.001;% 0.001; % those with more than a set threshold of NaN
Threshold = round(size(PtData,1)*SetThreshold);
%NaNSpace = isnan(PtData{:,:});
NaNSpace = ~isfinite(PtData{:,:});
ftrsToExclude = find(sum(NaNSpace(:,:))>Threshold);
ftrsToExclude = ftrsToExclude(ftrsToExclude>=firstdatacol);
str = ['Removing ', num2str(length(ftrsToExclude)), ' features with more than ', num2str(SetThreshold*100), '% of NaN data']; disp(str);
FeatureNames(ftrsToExclude,:)=[];
PtData(:,ftrsToExclude)=[];

%   %    ~ Ftrs removed
%   3       4
%   1       6
%   0.5     8
%   0.1     10

%% Remove features in "FtrsToExclude" list, using string search - startswith
FtrsToExclude = {...
    'Ali_NED',...       % this is unreliably NaN
    'AA_EFLI',...       % this is unreliably NaN
    'MinFlowChange'...  % this is always zero
    };
Ind = [];
for i=1:length(FtrsToExclude) %needs checking
    temp=find(startsWith(FeatureNames,FtrsToExclude(i)));
    if ~isempty(temp)
        Ind = [Ind;temp];
    end
end
str = ['Removing ', num2str(length(Ind)), ' features, manually excluded (first pass)']; disp(str);
FeatureNames(Ind,:)=[];
PtData(:,Ind)=[];

%% wipe out the double duration test data (if it exists) - endswith
FtrsToExclude = {...
    '_O_250'...
    };
Ind = [];
for i=1:length(FtrsToExclude) %needs checking
    temp=find(endsWith(FeatureNames,FtrsToExclude(i)));
    if ~isempty(temp)
        Ind = [Ind;temp];
    end
end
str = ['Removing ', num2str(length(Ind)), ' features, manually excluded (second pass)']; disp(str);
FeatureNames(Ind,:)=[];
PtData(:,Ind)=[];

%% Remove features in "FtrsToExclude" list, using string search - contains
% timing variant not appropriate
FtrsToExclude = {...
%     'Flut'...
    'AsymIndex_T'... % don't keep asym index in Ttran timing
    'TTran_i_Ti_O'... % don't keep Ttran measures in Orig timing
    'TTran_i_Ttot_O'... % don't keep Ttran measures in Orig timing
    'TTran_e_Te_O'... % don't keep Ttran measures in Orig timing
    'TTran_e_Ttot_O'... % don't keep Ttran measures in Orig timing
    };
Ind = [];
for i=1:length(FtrsToExclude) %needs checking
    temp=find(contains(FeatureNames,FtrsToExclude(i)));
    if ~isempty(temp)
        Ind = [Ind;temp];
    end
end
str = ['Removing ', num2str(length(Ind)), ' features, manually excluded (third pass)']; disp(str);
FeatureNames(Ind,:)=[];
PtData(:,Ind)=[];


%% remove any remaining breaths that are NaN for ANY features
BBwithNans = find(any(isnan(PtData{:,firstdatacol:end}),2));
str = ['Removing ', num2str(length(BBwithNans)), ' breaths with ANY NaN features']; disp(str);
PtData(BBwithNans,:)=[];

%% remove breaths with no drive data
if FlowDriveData
    if EdiDrive
        %nnz(isnan(PtData{:,13}))
        BBwithNanDrive = find(any(~isfinite(PtData{:,[11 13]}),2));
        PesNotFinite = find(any(~isfinite(PtData{:,[12 14]}),2));
        PtData{PesNotFinite,[12 14]}=0; % set ~isfinite PesDrive columns to zero
        driveStr = '_Edi';
    else % PesDrive 12 14
        BBwithNanDrive = find(any(~isfinite(PtData{:,[12 14]}),2));
        EdiNotFinite = find(any(~isfinite(PtData{:,[11 13]}),2));
        PtData{EdiNotFinite,[11 13]}=0; % set ~isfinite PesDrive columns to zero
        driveStr = '_Pes';
    end
    str = ['Removing ', num2str(length(BBwithNanDrive)), ' breaths with NaN drive']; disp(str);
    PtData(BBwithNanDrive,:)=[];
else
    % remove breaths that have no VE./Veup value, as this is req'd in OA
    VEVeup = PtData.VE./PtData.Veup;
    BBwithNanDrive = find(~isfinite(VEVeup));
    str = ['Removing ', num2str(length(BBwithNanDrive)), ' breaths with NaN VE./Veup']; disp(str);
    PtData(BBwithNanDrive,:)=[];
    PtData{:,[11 12 13 14]}=0; % set Drive columns to zero for the OA data
    driveStr = '_VEVeup';
end

%% plot NaN's in data
if ShowFigures
    %NaNSpace = isnan(PtData{:,:});
    NaNSpace = (~isfinite(PtData{:,:})*90);
    figure(1); clf(figure(1)); image(NaNSpace);
    ylabel('BB''s');
    xlabel('Breath info (1-17), Features (18-end)');
    title('Nonfinite data shown as yellow area');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    fig.Position = [1 1 12 6];
    ax = gca;
    ax.TickDir = 'out';
    box('off');
end

%% Outliers
% Find extreme outliers.
% for e.g. VTi_VTe should be around 1, so values of 100 are clearly abnormal.
firstdatacol = 18;  % confirm this during procesing
%SS version labels outliers e.g. if 2 fold above the 95-50 difference, e.g. median = 5, 95th centile = 10; outliers >15.
prctilethres=99; %e.g. 95
prctilefactor=3; %e.g. 2
centileupper = prctile(PtData{:,firstdatacol:end},50) + (prctile(PtData{:,firstdatacol:end},prctilethres)-prctile(PtData{:,firstdatacol:end},50))*prctilefactor;
centilelower = prctile(PtData{:,firstdatacol:end},50) - (prctile(PtData{:,firstdatacol:end},50)-prctile(PtData{:,firstdatacol:end},100-prctilethres))*prctilefactor;
outliers = (PtData{:,firstdatacol:end} > centileupper | ...
    PtData{:,firstdatacol:end} < centilelower);
upper_outliers = nnz(PtData{:,firstdatacol:end} > centileupper);
lower_outliers = nnz(PtData{:,firstdatacol:end} < centilelower);
filler = false(size(PtData,1),firstdatacol-1); % exclude the first n columns of PtData, which are not feature data
outliers = ([filler outliers]);
perBB = any(outliers,2);
% Modify extreme outliers.
%SS version replaces outliers with limit value, set as upper/lower boundary
for n=1:length(centileupper)
    PtData{:,n+firstdatacol-1}(PtData{:,n+firstdatacol-1} > centileupper(n)) = centileupper(n);
    PtData{:,n+firstdatacol-1}(PtData{:,n+firstdatacol-1} < centilelower(n)) = centilelower(n);
end
str = [ num2str(upper_outliers),' high outliers found, ', num2str(lower_outliers), ' low outliers found']; disp(str);
str = [ num2str(nnz(outliers)),' feature outliers were clipped in ', num2str(nnz(perBB)), ' breaths']; disp(str);

%% One last check, find features that contribute equivalent information
% make new table as the zscore for each column
% for c=1:length(columns)-1, compare column c with c+1:end
sameinfo=[];
str=['Searching for redundant features by comparing zscores.']; disp(str);
str=['This may take a couple of minutes...']; disp(str);
[PtData_z, ~, ~] = zscore(PtData{:,firstdatacol:end});
PtData_z = round(PtData_z,2);
Name = PtData.Properties.VariableNames(firstdatacol:end)';
progressbar(0,0);
for c1 = 1:(size(PtData_z,2))-1
   for c2 = c1+1:(size(PtData_z,2)) 
       progressbar(c1/(size(PtData_z,2)), c2/(size(PtData_z,2)));
       %if isequal(round(PtData_z(:,c1),2), round(PtData_z(:,c2),2)) % then compare c1 and c2
       % tried to speed up with sum == 0, this "could" sum to zero with right combination of non-zeros, but is unlikely
       if sum(PtData_z(:,c1) - PtData_z(:,c2)) == 0 
           str=['Warning - ', char(Name(c1)), ' and ', char(Name(c2)),...
               ' contain potentially equivalent information']; disp(str);
           sameinfo = [sameinfo; [c1 c2]];
       end
   end
end
progressbar(1,1)

if ~isempty(sameinfo) % no error checking on empty array etc
    samesame = sameinfo+firstdatacol-1; % remember to add the cols back
    if RemoveRedundantFtrs
        for RedunFtrs = 1:(size(sameinfo,1))
            str=['Removing ', char(Name(sameinfo(RedunFtrs,2)))]; disp(str);
        end
        PtData(:,(samesame(:,2)))=[];
    end
else
    str=['No obvious redundant features found.']; disp(str);
end

%% Make the final version of the required outputs
%   Amatrix - matrix of Features (per col) and breaths (per row)
Amatrix = PtData{:,firstdatacol:end};
%   FeatureNames - table of column titles (i.e. Features) in Amatrix
Name = PtData.Properties.VariableNames(firstdatacol:end)';
FeatureNames = table((1:1:size(Amatrix,2))', Name, 'VariableNames', {'Number', 'Name'});
%   PtData - Table of patient and breath data
PtData = PtData(:,1:firstdatacol-1);

%unique(PtData.PT)

%% Save the required outputs
if SaveCleanMat
savefilename = [filename(1:end-4), driveStr,  '_Clean_new.mat'];
str=['Saving ' savefilename]; disp(str);
try
save(savefilename, 'Amatrix', 'FeatureNames', 'PtData', 'RemovedBB_Apnoea', 'RemovedBB_LowFlow'); 
catch me
    disp(me.message);
end
end


%% setup the binary response variable for classification
% set up grouping and criteria data.
% could remove pts who are not contributing. See BBperPt_table for guidance
if 0 % currently not doing any of this
DriveThreshold_upper = 0.7;
DriveThreshold_lower = 0.7;
pt_exclude = [];
PtData(pt_exclude,:)=[]; %truncate the PtData
n_pts = size(BBperPt,1);
if FlowDriveData
        Gtest_continuous = PtData.g_Edi_Adj;
        
        Gtest = ones(size(PtData,1),1)*2; % set a vector of two
        Gtest(Gtest_continuous>=DriveThreshold_upper) = 1; % NFL
        Gtest(Gtest_continuous<DriveThreshold_lower) = 0; % FL
        
        exclude = Gtest==2; nnz(exclude)
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
else % oral appliance data
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
            nnz(VEVEup(Gtest==2)>0.9&Ar1(Gtest==2))
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
        
        
end
end