function [BreathDataTableE,BreathDataTableE_,BreathFLDataTableE]=GetNonOvlappedVE(varargin)

cut = 0.67;

%% Expand Tables
BreathDataTable=varargin{1};
BreathDataTableEtemp = expandtable(BreathDataTable);

% Clean flow shape table if variable is included as input
% ... (DLM) but, if it's not included in input, we can't do this fn anyway ?
if length(varargin) == 2
    BreathFLDataTable=varargin{2};
    BreathFLDataTableEtemp = expandtable(BreathFLDataTable);
    if ~isa(BreathFLDataTableEtemp, 'table') && all(size(BreathFLDataTableEtemp)==[1 1]) && isnan(BreathFLDataTableEtemp)
        % we get here if BreathFLDataTable is not a table, is size 1, and is NaN
        % if we are here, we can't do this analysis, so return
        BreathDataTableE = NaN;
        BreathDataTableE_ = NaN;
        BreathFLDataTableE = NaN;
        return
    end
    BreathDataTableFull = outerjoin(BreathDataTableEtemp,BreathFLDataTableEtemp,'Keys',{'UniqueID','UniqueID'},...
        'MergeKeys',1, 'Type', 'Left');
    BreathDataTableE = BreathDataTableFull(:,1:size(BreathDataTableEtemp,2));
    BreathFLDataTableE = BreathDataTableFull(:,...
        size(BreathDataTableEtemp,2)+1:size(BreathDataTableFull,2));
    BreathFLDataTableE.UniqueID = BreathDataTableE.UniqueID;
end

clear BreathDataTableEtemp BreathFLDataTableEtemp

%% Find dupilcate breaths
FDuplicated = FindDuplicated(BreathDataTable);
BreathDataTableE.FDuplicated = FDuplicated;

%% Remove Duplicates and Tidy Table (make this a function later)
BreathDataTableE_ = BreathDataTableE;
BreathDataTableE_(FDuplicated>cut,:) = [];

%% Find duplicate breaths part 2
FDuplicated2 = FindDuplicatedPt2(BreathDataTableE_);
BreathDataTableE_.FDuplicated2=FDuplicated2;

Temp = outerjoin(BreathDataTableE,BreathDataTableE_,'Keys',[1 2]);
FDuplicated2_ = Temp.FDuplicated2;

BreathDataTableE.FDuplicated2=BreathDataTableE.FDuplicated;
BreathDataTableE.FDuplicated2(~isnan(FDuplicated2_))=FDuplicated2_(~isnan(FDuplicated2_));

BreathDataTableE_(FDuplicated2==1,:) = [];

clear FDuplicated FDuplicated2 FDuplicated2_ Temp cut
%%

% figure(1); clf(1);
%
% I=find(diff(BreathDataTableE_.Time0)>0);
% gridheight=0.5;
% flowplotpos=0.5;
% xgridlist = BreathDataTableE_.Time_start(I+1)';
% plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^2);
% hold('on')
%
% plot(BreathDataTableE.Time_start,BreathDataTableE.VI,'r.-')
%
% stairs(BreathDataTableE_.Time_start,BreathDataTableE_.VI)
%


