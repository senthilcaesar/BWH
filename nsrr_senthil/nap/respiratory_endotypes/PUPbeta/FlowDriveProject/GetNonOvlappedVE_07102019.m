function [BreathDataTableE,BreathDataTableE_]=GetNonOvlappedVE(BreathDataTable)

cut = 0.67;


%% Find dupilcate breaths
%Expand Tables
BreathDataTableE = expandtable(BreathDataTable);

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

      
