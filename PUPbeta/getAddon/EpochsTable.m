function[Evts]=EpochsTable(SigT,Evts,AMasterSpreadsheet,settings,ii)

% this function makes Epochs level Table  -1/19/2021 RMA

%% setting up
dt = SigT{2,1} - SigT{1,1};
Fs = 1/dt;


    PosDatabaseDir = [settings.codedir filesep 'Position' filesep];
    PosDatabaseSpreadsheet = [PosDatabaseDir, 'PositionDatabase.xlsx']; %
    [~,~,settings.poscodesdatabase] = xlsread(PosDatabaseSpreadsheet,1,'B2:K55');

settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');

if (isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1)
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(ii,:));
else
    [~,~,settings.protocol] = xlsread(AMasterSpreadsheet,1,'AF4:AF10003');
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{ii});
end
%% make epochs table
Epochs = Evts.Hypnogram;
EpochsT = table(Epochs);
EpochsT.EpochStart = Evts.Hypnogram_t;

zerocheck = std(diff(EpochsT.EpochStart));

EpochsT.EpochNum=(1:size(EpochsT,1))';
EpochsT=[EpochsT(:,3) EpochsT(:,1:2)];

%% make events table | - only events meeting 3pa
try % DV added try catch for DISE data 8/4/21
RespT3pa = Evts.RespT;

ColsNotIncl = [2:20];
Ind0=RespT3pa.InclAHI3a==0;
temp = RespT3pa{:,ColsNotIncl};
temp(Ind0,:)=NaN;
RespT3pa{:,ColsNotIncl} = temp;


I = sum(RespT3pa.EventCodes==[2:6],2)>0;
RespT3pa = RespT3pa(I,:);


RespT3pa.EpochNumber = floor((RespT3pa.EventStart - EpochsT.EpochStart(1))/30 + 1);

temp = histc(RespT3pa.EpochNumber,(1:height(EpochsT))');
EpochsT.Nevents = temp(:); %SS edit, 6/22/21 to make upright and fix error
catch me
end
%% stages
% EpochsT.N1 = 1*(EpochsT.Epochs==2);
% EpochsT.N2 = 1*(EpochsT.Epochs==1);
% EpochsT.N3 = 1*(EpochsT.Epochs==0);
% EpochsT.REM = 1*(EpochsT.Epochs==3);
% EpochsWake=1*(EpochsT.Epochs==4);
% TST=size(EpochsT.Epochs,1)-sum(EpochsWake);   % TST as number of epochs
% sum(EpochsT.N1)+sum(EpochsT.N2)+sum(EpochsT.N3)+sum(EpochsT.REM)

%% Position : Flateral, Arousals

EpochsT.NArousals = NaN(size(EpochsT.Epochs,1),1);
ArNum= SigT.EventsAr;

% EpochsT.ODI3 = NaN(size(EpochsT.Epochs,1),1);
% EpochsT.ODI4 = NaN(size(EpochsT.Epochs,1),1);
% SpO2Ch = find(strcmp(ChannelsList,'SpO2'));
% % Remove SpO2 Artifacts
% SpO2Raw=SpO2ArtifactReject(DataEventHypnog_Mat(:,SpO2Ch),1/Fs);
% 

EpochsT.PosEpoch=NaN(size(EpochsT.Epochs,1),1);
EpochsT.Flateral = NaN(size(EpochsT.Epochs,1),1);
EpochsT.FSupine = NaN(size(EpochsT.Epochs,1),1);
try
    PosRaw = SigT.Position;
catch 
    PosRaw =NaN*SigT.Time;
end


for xx=1:size(EpochsT.Epochs,1)
    %code might be slow here, avoid using find 
    Idx1=find(SigT{:,1}>=EpochsT.EpochStart(xx),1,'first');
    Idx2=find(SigT{:,1}<=(EpochsT.EpochStart(xx)+30),1,'last')-1;
    
    temp = PosRaw(Idx1:Idx2);
    temp2 = PositionTranslator(positioncodes,settings.positioncodesout,temp);
    EpochsT.PosEpoch(xx) = mode(temp2);
    EpochsT.FSupine(xx) = nanmean(PositionSelector(temp2,'Supine'));
    EpochsT.Flateral(xx)=1-EpochsT.FSupine(xx); %technically lateral here means non-supine, and includes unknown and upright
    
    % Arousal
    clear temp I1
    temp=ArNum(Idx1:Idx2);
    I = diff([NaN;temp]);
    I1 = find(I==1);
    %     I2 = find(I==-1);
    %     [I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
    %     lengthsN = I2N-I1N;
    EpochsT.NArousals(xx)=size(I1,1);
    
    
%         SpO2temp=SpO2Raw(Idx1:Idx2);
%         Include=1-(SpO2temp==0); % include only when sao2 is non-zero
%         [ODI3temp,ODI4temp]=CalcODI(SpO2temp,dt,Include);
%         EpochsT.ODI3(xx) = ODI3temp;
%         EpochsT.ODI4(xx)=ODI4temp;
    
end

%% COMB 30 -- not used for now.

if 0
EpochsT2=EpochsT;

combW = 30;
I = 1:size(EpochsT2,1);
for j=1:3
    swap=j-1;
    clear I2 EpochT2Iteration
    
    
    %     combW = 30;
    %     epochnumber = [1:400]';
    Etime = (EpochsT.EpochStart - EpochsT.EpochStart(1))/60;
    %     swap = 1;
    temp2 = mod(Etime + combW*swap,2*combW)/(2*combW);
    criteriacomb = temp2<0.5;
    if j==3
        criteriacomb = logical(ones(size(EpochsT2,1),1));
    end
    %
    %         I2 = datasample(I,length(I))';
    EpochT2Iteration=EpochsT2(criteriacomb,:);
    
    EpochT2Iteration(EpochT2Iteration.Epochs>3|isnan(EpochT2Iteration.Epochs),:)=[];
    %stages
    TST(j)=size(EpochT2Iteration.Epochs,1); %TST is in number of epochs
    FN1temp(j,1)=(sum(EpochT2Iteration.N1==1))/TST(j);
    FN2temp(j,1)=(sum(EpochT2Iteration.N2==1))/TST(j);
    FN3temp(j,1)=(sum(EpochT2Iteration.N3==1))/TST(j);
    FRemtemp(j,1)=(sum(EpochT2Iteration.REM==1))/TST(j);
    %position
    FLateraltemp(j,1)=nanmean(EpochT2Iteration.Flateral==1);
    %AHI
    TotHr=TST(j)*30/3600;
    AHItemp(j,1)=(sum(EpochT2Iteration.Nevents))/TotHr;
    %ArI
    ArItemp(j,1)=(sum(EpochT2Iteration.NArousals))/TotHr;
    %TST Total
    TSTtemp(j,1)=TST(j)*30/60; % in minutes
    %ODI
    Timetemp(j,1)=nanmean(EpochT2Iteration.EpochStart);
end

%%
subj=ii*(ones(3,1));
J = [1 2 3]';
Type = {'Odd','Even','Total'}';
CombTtemp=table(subj,J,Type,FN1temp,FN2temp,FN3temp,FRemtemp,FLateraltemp,AHItemp,ArItemp,TSTtemp,Timetemp);
CombTtemp.Properties.VariableNames={'Subject','TypeID','Type','FN1','FN2','FN3','FREM','FLateral','AHI',...
    'ArI','TST','MeanTime'};

if ~isempty(CombTtemp)
    if ~exist('CombTPSGParam')
        CombTPSGParam = CombTtemp;
    else
        CombTPSGParam = [CombTPSGParam;CombTtemp];
    end
end
clear CombTtemp FN1temp FN2temp FN3temp FRemtemp FLateraltemp AHItemp ArItemp;
% end subject analysis


end

Evts.EpochT=EpochsT;
