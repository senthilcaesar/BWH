%% Figure 5
% ExpFlutPowOrig
% InspFlutPowOrig
% QuadE (SinE)
% QuadI50 (SinI50)
% SS_Area


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Load PtData <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
close; clear; clc
addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO'); % code
datadir = 'C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\'; % data
experimentnumber = '_n08';
load([datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TrainFlow_TestFlowAndPnasal', experimentnumber, '.mat'], 'PtData_flow');

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Settings  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subj = 32; % set patient/subject
savefigas = '';%'saveasPNG'; %   ''; % 'saveasTIFF'; % options are saveasPNG, saveasFIG, and saveasTIFF
closefigs = 0;
addNaNgaps=1;
useaxes = 1;
TickFntSz = 12;
LabelFntSz = 18;
FntSz = 18;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Load Signals <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:G56');
studyname = char(patients{subj});
studynameNoExt = char(patients{subj}(1:end-4));
load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\',studyname], 'Edi', 'Flow', 'StarttimeSpike');
try
    load(['C:\PSG_Data\FlowDrive\Converted\',studynameNoExt,'_XHz.mat']);
catch GetEvtsFail
    disp(GetEvtsFail.getReport);
end
EventsAr = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1));
EventsResp = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsResp')==1));
% arousal   = 1
% ap-O      = 2
% ap-C      = 3
% hyp-O     = 4
% mixed     = 5
% hyp-C     = 6

if ~(exist('StarttimeSpike', 'var') == 1)
    StarttimeSpike = 0;
end
Time = StarttimeSpike:0.008:StarttimeSpike+(length(Flow.values)-1)*0.008;
Isubj = PtData_flow.PT==subj;
Data1 = [PtData_flow.BB_time(Isubj) PtData_flow.BB_Ttot(Isubj)];
cols = [NaN 1 2];
if addNaNgaps
    tol2=0.1; i=1; M=size(Data1,2);
    while i<(size(Data1,1)-1)
        if (Data1(i,cols(2))+Data1(i,cols(3))+tol2)<Data1(i+1,cols(2))
            Data1 = [Data1(1:i,:); NaN*ones(1,M); Data1((i+1):size(Data1,1),:)];
            i=i+1;
        end
        i=i+1;
    end
end

if Flow.length<=length(Time); Flow.values(end:length(Time))=0; else; Flow.values(length(Time):end)=[]; end
if Edi.length<=length(Time); Edi.values(end:length(Time))=0; else ;Edi.values(length(Time):end)=[]; end

dsf=5; dt=Flow.interval; FlowF=Flow.values;
EdiF = Edi.values;
if 1
    filter_HFcutoff_butter0 = 12.5;
    filter_order0 = 1;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
    EdiF = filtfilt(B_butter0,A_butter0,EdiF); %filtfilt, otherwise flow signal is right-shifted
end

% use the pt timing in PtData to locate breaths in pt signals
Isubj = PtData_flow.PT==subj;
BB_time_pt = PtData_flow.BB_time(Isubj);
BB_Ttot_pt = PtData_flow.BB_Ttot(Isubj);

%% all breaths
figure(500); clf(figure(500)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  1   15   7];
ax(1) = subplot(2,1,1);
plot(Time, FlowF); hold on;
ax(2) = subplot(2,1,2);
plot(Time, EdiF); hold on;
plot(BB_time_pt, -5, 'b+', 'Markersize', 2);
linkaxes(ax, 'x');

%% breath number
%n=find(BB_time_pt>9.07e4,1,'first');
%n=483;
%for n=1:100
for n = 4
n2 = 31;

% Figure 5 - Top features
% set up the data for breath 1
CB_Start = BB_time_pt(n);
CB_End = CB_Start + BB_Ttot_pt(n);
CB_Start_i = find(Time>CB_Start, 1,'first');
CB_End_i = find(Time>CB_End, 1,'first');
BT = Time(CB_Start_i:CB_End_i);
BF = (FlowF(CB_Start_i:CB_End_i))/max(FlowF(CB_Start_i:CB_End_i));
BV = cumsum(BF);
[~,BB_mid] = max(BV);
BV = BV/max(BV);
InspTime = BT(1:BB_mid)'; ExpTime = BT(BB_mid:end)';
InspFlow = BF(1:BB_mid)'; ExpFlow = BF(BB_mid:end)';
dt=(Time(end)-Time(1))/(length(Time)-1);
Fs = round(1/dt);
[Pxx1,f1] = periodogram(InspFlow,[],[],Fs);
[Pxx2,f2] = periodogram(ExpFlow,[],[],Fs);
Pow1_Orig=bandpower(Pxx1,f1,[5,round(Fs/2)-1],'psd');
Pow2_Orig=bandpower(Pxx2,f2,[5,round(Fs/2)-1],'psd');
[fft_x_insp_b1, fft_y_insp_b1] = CheckSignalForNoisewFig(BT(1:BB_mid),InspFlow, Fs);
[fft_x_exp_b1, fft_y_exp_b1] = CheckSignalForNoisewFig(BT(BB_mid:end),ExpFlow, Fs);

% set up the data for breath 2
CB_Start_b2 = BB_time_pt(n2);
CB_End_b2 = CB_Start_b2 + BB_Ttot_pt(n2);
CB_Start_i_b2 = find(Time>CB_Start_b2, 1,'first');
CB_End_i_b2 = find(Time>CB_End_b2, 1,'first');
BT_b2 = Time(CB_Start_i_b2:CB_End_i_b2);
BF_b2 = (FlowF(CB_Start_i_b2:CB_End_i_b2))/max(FlowF(CB_Start_i_b2:CB_End_i_b2));
BV_b2 = cumsum(BF_b2);
[~,BB_mid_b2] = max(BV_b2);
BV_b2 = BV_b2/max(BV_b2);
InspTime_b2 = BT_b2(1:BB_mid_b2)'; ExpTime_b2 = BT_b2(BB_mid_b2:end)';
InspFlow_b2 = BF_b2(1:BB_mid_b2)'; ExpFlow_b2 = BF_b2(BB_mid_b2:end)';
dt=(Time(end)-Time(1))/(length(Time)-1);
Fs = round(1/dt);
[Pxx1_b2,f1_b2] = periodogram(InspFlow_b2,[],[],Fs);
[Pxx2_b2,f2_b2] = periodogram(ExpFlow_b2,[],[],Fs);
Pow1_Orig_b2=bandpower(Pxx1_b2,f1_b2,[5,round(Fs/2)-1],'psd');
Pow2_Orig_b2=bandpower(Pxx2_b2,f2_b2,[5,round(Fs/2)-1],'psd');
[fft_x_insp_b2, fft_y_insp_b2] = CheckSignalForNoisewFig(BT_b2(1:BB_mid),InspFlow_b2, Fs);
[fft_x_exp_b2, fft_y_exp_b2] = CheckSignalForNoisewFig(BT_b2(BB_mid:end),ExpFlow_b2, Fs);


%% smooth a little bit for demo figure - it's a bit busy
MA = 4;
fft_y_insp_b1 = smooth(fft_y_insp_b1, MA )';
fft_y_exp_b1 = smooth(fft_y_exp_b1, MA )';
fft_y_insp_b2 = smooth(fft_y_insp_b2, MA )';
fft_y_exp_b2 = smooth(fft_y_exp_b2, MA )';


%% do the plot (more data calculated as plot is made)
figure(500+subj); clf(figure(500+subj));fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [22  1   10   8];

% -----------------------------------------------------------------------
% FFT plot - insp 
if useaxes
    % inset above two panel
    %axes('Position',[0.356 0.8 0.135 0.18]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]  
    % to left of two panel
    %axes('Position',[0.05 0.55 0.12 0.35]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]  
    % to righ of two panel
    axes('Position',[0.82 0.55 0.12 0.35]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]  
else
    subplot(2,5,1); 
end
ax = gca; 
ind = find(fft_x_insp_b1>5, 1, 'first'); fft_y_insp_b1(1:ind) = NaN;
ind = find(fft_x_insp_b1>12.5, 1, 'first'); fft_y_insp_b1(ind:end) = NaN;
fft_x_insp_b1(fft_x_insp_b1<0) = NaN; fft_y_insp_b1(fft_y_insp_b1<0) = NaN;

% shade area under curve for good breath blue
Lowercurve = zeros(1, length(fft_y_insp_b1));
idx = ~isnan(fft_y_insp_b1); 
jbfill(fft_x_insp_b1(idx), fft_y_insp_b1(idx), Lowercurve(idx), [0 0 1], [0 0 1], 0, 0.1); 
hold on
plot(fft_x_insp_b1(idx), fft_y_insp_b1(idx), 'b','linewidth', 1.5);

ind = find(fft_x_insp_b2>5, 1, 'first'); fft_y_insp_b2(1:ind) = NaN;
ind = find(fft_x_insp_b2>12.5, 1, 'first'); fft_y_insp_b2(ind:end) = NaN;
fft_x_insp_b2(fft_x_insp_b2<0) = NaN; fft_y_insp_b2(fft_y_insp_b2<0) = NaN;

% shade area under curve for bad breath red
Lowercurve = zeros(1, length(fft_y_insp_b2));
idx = ~isnan(fft_y_insp_b2); 
jbfill(fft_x_insp_b2(idx), fft_y_insp_b2(idx), Lowercurve(idx), [1 0 0], [1 0 0], 0, 0.1);
plot(fft_x_insp_b2(idx), fft_y_insp_b2(idx), 'r','linewidth', 1.5);

xlim([4 20]); yticks([]);
ylabel('power (inspiration)','FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Hz','FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); %xlabel({'frequency (Hz)','in expiration'});
set(ax, 'XScale', 'log'); % set(ax, 'YScale', 'log') % log y kills the shading
set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);

% -----------------------------------------------------------------------
% FFT plot - exp
if useaxes
    % inset above two panel
    % axes('Position',[0.769 0.8 0.135 0.18]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]
    % to left of two panel
    %axes('Position',[0.05 0.15 0.12 0.35]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]
    % to right of two panel
    axes('Position',[0.82 0.15 0.12 0.35]); % use axes instead. [lowerleftX, lowerleftY, boxX, boxY]
else
    subplot(2,5,6); 
end
ax = gca; 
ind = find(fft_x_exp_b1>5, 1, 'first'); fft_y_exp_b1(1:ind) = NaN;
ind = find(fft_x_exp_b1>12.5, 1, 'first'); fft_y_exp_b1(ind:end) = NaN;
fft_x_exp_b1(fft_x_exp_b1<0) = NaN; fft_y_exp_b1(fft_y_exp_b1<0) = NaN;

% shade area under curve for good breath blue
Lowercurve = zeros(1, length(fft_y_exp_b1));
idx = ~isnan(fft_y_exp_b1); 
jbfill(fft_x_exp_b1(idx), fft_y_exp_b1(idx), Lowercurve(idx), [0 0 1], [0 0 1], 0, 0.1); 
hold on
plot(fft_x_exp_b1(idx), fft_y_exp_b1(idx), 'b','linewidth', 1.5);

ind = find(fft_x_exp_b2>5, 1, 'first'); fft_y_exp_b2(1:ind) = NaN;
ind = find(fft_x_exp_b2>12.5, 1, 'first'); fft_y_exp_b2(ind:end) = NaN;
fft_x_insp_b2(fft_x_exp_b2<0) = NaN; fft_y_exp_b2(fft_y_exp_b2<0) = NaN;

% shade area under curve for bad breath red
Lowercurve = zeros(1, length(fft_y_exp_b2));
idx = ~isnan(fft_y_exp_b2); 
jbfill(fft_x_exp_b2(idx), fft_y_exp_b2(idx), Lowercurve(idx), [1 0 0], [1 0 0], 0, 0.1);
plot(fft_x_exp_b2(idx), fft_y_exp_b2(idx), 'r','linewidth', 1.5);

xlim([4 20]); yticks([]);
ylabel('power (expiration)','FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Hz','FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); %xlabel({'frequency (Hz)','in expiration'});
set(ax, 'XScale', 'log'); % set(ax, 'YScale', 'log') % log y kills the shading
set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);

% -----------------------------------------------------------------------
% good breath
subplot(2,5,[1:2,6:7]); ax = gca; 
%subplot(2,5,[2:3,7:8]); ax = gca; 

plot(BT, BF, 'b','linewidth', 1.5); hold on;
ylim([-1.5 1.5]);

addSSArea(InspFlow, InspTime, [0, 0, 1], [0, 0, 1]);
addCurveFits(InspTime, InspFlow, [0, 0, 1], ExpTime, ExpFlow,  [0, 0, 1]);

set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
%ylim([-1.5 1.1]);
ylabel('Flow (normalized)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
%xlabel('time (seconds)');
rline1 = refline(0,0); rline1.Color = [0.8 0.8 0.8];
set(gca,'xtick',[],'box','off');
set(gca,'xcolor',[1 1 1]);

titlestr=num2str(n);
title(titlestr);


% -----------------------------------------------------------------------
% bad breath
subplot(2,5,[3:4,8:9]); ax = gca; 
%subplot(2,5,[4:5,9:10]); ax = gca; 

plot(BT_b2, BF_b2, 'r','linewidth', 1.5); hold on;
ylim([-1.5 1.5]);

addSSArea(InspFlow_b2, InspTime_b2, [1, 0, 0], [1, 0, 0]);
addCurveFits(InspTime_b2, InspFlow_b2, [1, 0, 0], ExpTime_b2, ExpFlow_b2,  [1, 0, 0]);

set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
rline2 = refline(0,0); rline2.Color = [0.8 0.8 0.8];
set(gca,'xtick',[],'box','off');
set(gca,'xcolor',[1 1 1]);
set(gca,'ytick',[],'box','off');
set(gca,'ycolor',[1 1 1]);

titlestr=num2str(n2);
title(titlestr);

% need to add scale bar for time

str = ['..\Figures\Figure5\TopFtrs_Pt',num2str(subj),'_BB',num2str(n)];

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%pause(2);
end
