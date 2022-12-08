function [BrDataTable,DriveExponent,Drive] = VdriveFromEdi(BrDataTable,VESignal,DriveSignal,EdiLinearization)


%% processing
%Set up signal names

Drive = abs(BrDataTable.(DriveSignal));
if sum(~isnan(Drive))<10
    DriveExponent=1;
    disp('No Drive Data');
    return
end
    
VE = BrDataTable.(VESignal);

if strcmp(DriveSignal,'DeltaEdi')
    DriveOutStr = 'VdriveEdi';
elseif strcmp(DriveSignal,'DeltaPes')
    DriveOutStr = 'VdrivePes';
elseif strcmp(DriveSignal,'DeltaPmus')
    DriveOutStr = 'VdrivePmus';
elseif strcmp(DriveSignal,'DeltaPepi')
    DriveOutStr = 'VdrivePepi';
else
    DriveOutStr = ['Vdrive' DriveSignal];
end

%% Identify breaths clearly away from sleep, currently relies entirely on LH arousal scoring here
BrDataTable=howfarawayfromsleep2(BrDataTable); %replace a with (a==1)&(hyp==4) if arousals not scored in wake...
ad = BrDataTable.nawayfromsleep;
%check enough breaths, change thresholds if not
if 1
    minNwakebreaths = 50;
    hh=hist(ad,[1:11]);
    hh(end)=NaN;
    th=find(hh>minNwakebreaths,1,'last');
    threshold = min([4 th]);
else
    %find breaths during wakefulness and arousals
    threshold = 2;
end
a1 = ad>threshold; % threshold # breaths away from sleep, usually 4 can vary, see above
a2 = ad>=2; % at least two breaths away from sleep
if 1 %added 20200528
    a2(VE==0)=0; %apneas are either not wake or erroneous data that are not useful.
end
ClearWakeorArousal = a2;

%% Transform/Linearize if requested
if EdiLinearization
    X=VE;
    Y=Drive; 
    x=X(ClearWakeorArousal==1);
    y=Y(ClearWakeorArousal==1);
    Brange=[1 4];
    [fitresult, gof] = fitAxPowerB(x, y,Brange);
    %Replace
    temp = Drive.^(1./fitresult.b);
    Drive = temp.*(nanmedian(Drive)./nanmedian(temp)); %transformed and untransformed Drive have same median value
    DriveExponent = fitresult.b;
else
    DriveExponent = 1;
end

%% Find what is normal mechanics (G = flow:drive, L/s/DriveUnits) over course of the night
%Start to identify breaths and times for estimating normal respiratory mechanics 
Tmta = BrDataTable.Time_start(a2); %Tmta is the time signal
    [Tmta,ia,ic] = unique(Tmta);
data_edi = VE(a2)./Drive(a2);
 data_edi = data_edi(ia);
 if 1
figure(201);clf
subplot(3,1,1)
histogram(VE(a2),20)
subplot(3,1,2)
histogram(Drive(a2),20)
subplot(3,1,3)
histogram(data_edi,20)
end
data_edi(data_edi>prctile(data_edi,75))=NaN; %%changed from 99,1 LKG test 0106
data_edi(data_edi<prctile(data_edi,25))=NaN;
% data_edi(data_edi>prctile(data_edi,99))=NaN; %%changed from 99,1 LKG test 0106
% data_edi(data_edi<prctile(data_edi,1))=NaN;
Tmta(isnan(data_edi))=[]; data_edi(isnan(data_edi))=[];
maxt = 1800;
Gmta_edi=NaN*Tmta;
for j=1:length(Tmta)
    temp2 = abs(Tmta(j)-Tmta);
    weights = maxt-temp2; weights(weights<0)=0; weights = weights/nansum(weights);
    Gmta_edi(j)=nansum(data_edi(weights>0).*weights(weights>0));
end
G_edi_w=sortrows([Tmta,Gmta_edi]);
G_edi_w_indiv=sortrows([Tmta,data_edi]);
str = [num2str(length(G_edi_w)),' ',DriveSignal,' reference breaths']; disp(str);
a2_locs = find(a2==1);
refBBintoftr = a2_locs(ia); % set the index of ref breaths, mainly for plotting

%Alternative, simple median:
G_edi_median = nanmedian(data_edi);

%Add in not mta version if preferred
G_edi = G_edi_w;

    
%% all breaths are then normalised, i.e. divided by their corresponding VE/Vdrive ratio
% Edi first
g_Edi_Adj = NaN(length(VE),1);

    G=interp1(G_edi(:,1),G_edi(:,2),BrDataTable.Time_start); %G_edi(:,1) is time, G_edi(:,2) is flow:drive in L/s/DriveUnits
    G(BrDataTable.Time_start<G_edi(1,1))=G_edi(1,2);
    G(BrDataTable.Time_start>G_edi(end,1))=G_edi(end,2);
    
    if 0
        figure(1); clf(1);
        plot(Tmta,60*data_edi,'.'); hold('on')
        plot(BrDataTable.Time_start,60*G,'r-');
        ylim([0 60*1.5*max(G)]);
    end
    
    FlowDrive = VE./(Drive.*G); %Flow drive in L/s per L/s (unitless)

    %%
            BrDataTable.(DriveOutStr) = Drive.*G; %L/s
      %      BrDataTable.(DriveOutStr) = Drive.*G_edi_median; %L/s
            BrDataTable.([DriveSignal 'mtaCal']) = Drive.*G/G_edi_median; %L/s
            
            if 1 %remove outliers
                tempcol = BrDataTable.(DriveOutStr);
                temp = prctile(tempcol,[50 90]);
                upper = diff(temp)*3+temp(1);
                tempcol(tempcol>2*upper)=NaN;
                tempcol(tempcol>upper)=upper;
                BrDataTable.(DriveOutStr) = tempcol;
            end
            
            %Normalize: Convert drive from L/s to eupnea
            BrDataTable.([DriveOutStr 'Norm']) = BrDataTable.(DriveOutStr)./BrDataTable.Veup; %Fraction of moving average eupnea value
  
        