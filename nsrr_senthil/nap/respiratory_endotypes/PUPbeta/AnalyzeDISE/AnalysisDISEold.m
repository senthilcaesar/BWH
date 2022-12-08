global AMasterSpreadsheet settings ChannelsList ChannelsFs

settings1 = ImportSettingsAnalysis(settings,AMasterSpreadsheet); 
settings=settings1;
[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));
settings.invertflowlist = logical(num(:,1));
settings.protocol = patients(:,3);
TotalNumPts=size(patients,1);
settings.patients = patients;

PtRangeTemp = 1:1:TotalNumPts; %normally
PtRange = PtRangeTemp(analyzelist==1);
%% loop
for ptnum=PtRange
    %% Load converted data
    settings.filename=[char(settings.patients(ptnum,1))]; %seems to be unused. 
    % DLM says, settings.filename is used in LGfromFlowBetaPart1 to
    % save flowdrive plot, and Part2 to save loop gain plot
    directoryn=char(settings.patients(ptnum,2));
    MATfilename=[directoryn char(settings.patients(ptnum,1))];
    Evts=struct();
    temp=[];
    temp=load(MATfilename);
    DataEventHypnog_Mat=temp.DataEventHypnog_Mat;
    Evts=temp.Evts;
    %WakeSleepInfo=temp.WakeSleepInfo; % removed after modifying 'Info' in Convert 
    ChannelsFs=temp.ChannelsFs;
    ChannelsList=temp.ChannelsList;
    
    % load snore channels list
    SnoreMATfilename = [MATfilename(1:end-4), '_snore'];
    temp2=load(SnoreMATfilename);
    SnoreChannelsList=temp2.SnoreChannelsList;

    %% Sample data
    Time = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Time'));
    Flow = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Flow'));
    Snore = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Snore'));
    SnoreDB = DataEventHypnog_Mat(:,strcmp(ChannelsList,'SnoreDB'));
    fsFlow = ChannelsFs(strcmp(ChannelsList,'Flow'));
    vol=cumsum(Flow(:))*(1/fsFlow);
    
    %% Manually select breath start stop
    % Find initial peaks and troughs
    % Call the function "get_vmin2" to find the troughs of the volume signal
    ti = 2; te=2;
    [vmin,ivmin] = get_vmin2(vol,ti,te,1.4,fsFlow);
    tvmin = ivmin/fsFlow;

    % Find the peaks of the volume signal
    vmax = [];
    ivmax = [];
    for n = 2:length(ivmin)
        [maxvol,indmax] = max(vol(ivmin(n-1):ivmin(n)));
        indmax = (ivmin(n-1)+indmax)-1;
        vmax = [vmax maxvol];
        ivmax = [ivmax indmax];
    end
    vmax = vmax';
    ivmax = ivmax';
    tvmax = ivmax/fsFlow;
    clear indmax maxvol n;
    % Output = BB_Times
    
    %% Adjust vmin
    [vmin, ivmin] = BBSelectManual(Time,vol,Flow,Snore,SnoreDB,ivmin,vmin,ivmax,fsFlow);
    
    %% Find Vmax again now that we have a better 
    tvmin = ivmin/fsFlow;

    clear ivmax tvmax vmax

    % re-find the peaks of the volume signal after you have corrected the
    % troughs
    vmax = [];
    ivmax = [];
    for n = 2:length(ivmin)
        [maxvol,indmax] = max(vol(ivmin(n-1):ivmin(n)));
        indmax = (ivmin(n-1)+indmax)-1;
        vmax = [vmax maxvol];
        ivmax = [ivmax indmax];
    end
    tvmax = ivmax/fsFlow;
    clear indmax maxvol n;

    % Correct any tvmax's that equal tvmin's
    for n = 1:length(ivmax)
        if ivmax(n) == ivmin(n)
            ivmax(n) = ivmin(n) + round((ivmin(n+1)-ivmin(n))/2); %put it halfway between the ivmins
        end
    end

    ivmax = ivmax';
    vmax = vmax';
    
    %% Adjust vmax
    [vmax, ivmax] = BBSelectManual(Time,vol,Flow,Snore,SnoreDB,ivmax,vmax,ivmin,fsFlow);

    %% One more time makes sure Vmin looks good
    while 1
        answer = questdlg('Review breath times again?', ...
        'Breath times review', ...
        'Yes, insp start times','Yes, insp stop times','No, all done','No, all done');
        % Handle response
        switch answer
            case 'Yes, insp start times' 
                [vmin, ivmin] = BBSelectManual(Time,vol,Flow,Snore,SnoreDB,ivmin,vmin,ivmax,fsFlow);

            case 'Yes, insp stop times'
                [vmax, ivmax] = BBSelectManual(Time,vol,Flow,Snore,SnoreDB,ivmax,vmax,ivmin,fsFlow);

            case 'No, all done' % now check that everything looks right, or else re-run loop
                %% Check times
                % vmin(1) has to be lowest number (has to start with inspration)
                if sum(ivmin(1) > ivmax(1)) > 0
                    ivmax(1) = [];
                end

                % vmin(end) has to be highest number (has to end with end-expiration,
                % i.e start inspiration)
                if sum(ivmin(end) > ivmax(end)) > 0
                    ivmax(end) = nan;
                end

                if length(ivmin) - length(ivmax) == 1
                    ivmax(end+1) = nan;
                elseif length(ivmin) ~= length(ivmax)
                    disp('Start and stop times vectors are different sizes - investigate this')
                    % make same size
                    if length(ivmin) < length(ivmax)
                        difflength = length(ivmax)-length(ivmin);
                        ivmin(end+1:end+difflength) = nan;
                    else
                        difflength = length(ivmin)-length(ivmax);
                        ivmax(end+1:end+difflength) = nan;
                    end
                end

                % check breath times        
                breathtimes = (ivmax - ivmin)./fsFlow;
                abberantBB = find(breathtimes < 0 | breathtimes > 4);
                nanBB = find(isnan(breathtimes(2:end-1)));
                if ~isempty(abberantBB) || ~isempty(nanBB)
                    disp('Breath(s) with unrealistic times located at:')
                    abberantBBtimes = ivmin(abberantBB)./fsFlow

                    disp('NaN breath(s) located at:')
                    nanBBtime = ivmin(nanBB+1)./fsFlow
                    
                    % find issue areas
                    interactive_vmin3(Time,vol,Flow,Snore,SnoreDB,ivmin,ivmax,1,length(Flow),1);
                    pause; disp('script paused to explore data for abberant breaths')
                else
                    break % while loop ends once no abberant breaths found
                end
        end
    end

    %% Compute ventilation - can use this to decide if flow shape is reliable
    BB_i_start = ivmin(1:end-1); % start inspiration
    BB_i_mid = ivmax(1:end-1);% end inspiration/start expiration
    BB_i_end = ivmin(2:end);
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)*(1/fsFlow);
    Te = (BB_i_end-BB_i_mid)*(1/fsFlow);
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    VE = VT./Ttot; %force these to be the same for now until LGfromFlow is updated
    
    %% Construct breath data table
    Apnea_B = zeros(1,length(BB_Times)); % can also base this on ventilation
    BB_Ttrans = []; TiTrans = nan(length(BB_i_start),1); 
    TeTrans = nan(length(BB_i_start),1); TAA_ = [];
    Time0 = zeros(length(BB_i_start),1);
    
    BB_Times(:,1) = BB_i_start; 
    BB_Times(:,2) = BB_i_mid; 
    BB_Times(:,3) = BB_i_end; 
    
    BreathDataTable=table(Time0, Time(BB_i_start), Time(BB_i_mid),...
        Time(BB_i_end), BB_i_start, BB_i_mid, BB_i_end, VE(:));
    
    [BreathFLDataTable] = ComputeBreathFeatures(Time, Flow, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[0 1 0]); %[0: downsampleHz; 1=original timing, no resample]
    
    [BreathSnoreTable] = ComputeSnoreFeatures(Time,DataEventHypnog_Mat,SnoreChannelsList,BB_Times);

    %% Save data
    datatosave = struct('ivmin',ivmin, 'vmin',vmin,'ivmax',ivmax,'vmax',vmax,...
        'BreathDataTable',BreathDataTable,'BreathFLDataTable',BreathFLDataTable,...
        'BreathSnoreTable',BreathSnoreTable);
    cd(settings.OutputDataDirectory);
    savename = [patients{ptnum,1}(1:end-4) '_A.mat'];
    save(savename,'-struct','datatosave')
end
