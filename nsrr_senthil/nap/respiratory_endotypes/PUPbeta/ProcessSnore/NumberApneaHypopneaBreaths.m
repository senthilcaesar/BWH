function BreathDataTable = NumberApneaHypopneaBreaths(BreathDataTable)

    if sum(strcmp(fieldnames(BreathDataTable), 'HypopNum'))
        return
    end

    %% Isolate Hypopnea Breaths
    %Generate large arrays with only hypopnea breath
    %Add new column indicating hypopnea number
    %Can later filter by hypopnea number
    Etype4 = BreathDataTable.Etype;
    Etype4(Etype4 ~= 4 & Etype4 ~= 6) = 0;
%     Etype4(Etype4 ~= 4) = 0;
    
    HypopShiftDown = circshift(Etype4,1);
    HypopShiftUp = circshift(Etype4,-1);
    HypopDiffDown = Etype4 - HypopShiftDown;
    HypopDiffUp = Etype4 - HypopShiftUp;

    HypopFirstIdx = find(HypopDiffDown == 4 | HypopDiffDown == 6);
    HypopLastIdx = find(HypopDiffUp == 4 | HypopDiffDown == 6);
    
%     HypopFirstIdx = find(HypopDiffDown == 4);
%     HypopLastIdx = find(HypopDiffUp == 4);

    %Generate variable indicating hypopnea number to be stored in table
    HypopNum = zeros(size(Etype4,1),1);
    NumHypop = zeros(size(Etype4,1),1);
    %Loop through hypopnea indices and fill hypopnea number variable
    %Note: 1 is the terminal hypopnea, 2 is second last, etc
    for ii = 1:length(HypopFirstIdx)
        totalHypop = HypopLastIdx(ii) - HypopFirstIdx(ii) + 1;
        HypopNum(HypopFirstIdx(ii):HypopLastIdx(ii),1) = totalHypop:-1:1;
        NumHypop(HypopFirstIdx(ii):HypopLastIdx(ii),1) = totalHypop;
    end

    % Add to table
    BreathDataTable.HypopNum = HypopNum;
    BreathDataTable.NumHypop = NumHypop;

    %% Isolate Apnea breaths
    Etype2 = BreathDataTable.Etype;
%     Etype2(Etype2 ~= 2 & Etype2 ~=3 & Etype2 ~= 5) = 0;
    Etype2(Etype2 ~= 2 & Etype2 ~= 5) = 0; %exclude central events
    
    ApneaShiftDown = circshift(Etype2,1);
    ApneaShiftUp = circshift(Etype2,-1);
    ApneaDiffDown = Etype2 - ApneaShiftDown;
    ApneaDiffUp = Etype2 - ApneaShiftUp;
% 
%     ApneaFirstIdx = find(ApneaDiffDown == 2 | ApneaDiffDown == 3 | ApneaDiffDown == 5);
%     ApneaLastIdx = find(ApneaDiffUp == 2 | ApneaDiffUp == 3 | ApneaDiffUp == 5);
    
    ApneaFirstIdx = find(ApneaDiffDown == 2 | ApneaDiffDown == 5); % exclude central
    ApneaLastIdx = find(ApneaDiffUp == 2 | ApneaDiffUp == 5);

    %Generate variable indicating apnea number to be stored in table
    ApneaNum = zeros(size(Etype2,1),1);
    NumApnea = zeros(size(Etype2,1),1);
    %Loop through apnea indices and fill apnea number variable
    %Note: 1 is the terminal hypopnea, 2 is second last, etc
    for ii = 1:length(ApneaFirstIdx)
        totalApnea = ApneaLastIdx(ii) - ApneaFirstIdx(ii) + 1;
        ApneaNum(ApneaFirstIdx(ii):ApneaLastIdx(ii),1) = totalApnea:-1:1;
        NumApnea(ApneaFirstIdx(ii):ApneaLastIdx(ii),1) = totalApnea;
    end

    % Add to table
    BreathDataTable.ApneaNum = ApneaNum;
    BreathDataTable.NumApnea = NumApnea;
    