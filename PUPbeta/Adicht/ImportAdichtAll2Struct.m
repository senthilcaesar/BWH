
%addpath(genpath([pathpiece 'Luos Data\ADIMAT'])); %read ADI toolset
%filen = [pathpiece 'Luos Data\Source\' filename '\' filename '.adicht'];

filen = [directory fnameadicht];

f = adi.readFile(filen);
channelNames = f.channel_names;

%%

clear LCchannelnameoptions
        LCchannelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal','Nasal Pet CO¤2'};
        LCchannelnameoptions.PO2={'pO2','O2_Ana','O2_anal','Pet O2'};
        LCchannelnameoptions.Pepi={'S. Glot Pr'};
        LCchannelnameoptions.Pmask={'Mask Pr'};
        LCchannelnameoptions.SpO2ear={'Ear SaO¤2'};
        LCchannelnameoptions.SpO2finger2={'Finger SaO2'};
        LCchannelnameoptions.FlowSync={'Flow'}; %
        LCchannelnameoptions.TcPCO2={'TCO2'}; %
        
        channelnamestemp=fieldnames(LCchannelnameoptions);
        
        channelfound=nan(1,length(channelnamestemp));
        for i=1:length(channelnamestemp)
            temp=eval(['LCchannelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            Amatch = string(temp)==string(channelNames(:));
            temp3 = find(any(Amatch,1)==1);
            temp2 = find(any(Amatch,2)==1);
            if ~isempty(temp2)
                disp(['found channel' channelnamestemp{i} ' in ' channelNames{temp2}])
                channelfound(i) = temp2;
            end
        end

I = isnan(channelfound) ;
channelnamestemp(I)=[];
channelfound(I)=[];

%%
channelSpecs = f.channel_specs;

for i=1:length(channelNames)
    FsList(1,i) = mode(channelSpecs(i).fs);
    FsSDList(1,i) = nanstd(channelSpecs(i).fs); %check zero, else problem
end

% channelNamesNew=channelNames;
% for i=1:length(channelNamesNew)
%     channelNamesNew{i}(channelNamesNew{i}==' ')=[];
% end

T=[];
for i=1:length(channelnamestemp)
    temp = importLabChartData(f,channelNames{channelfound(i)},f.n_records);
    T = setfield(T,channelnamestemp{i},temp);
end

%find start time of LC file (handle gaps; largest segment has exact time)
nRecords = length(f.records);
    TimeData=[];
    clear temp
    for ii=1:nRecords
        temp(ii,1) = f.records(ii).record_start;
        temp3(ii,1) = f.records(ii).n_ticks;
    end
    temp2 = temp - floor(temp(1));
    if temp2(1)<0.5
        temp2 = temp2+1;
    end
    temp2 = temp2*86400;
    [~,idx] = max(temp3);
    Fs = f.records(idx).tick_fs;
    temp5 = cumsum(temp3);
    temp4 = [temp2 temp3 temp5];
    StartTimeLCeq = temp2(idx)-temp5(idx-1)/Fs;
    
