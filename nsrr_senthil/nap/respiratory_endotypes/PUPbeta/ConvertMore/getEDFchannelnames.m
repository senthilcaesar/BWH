function getEDFchannelnames(Mrange)
global AMasterSpreadsheet handletext settings HChannelSpreadsheet
firstcol='BS';

isopen = xls_check_if_open(AMasterSpreadsheet); %the close function doesn't work, thus not called.
if isopen==1
    displaytext=['Warning: Spreadsheet must be closed for MATLAB to write data to it'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

[~,txt] = xlsread(AMasterSpreadsheet,1,'U4:AB10003');
maxN = 99; % DLM increased from 80 to 99, some UoW data with 93 channels

% Note from SS: "~isfield(settings,'RunPUPA2Z') && ~settings.RunPUPA2Z"
% code fails if ~isfield(settings,'RunPUPA2Z') is true, thus commented out.

N=size(txt,1);
 if ~exist('Mrange')==1
    Mrange = 1:N;
 end

Labels = cell(length(Mrange),maxN); % Assuming max maxN channels
Transducers = cell(length(Mrange),maxN);
Fss = nan(length(Mrange),maxN);

for n=1:length(Mrange)%1:N %size(txt,1)
    directory = char(txt(Mrange(n),4));
    fname=char(txt(Mrange(n),1));
    try
        if 0
            tic
            [Label,Transducer,Fs] = EDFChannelLabels([directory '\' fname]); %now fails in some MESA EDFs [2/2021]
            toc
        else
            tic
            [~,signalHeader] = blockEdfLoad(fname);
            signalHeaderT = struct2table(signalHeader);
            Label = signalHeaderT.signal_labels';
            Transducer = signalHeaderT.transducer_type';
            Fs = signalHeaderT.samples_in_record';
            toc
        end
        maxi = min([maxN size(Label,2)]);
        Labels(n,1:maxi) = Label(1:maxi);
        Transducers(n,1:maxi) = Transducer(1:maxi);
        Fss(n,1:maxi) = Fs(1:maxi);
        displaytext=['Reading n=' num2str(n) '/' num2str(N) ', ' fname];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
    catch me
        
        
        displaytext=['Reading n=' num2str(n) ', Error: ' me.message];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
    end
end

%Labels((maxN+1):end)=[];
displaytext=['Writing data to Xls'];
disp(displaytext); set(handletext,'String',displaytext); drawnow;


if isfield(settings,'LabelsOnly') && settings.LabelsOnly==1  %use Labels not Transducers
    xlswrite(AMasterSpreadsheet,Labels,1,[firstcol num2str(3+1+ Mrange(1)-1)]);
else
    LabelsTransducers=Labels;
    for i=1:size(Labels,1)
        for j=1:size(Labels,2)
            if ~isempty(LabelsTransducers{i,j})
                LabelsTransducers{i,j}=[Labels{i,j} '|' Transducers{i,j} '|' num2str(Fss(i,j))];
                LabelsTransducers{i,j} = regexprep(LabelsTransducers{i,j},' ','');
            end
        end
    end
    xlswrite(AMasterSpreadsheet,LabelsTransducers,1,[firstcol num2str(3+1+ Mrange(1)-1)]);
end


try
    if isfield(settings,'UseHarmonizedChannelNumbers')&& settings.UseHarmonizedChannelNumbers==1
        HChannelSpreadsheet=[settings.workdir 'PUPStart\HarmonizedChannelLabels.xlsx'];
        ChannelNumbers = nan(length(Mrange),24);
        for n=1:length(Mrange)
            temp = getChannelNumbersHarmonizedLabels([],HChannelSpreadsheet,Labels(n,:));
            temp(25:end)=[];
            ChannelNumbers(n,1:length(temp)) = temp;
            %write ChannelNumbers to Excel start from AJ..stop at BE--no pulse/ppg
        end
        ChannelNumbers2 = ChannelNumbers(:,1:22);
        ChannelNumbers3 = num2cell(ChannelNumbers2);
        ChannelNumbers3(isnan(ChannelNumbers2)) ={'NaN'};
        xlswrite(AMasterSpreadsheet,ChannelNumbers3,1,['AJ' num2str(4+Mrange(1)-1)]);
        
    end
end

winopen(AMasterSpreadsheet);

end

