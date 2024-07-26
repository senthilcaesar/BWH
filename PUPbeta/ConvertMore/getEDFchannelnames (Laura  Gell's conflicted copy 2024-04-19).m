function getEDFchannelnames(Mrange)
global AMasterSpreadsheet handletext settings 
firstcol='BS';
%Mrange presumably needs to be a contiguous block e.g. 2:5, not [2 5 7].
%thus do this:


isopen = xls_check_if_open(AMasterSpreadsheet); %the close function doesn't work, thus not called.
if isopen==1
    displaytext=['Warning: Spreadsheet must be closed for MATLAB to write data to it'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

[~,txt] = xlsread(AMasterSpreadsheet,1,'U4:AB15000');%AB10003
maxN = 112; % JEO increased to 115, Nox data with 112 channels

% Note from SS: "~isfield(settings,'RunPUPA2Z') && ~settings.RunPUPA2Z"
% code fails if ~isfield(settings,'RunPUPA2Z') is true, thus commented out.

 if exist('Mrange')==1
     Mrange = [min(Mrange):1:max(Mrange)];
 end

 N=size(txt,1);
 if ~exist('Mrange')==1
    Mrange = 1:N;
 end

 Mrange(Mrange>length(settings.ConvertedFilename))=[];

Labels = cell(length(Mrange),maxN); Labels(:)={''};
% Assuming max maxN channels
Transducers = cell(length(Mrange),maxN); Transducers(:)={''};
Fss = nan(length(Mrange),maxN);

%% RMA-added partial path for edf channel check-9/12/2022
% earlier if only partial path was given edf channel check was failing

%detect if Source data directories are full or partial paths
try
    temp = nan(length(txt),1);
    for i=1:size(txt,1)
        temp(i)=any(txt{i,4}==':');
    end
    SourceDataIsPartialPath=nanmean(temp)<0.5;
    if SourceDataIsPartialPath==1 && ~(isfield(settings,'FileNameNoOverwrite')&& settings.FileNameNoOverwrite==1) %if not a full path (missing ":"), then do this:
        for j=4:6
            for i=1:size(txt,1)
                if isnan(txt{i,j})
                    txt{i,j} = settings.SourceDirectory; %write sourceDirectory to Filenames cols 4 to 6
                    if ~isfile([settings.SourceDirectory txt{i,1}]) && isfield(settings,'SourceSubfolders') && settings.SourceSubfolders==1 % Dan Vena added for edf files in individual folders
                        txt{i,j} = [settings.SourceDirectory,settings.patients{i,1}(1:end-4), filesep]; % Samer: 2023-09-14: used filesep instead of '\'
                    end
                else
                    txt{i,j} = [settings.SourceDirectory txt{i,j}]; %concat sourceDirectory to Filenames cols 4 to 6
                end
            end
        end
    end
end
%% run EDFChannelLabels

    disp('Warning: Using EDFChannelLabels function [fast], but has failed in read of MESA EDFs; switch to blockEDFLoad method if needed [slow]');
    
for n=1:length(Mrange)%1:N %size(txt,1)
    directory = char(txt(Mrange(n),4));
    fname     = char(txt(Mrange(n),1));
    try
        if 1 %this way needs to be default; other is too slow; add setting if needed or can work to revise the code to make the fread perform properly
            %tic
            [Label,Transducer,Fs] = EDFChannelLabels([directory fname]); %confirmed ok in MESA (20230728) using R2021a
            %toc
        else %backup but slow
            %tic
            [~,signalHeader] = blockEdfLoad([directory fname]);
            signalHeaderT = struct2table(signalHeader);
            Label = signalHeaderT.signal_labels';
            Transducer = signalHeaderT.transducer_type';
            Fs = signalHeaderT.samples_in_record'; %%%% this line is incorrect
            %toc
        end
        
        %Tidy up (20230727)
        Label = strtrim(Label);
        Transducer = strtrim(Transducer);
        try %JEO strtrim not consistently removing all whitespaces
            Label = deblank(Label);
            Transducer = deblank(Transducer);
        end
        try
            Transducer=strip(Transducer,'.'); %SS what does this do? It is breaking.
        end
        
        maxi = min([maxN size(Label,2)]);
        Labels(n,1:maxi) = Label(1:maxi);
        Transducers(n,1:maxi) = Transducer(1:maxi);
        Fss(n,1:maxi) = Fs(1:maxi);
        displaytext=['Reading n=' num2str(Mrange(n)) '/' num2str(N) ', ' fname];
        % disp(displaytext); set(handletext,'String',displaytext); drawnow;
        disp(displaytext); % set(handletext,'String',displaytext); drawnow;
        % Samer Bou Jawde: Oct. 23 2023 - code was breaking at
        %   set(handle..) so commented it out. Is it necessary?
        
    catch me
        
        displaytext=['Failed n=' num2str(Mrange(n)) ', Error: ' me.message];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
    end
end

%% Writing data to Xls
%Labels((maxN+1):end)=[];
displaytext = ['Writing data to Xls'];
disp(displaytext); set(handletext,'String',displaytext); drawnow;

if isfield(settings,'LabelsOnly') && settings.LabelsOnly==1  %use Labels not Transducers
    %xlswrite(AMasterSpreadsheet,Labels,1,[firstcol num2str(3+1+ Mrange(1)-1)]);
    writecell(Labels,AMasterSpreadsheet,'Sheet',1,'Range',[firstcol num2str(3+1+ Mrange(1)-1)])
else
    LabelsTransducers=Labels;
    for i=1:size(Labels,1)
        for j=1:size(Labels,2)
            if ~isempty(LabelsTransducers{i,j})
                LabelsTransducers{i,j}=[Labels{i,j} '|' Transducers{i,j} '|' num2str(Fss(i,j))];
                %LabelsTransducers{i,j} = regexprep(LabelsTransducers{i,j},' ','');
            end
        end
    end
   try
        writecell(LabelsTransducers,AMasterSpreadsheet,'Sheet',1,'Range',[firstcol num2str(3+1+ Mrange(1)-1)]) %% RMA-this is not working in somecases (11062023)
   catch
        'using backup write method'
        xlswrite(AMasterSpreadsheet,LabelsTransducers,1,[firstcol num2str(3+1+ Mrange(1)-1)]);
   end
end

%% UseHarmonizedChannelNumbers = 1 or 2
try
    if isfield(settings,'UseHarmonizedChannelNumbers') && settings.UseHarmonizedChannelNumbers>0
       
        Nmax=size(settings.ChannelNumbers,2);
        ChannelNumbers = getChannelNumbersFromHarmonizedLabels(settings.HChannelSpreadsheet,Mrange,Labels,Transducers,Nmax)
        
        %xlswrite(AMasterSpreadsheet,ChannelNumbers3,1,['AJ' num2str(4+Mrange(1)-1)]);
        disp('Saving Harmonized Channel Label Results');
        ChannelNumbers3 = num2cell(ChannelNumbers);
        ChannelNumbers3(isnan(ChannelNumbers))={'NaN'};
        try
         writecell(ChannelNumbers3,AMasterSpreadsheet,'Sheet',1,'Range',['AJ' num2str(4+Mrange(1)-1)])
        catch
         xlswrite(AMasterSpreadsheet,ChannelNumbers3,1,['AJ' num2str(4+Mrange(1)-1)]);
        end
        MissingFlowList = string(settings.ConvertedFilename(find(isnan(ChannelNumbers(:,1)))))
        LTerrFlow = Labels(find(isnan(ChannelNumbers(:,1))),:);
        LTerrThorax = Labels(find(isnan(ChannelNumbers(:,2))),:);
        LTerrAbdomen = Labels(find(isnan(ChannelNumbers(:,3))),:);
        LTerrEKG = LabelsTransducers(find(isnan(ChannelNumbers(:,19))),:);
        LTerrC3 = LabelsTransducers(find(isnan(ChannelNumbers(:,5))),:);
        LTerrC4 = LabelsTransducers(find(isnan(ChannelNumbers(:,7))),:);
        LTerrO1 = LabelsTransducers(find(isnan(ChannelNumbers(:,9))),:);
        LTerrO2 = LabelsTransducers(find(isnan(ChannelNumbers(:,9))),:);
        LTerrF3 = LabelsTransducers(find(isnan(ChannelNumbers(:,11))),:);
        LTerrF4 = LabelsTransducers(find(isnan(ChannelNumbers(:,13))),:);
        LTerrSpO2 = LabelsTransducers(find(isnan(ChannelNumbers(:,4))),:);
        MissingSpO2List = string(settings.ConvertedFilename(find(isnan(ChannelNumbers(:,4)))))
    end
end
% to do- concaatenate label transducer

%winopen(AMasterSpreadsheet);

end

% to do, rerun ImportSettings so that you don't have to rerun StartHere all
% over again after a Channel Check