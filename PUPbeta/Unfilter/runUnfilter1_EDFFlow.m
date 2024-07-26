%function tryUnfilterFlow()

%% global variables and settings
global AMasterSpreadsheet handletext F_samp settings

%% master spreadsheet
convertfilename = [settings.workdir, filesep, 'AMasterSpreadsheet.xlsx'];
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = convertfilename;
end

%% Start processing

% read spreadsheet (Files worksheet)
[~,~,raw] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
disp('Warning: XLS cols BD, BE now point to LOC, ROC (eyes); BF=Every2ndEEGref; update xls if needed');

InvertFlow = cell2mat(raw(:,14));

ChannelNumbers_ = raw(:,17:end-1);
ChannelNumbers = NaN*ones(size(ChannelNumbers_));
for i = 1:size(ChannelNumbers_,1)
    for j = 1:size(ChannelNumbers_,2)
        if ischar(ChannelNumbers_{i,j})
            ChannelNumbers(i,j)=NaN;
            %ChannelNumbers_{i,j}=NaN;
        else
            ChannelNumbers(i,j)=ChannelNumbers_{i,j};
        end
    end
end

%last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
lastrow = find(1*(~isnan(ChannelNumbers(:,1))),1,'last');
raw(lastrow+1:end,:)=[];
ChannelNumbers(lastrow+1:end,:)=[];

Filenames = raw(:,2:8);

% read spreadsheet (Options worksheet)
[~,~,options] = xlsread(AMasterSpreadsheet,2,'C3:C11');
runrange=(1:size(Filenames,1));

%%
clear T
AlternativeFound = nan(max(runrange),1);
T = table(AlternativeFound);
T.ChannelNumberA = nan(max(runrange),1);
T.ChannelNumberB = nan(max(runrange),1);
T.ChannelNumberC = nan(max(runrange),1);
%%
for n=runrange
    errorlist=[];
   
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    EDFfilenamedir = [directory fname];
    %%
    if ~exist(EDFfilenamedir)
        continue
    end
    %%   
        fname = EDFfilenamedir;
fid = fopen(fname,'r');
fseek(fid,252,-1);
clear Label Transducer Fs
M  = sscanf(char(fread(fid,4,'char')'),'%d');
for m=1:M
    fseek(fid,(m-1)*16+256,-1);
    Label{m} = char(fread(fid,16,'char')'); %%%%%%%%%
end

for m=1:M
    fseek(fid,(m-1)*80+256+M*16,-1); 
    Transducer{m} = fread(fid,80,'*char')';
end

fseek(fid,244,-1); 
Block = str2double(fread(fid,8,'*char')');
for m=1:M
    fseek(fid,(m-1)*8+256+M*16+M*80+M*120,-1); 
    Fs(m) = str2double(fread(fid,8,'*char')')/Block;
end
fclose(fid); % Close file
   
                
        T.ChannelNumberA(n) = ChannelNumbers(n,1);
        T.FlowLabelA{n} = Label{T.ChannelNumberA(n)};
        T.TransducerLabelA{n} = Transducer{T.ChannelNumberA(n)};
        I = find(strcmp(Label{T.ChannelNumberA(n)},Label));
        I(I==T.ChannelNumberA(n))=[];
        %%
        if isempty(I)
            disp(['No alternative Channel labelled: ' T.FlowLabelA{n}]);
            T.AlternativeFound(n)=0; 
            continue
        else
            T.AlternativeFound(n)=1;
            T.ChannelNumberB(n) = I(1);
            T.FlowLabelB{n} = Label{T.ChannelNumberB(n)};
            T.TransducerLabelB{n} = Transducer{T.ChannelNumberB(n)};
            if length(I)>1
                T.ChannelNumberC(n) = I(2);
                T.FlowLabelC{n} = Label{T.ChannelNumberC(n)};
                T.TransducerLabelC{n} = Transducer{T.ChannelNumberC(n)};
            end
        end
        
        
end

%%
      testrange  = 11;
for n=testrange
        
    
    %%
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    EDFfilenamedir = [directory fname];
    
    %%
    if ~exist(EDFfilenamedir)
        continue
    end
    %%
    figure(899); clf(899); set(gcf,'color',[1 1 1]);
    
    [x,fs,~,~,label,~,~,~,~] = readedfrev3(EDFfilenamedir,T.ChannelNumberA(n)-1,0,Inf);
         t=(0:(1/fs):0+(length(x)-1)*(1/fs))'; % This is the time vector associated with the _XHz Flow data.
    if InvertFlow(n)
        x=-x;
    end
    
    %attempt unfilter...
    x = unfilter1(x,fs,1/30,2,1,[1/60],1);
    
    %
    %xlims = get(gca,'xlim');
    %set(gca,'xlim',xlims);
    
    %%
    %x = unfilter1(x,fs,1/2,1,1,[]);
    tminmax = get(gca,'xlim');
    
    Fstart = 0.5;
    xlim([tminmax(2)*Fstart tminmax(2)*Fstart+600]);
    
    %%
end

%%

