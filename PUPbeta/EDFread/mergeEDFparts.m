clear all

filepath = ('C:\Users\LauraGell\Apnimed Dropbox\Apnimed\F -- R&D\Other combos\VicTor\08 Data and Reports\VicTor\Source\VICTOR-0011-023-PSG3\00100790-179440\');

filepath = 'C:\Users\LauraGell\Downloads\'
%load sections of edf
filetoload = 'LUN-1321-004-PSG1.edf';
filename = [filepath, filetoload ];
data1 = edfread(filename,'TimeOutputType','datetime');
info1 = edfinfo(filename);
info=info1


filetoload = 'SYN-1391-027-PSG0.edf';
filename = [filepath, filetoload ];
data2 = edfread(filename,'TimeOutputType','datetime');
info2 = edfinfo(filename);


% Rename duplicate labels
labels1 = matlab.lang.makeValidName(renameDuplicates(info1.SignalLabels));
labels2 = matlab.lang.makeValidName(renameDuplicates(info2.SignalLabels));

% Update variable names in data1 to match new labels
data1.Properties.VariableNames = (labels1);

% Update variable names in data2 to match new labels
data2.Properties.VariableNames = (labels2);

% Find extra channels in data1
extraLabels = setdiff(labels1, labels2);

% Add empty channels to data2 for extra labels
for i = 1:length(extraLabels)
    data2.(matlab.lang.makeValidName(extraLabels{i})) = repmat({zeros(100, 1)}, height(data2), 1);
end

% Ensure data2 has the same order of columns as data1
[~, order] = ismember(labels1, data2.Properties.VariableNames);
data2 = data2(:, order);


% Check and align data formats
for i = 1:size(data1,2)
    colName = data1.Properties.VariableNames{i};
    if ~isequal(class(data1.(colName)), class(data2.(colName)))
        % Convert data formats to match
        if iscell(data1.(colName)) && ~iscell(data2.(colName))
            data2.(colName) = num2cell(data2.(colName)); % Convert numeric to cell array
        elseif ~iscell(data1.(colName)) && iscell(data2.(colName))
            data2.(colName) = cell2mat(data2.(colName)); % Convert cell array to numeric
        end
    end
end


% Merge the data
datatowrite = [data1; data2];



clear data1 data2 data3 data4 data5

%%
%this adds missing data to timetable, fills with previous values repeated
%at base sample rate (so here get repeating last 1sec of data - replace
%with zeros
TT2 = retime(datatowrite,'secondly','previous') %'regular','SampleRate',1) 

TTbackup = TT2
TT2=data1

%fills any gaps with zeros - this may be an issue for some signals if far
%out of range
clear gapflag
for idx=1:height(TT2)
    gapflag = sum((TT2.("Record Time")(idx)==[datatowrite.("Record Time")]))==0;
    if gapflag
        for j =1:size(TT2,2)
            if iscell(TT2.(j)(idx))
                TT2.(j){idx} = zeros(size(TT2.(j){idx}));
            else
                TT2.(j)(idx) = zeros(size(TT2.(j)(idx)));
            end
        end
    end
end
TT2=TTbackup

clear sig temp
for i=1:size(TT2,2)
    temp=TT2{:,i};
    if iscell(temp)
    sig{i}=cat(1,temp{:});
    else
    sig{i}=temp(:);
    end
      
end


if 0
figure;
hold on
%plot(sig{1})
%plot(sig{15}(1:20:end))
plot(sig{19}(1:20:end))
TidyPlot()
end

test=sig{19};
for jj=2:length(test)
if test(jj)==0
test(jj)=test(jj-1);
end
end

sig{19}=test

%% 
figure;
hold on
%plot(sig{1})
%plot(sig{15}(1:20:end))
plot(sig{19}(1:20:end))
plot(test(1:20:end))
TidyPlot()
%%

%load the header info from first edf - min/max channel ranges seemed to be the same here so worked ok   
hdr = edfheader("EDF")
labels = fieldnames(hdr)
for i = 1:length(fieldnames(hdr))
 hdr.(labels{i})=info1.(labels{i})
end
hdr.NumDataRecords=height(TT2); %assign new file length
hdr.Reserved="";

new.PhysicalMin=info1.PhysicalMin;
new.PhysicalMax=info1.PhysicalMax;

for i=1:length(info1.PhysicalMin)
    if info1.PhysicalMin(i)>info1.PhysicalMax(i)
temp1=info1.PhysicalMin(i);
new.PhysicalMin(i)=info1.PhysicalMax(i);
new.PhysicalMax(i)=temp1;
    end
end
hdr.PhysicalMin=new.PhysicalMin
hdr.PhysicalMax=new.PhysicalMax


edfwrite([filepath 'newEDF2.edf'],hdr,sig,'InputSampleType' , "physical")

data1 = edfread([filepath 'newEDF2.edf'],'TimeOutputType','datetime');


