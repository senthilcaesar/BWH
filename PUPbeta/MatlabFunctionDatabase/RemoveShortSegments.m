function Out = RemoveShortSegments(In,MinNoiseTime,dT,enableremovegood)
%Signal input is a series of 1s and 0s; 1 indicates artifact (bad); 0 indicates good signal
%Goal is to remove short segments of bad or good signals
%Method removes segments progressively in order from shortest to longest, not time order, preserving length

removeshorterthani = MinNoiseTime/dT; % set the duration threshold to remove

% turn enableremovegood on by default (e.g. if missing)
if ~exist('enableremovegood')
    enableremovegood=1;
end

if enableremovegood==0
    disp('enableremovegood==0, brief testing only, please check data')
end

I = diff([NaN;In]);
I1 = find(I==1);
I2 = find(I==-1);
[I1Ones,I2Ones]=TidyStartEndEventList(I1,I2,length(In));

if isempty(I1Ones) %no bad data case
    Out=0*In;
    return
end

% remove extremely small gaps, bad only
I = find((I2Ones-I1Ones)<=ceil(removeshorterthani/1000)); %threshold should be at least one sample to be useful.
I1Ones(I)=[];
I2Ones(I)=[];

Tindex = table(I1Ones,ones(length(I1Ones),1)); Tindex.Properties.VariableNames={'Index','Art'};
TindexA = table(I2Ones,zeros(length(I2Ones),1)); TindexA.Properties.VariableNames={'Index','Art'};
Tindex = [Tindex;TindexA];
Tindex=sortrows(Tindex);

if Tindex.Index(1)~=1 % edge case - first entry in table
    Tindex = [Tindex(1,:);Tindex];
    Tindex{1,:} = [1 0];
    %add row with 1st row at state [1 0]
end

if Tindex.Index(end)~=length(In) % edge case - last entry in table
    Tindex = [Tindex;Tindex(end,:)];
    Tindex{end,:} = [length(In) NaN];
    %add row with 1st row at state [1 0]
end
Tindex.Art(end)=NaN;
Tindex.Length=[diff(Tindex.Index);NaN]*dT;
if ~enableremovegood
    Tindex.Length(Tindex.Art==0)=Inf;
end

while 1
    %if enableremovegood
    if height(Tindex)<=2 %if only one long segment, break
        break
    end
    if min(Tindex.Length)>MinNoiseTime %if no short segments, break
        break
    end
    
    [~,i]=min(Tindex.Length);
    if i>1 && i<(height(Tindex)-1)
        Tindex([i:(i+1)],:)=[];
    elseif i==1
        Tindex(i,:)=[];
        Tindex.Index(1)=1;
    elseif i==(height(Tindex)-1)
        Tindex(i,:)=[];
    end
    Tindex.Length=[diff(Tindex.Index);NaN]*dT;
    if ~enableremovegood
        Tindex.Length(Tindex.Art==0)=Inf;
    end
end

Out=0*In;
I = find(Tindex.Art==1);
for i=1:length(I)
    Out(Tindex.Index(I(i)):Tindex.Index(I(i)+1)-1)=1;
end






