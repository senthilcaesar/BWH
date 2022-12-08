function BrDataTableOut=howfarawayfromsleep2(BrDataTable)
%distance from sleep is considered a minimum, since could also be next to a border. 
%%
BrDataTable.TimeJumpNext = abs([BrDataTable.Time_start(2:end);NaN] - BrDataTable.Time_end)>0.1;
BrDataTable.TimeJumpPast = abs([NaN;BrDataTable.Time_end(1:end-1)] - BrDataTable.Time_start)>0.1;
BrDataTable.TimeJumpPast(1)=1;
BrDataTable.TimeJumpNext(end)=1;

if sum(strcmp(BrDataTable.Properties.VariableNames, 'AR')) == 1
    a = BrDataTable.AR(:);  
else %legacy
    a = BrDataTable.ARei(:);  
end
A = cumsum(a);
B = [NaN;diff(a)==-1]; %reset at sleep boundary
Bwin = [BrDataTable.TimeJumpPast==1]; %reset at window boundary
I=find(B==1);
Iwin=find(Bwin==1);
Iboth = [I(:);Iwin(:)];
Itype = [ones(length(I),1);zeros(length(Iwin),1)];
array = [Iboth Itype]; 
array = sortrows(array);
Iboth = array(:,1);
Itype = array(:,2);
nleft=A;
for i=1:length(Iboth)
    if Itype(i)==1
        nleft(Iboth(i):end)=nleft(Iboth(i):end)-nleft(Iboth(i)); %reset to zero at i
    else
        nleft(Iboth(i):end)=nleft(Iboth(i):end)-(nleft(Iboth(i))-1);
    end
end
BrDataTable.nleft = nleft;

A2 = cumsum(a,'reverse');
B2 = [diff(a)==1;NaN]; %reset
B2win = [BrDataTable.TimeJumpNext==1]; %reset at window boundary
I=find(B2==1);
Iwin=find(B2win==1);
Iboth = [I(:);Iwin(:)];
Itype = [ones(length(I),1);zeros(length(Iwin),1)];
array = [Iboth Itype]; 
array = sortrows(array);
Iboth = array(:,1);
Itype = array(:,2);
nright=A2;
for i=length(Iboth):-1:1
    if Itype(i)==1
        nright(1:Iboth(i))=nright(1:Iboth(i))-nright(Iboth(i));
    else
        nright(1:Iboth(i))=nright(1:Iboth(i))-(nright(Iboth(i))-1);
    end
end
BrDataTable.nright = nright;
BrDataTable.nawayfromsleep = min([nleft nright]')';

BrDataTableOut = BrDataTable;

