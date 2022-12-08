function [n,nleft,nright]=howfarawayfromsleepInf(a,win)
%distance from sleep is considered a minimum, since could also be next to a border. 

if 0 %workable example
a=[1 0 0 1 1 1 0 0 0 1 1 1 1 1 1 1 0]';
win=[1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3]';

a=[1 1 1 1 0 0 1 1 1 0 0 0 1 1 1 1 1 1 1 0 1 1 1 1]';
win=a*0;
end

A = cumsum(a);
B = [NaN;diff(a)==-1]; %reset at sleep boundary
Bwin = [NaN;diff(win)~=0]; %reset at window boundary
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
        nleft(Iboth(i):end)=nleft(Iboth(i):end)-nleft(Iboth(i));
    else
        nleft(Iboth(i):end)=nleft(Iboth(i):end)-(nleft(Iboth(i))-1);
    end
end
if ~isempty(Iboth)
    nleft(1:(Iboth(1)-1))=Inf;
end

A2 = cumsum(a,'reverse');
B2 = [diff(a)==1;NaN]; %reset
B2win = [diff(win)~=0;NaN]; %reset at window boundary
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
if ~isempty(Iboth)
    nright((Iboth(end)+1):end)=Inf;
end

n=min([nleft nright]')';


