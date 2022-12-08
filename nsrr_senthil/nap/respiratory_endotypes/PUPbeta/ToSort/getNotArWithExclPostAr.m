function nota = getNotArWithExclPostAr(a,win,NExcludepostAR)

%used to select breaths during sleep, i.e. nota==1

nota=1-a;

newwinstart = 1*([NaN;diff(win)]~=0);
Inewwin = [find(newwinstart==1)];
if length(Inewwin)>0
    for i=1:length(Inewwin)
        li=Inewwin(i);
        ri=li+NExcludepostAR-1;
        if ri>length(a), ri=length(a); end
        nota(li:ri)=0;
    end
end


aoffset = -[NaN;diff(a)];
I=find(aoffset==1);
if length(I)>0
    for i=1:length(I)
        li=I(i);
        ri=li+NExcludepostAR-1;
        if ri>length(a), ri=length(a); end
        newwindata = newwinstart(li:ri);
        if sum(newwindata)>1
            ri=find(newwindata==1)-2+li;
        end
        nota(li:ri)=0;
    end
end

temp = [a nota aoffset win newwinstart];