function [t x y a win nota veup hyp e x2] = VEVdriveArray(DataOut,n,criteria)
%%
%build array of ventilation and ventilatory drive data specific to PUP format

NExcludepostAR = 2;

VEcol = 4;
Veupcol = 11;
Errorcol = 9;
ARcol = 10;
tcol = 1;
hypcol = 14;
Vchemcol = 5;
Ecol = 6;
Vdrivetotal = 12;

x=[]; %ventilatory drive, normalized locally
y=[]; %ventilation, normalized locally
t=[];
a=[];
nota=[]; 
win=[];
veup=[];
hyp=[];
e=[];
y2=[];

subject = n;
for w=1:length(DataOut{n})
    if (size(DataOut{n}{w},1)==1&&isnan(DataOut{n}{w}))||isempty(DataOut{n}{w})||criteria(w)==0
        continue
    end
    y1=DataOut{n}{w}(:,VEcol);
    y21=DataOut{n}{w}(:,Vdrivetotal);
    e1=DataOut{n}{w}(:,Ecol);
    x1=DataOut{n}{w}(:,Vchemcol);
    t1=DataOut{n}{w}(:,tcol);
    a1=ceil(DataOut{n}{w}(:,ARcol));
    hyp1=DataOut{n}{w}(:,hypcol);
    nota1=1-a1;
        nota1(1:NExcludepostAR)=0;
        aoffset = -[NaN;diff(a1)];
        I=find(aoffset==1);
        if length(I)>0
            for i=1:length(I)
                li=I(i);
                ri=I(i)+NExcludepostAR-1;
                if ri>length(x1), ri=length(x1); end
                nota1(li:ri)=0;
            end
        end
        
    win1 = w+0*x1;
    veup1 = DataOut{n}{w}(:,Veupcol);
    
    x = [x;x1];
    e = [e;e1];
    y = [y;y1];
    y2 = [y2;y21];
    t = [t;t1];
    a = [a;a1];
    win = [win;win1];
    nota=[nota;nota1];
    veup=[veup;veup1];
    hyp = [hyp;hyp1];
if 0
figure(100)
stairs(t1,[x1 y1 2+0.3*a1]);
hold('on')
end


end
x=x+1; %add eupnea back
x2=x+y2; %vdrivetotal

