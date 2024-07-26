function [t,x,y,a,win,nota,veup,hyp,pcw,t2,e,GGp,GGt] = VEPesArray(DataOut,n,criteria,drivecol)
%%
%build array of ventilation and ventilatory drive data specific to PUP format
n=1;

NExcludepostAR = 2;

VEcol = 15;
Veupcol = 11;
Errorcol = 9;
ARcol = 10;
tcol = 1;
t2col = 3;
hypcol = 14; 
pescol = drivecol; %16 [16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak]
pmuscol = 17; %17
Ecol = 6;
GGpcol = 21;
GGtcol = 22;


x=[]; %ventilatory drive, normalized locally
y=[]; %ventilation, normalized locally
t=[];
t2=[];
a=[];
nota=[]; 
win=[];
veup=[];
hyp=[];
pcw=[];
e=[];
GGp=[];
GGt=[];
%subject = n;

for w=1:length(DataOut{n})
    if (size(DataOut{n}{w},1)==1&&isnan(DataOut{n}{w}))||isempty(DataOut{n}{w})||criteria(w)==0
        continue
    end
    GGp1=DataOut{n}{w}(:,GGpcol);
    GGt1=DataOut{n}{w}(:,GGtcol);
    y1=DataOut{n}{w}(:,VEcol);
    e1=DataOut{n}{w}(:,Ecol);
    x1=DataOut{n}{w}(:,pescol);%-DataOut{n}{w}(:,Errorcol);
    pcw1=DataOut{n}{w}(:,pmuscol) - DataOut{n}{w}(:,pescol);%-DataOut{n}{w}(:,Errorcol);
    t1=DataOut{n}{w}(:,tcol);
    t21=DataOut{n}{w}(:,t2col);
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
    
    GGp = [GGp;GGp1];
    GGt = [GGt;GGt1];
    x = [x;x1];
    e = [e;e1];
    y = [y;y1];
    t = [t;t1];
    t2 = [t2;t21];
    a = [a;a1];
    win = [win;win1];
    nota=[nota;nota1];
    veup=[veup;veup1];
    hyp = [hyp;hyp1];
    pcw = [pcw;pcw1];
if 0
figure(100)
stairs(t1,[x1 y1 2+0.3*a1]);
hold('on')
end
end
t2 = t2-t;

