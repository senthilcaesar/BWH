function [t x y z a win nota veup hyp pcw] = GGPesArray(DataOut,n,criteria,drivecol)
%%
%build array of ventilation and ventilatory drive data specific to PUP format

NExcludepostAR = 2; %nota is ~a with further exclusion for immediate post sleep onset breaths.

GGpeakcol = 21; %this is GG peak...
GGtoniccol = 22; %this is GG tonic...
Veupcol = 11;
Errorcol = 9;
ARcol = 10;
tcol = 1;
hypcol = 14; 
pescol = drivecol; %16 [16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak]
pmuscol = 17; %17

x=[]; %ventilatory drive, normalized locally
y=[]; %ventilation, normalized locally
z=[]; %ventilation, normalized locally
t=[];
a=[];
nota=[]; 
win=[];
veup=[];
hyp=[];
pcw=[];

subject = n;
for w=1:length(DataOut{n})
    if (size(DataOut{n}{w},1)==1&&isnan(DataOut{n}{w}))||isempty(DataOut{n}{w})||criteria(w)==0
        continue
    end
    y1=DataOut{n}{w}(:,GGpeakcol);
    z1=DataOut{n}{w}(:,GGtoniccol);
    x1=DataOut{n}{w}(:,pescol);%-DataOut{n}{w}(:,Errorcol);
    pcw1=DataOut{n}{w}(:,pmuscol) - DataOut{n}{w}(:,pescol);%-DataOut{n}{w}(:,Errorcol);
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
    y = [y;y1];
    z = [z;z1];
    t = [t;t1];
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

