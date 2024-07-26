%% UA phenotype info using Edi/Pes
if ~exist('T')
    T = LargeTableFromCellsOfTables(BreathDataTable);
end
%%

dividedrivebyttot=0;
normalizeyusingveupwindow=1;
normalizeyusingknownVeupnea=0;
G=[];
gobackN2 = 5; %3 for application
gobackN = 10; %3 for application
%get wake data next (27,23,17) 
EdiPrearArray_N=[];
clear Delta
EdiPrearSubj_N=[];
PrearArray_Pbeta_N=[];
PrearArray_Palpha_N=[];
PrearArray_Ptheta_N=[];
PrearArray_Pdelta_N=[];
PrearArray_SpO2_N=[];
PrearArray_VE_N=[];


PreVEsleepArray_Edi_N=[];
PreVEsleepArray_VE_N=[];
PreVEsleepArray_Subj_N=[];
for n=1:M
    if Exclude(n)==1
        continue
    end
    %%
    criteria=ones(length(DataOut{n}),1);
    
    %window-based "criteria" are not being applied here yet...
    T_n = T(T.subj==n,:);
    t = T_n.Time_start;
	x = T_n.DeltaEdi;
    y = T_n.VE; %check what is VI
    veup = T_n.Veup;
    win = T_n.win;
    a = T_n.ARei; %Is actually ARie I think
    hyp = T_n.hypnog_B;
    t2 = T_n.Time_end-t; 
    e = T_n.Etype;
    
    s = T_n.spo2;
    Pbeta = T_n.Pbeta;
    Palpha = T_n.Palpha;
    Pdelta = T_n.Pdelta;
    Ptheta = T_n.Ptheta;
    
    if 1
        a(e>1)=0;
    end
    nota = getNotArWithExclPostAr(a,win,2);
    
    %[t,x,y,a,win,nota,veup,hyp,pcw,t2] = VEPesArray(DataOut,n,criteria,19);
    [tL,xL,yL,aL,winL,notaL,veupL,hypL,pcwL,t2L,eL] = VEPesArray(DataOut,n,criteria,19);
    
    temp = [a aL];
    
    %checklengths(n)=length(t)==length(tL);
    checkcontent(n)=sum(a(~isnan(a))~=aL(~isnan(aL)));
    
    if dividedrivebyttot
        x=x./t2;
    end
    [ad]=howfarawayfromsleep(a,win); %replace a with (a==1)&(hyp==4) if arousals not scored in wake...
    
    %find breaths during wakefulness and arousals
    if 1
        minNwakebreaths = 50;
        hh=hist(ad,[1:11]); hh(end)=NaN; th=find(hh>minNwakebreaths,1,'last');
        threshold(n) = min([4 th]);
    else
        threshold(n) = 2;
    end
    a1 = ad>threshold(n);
    G(n) = nanmedian(y(a1)./x(a1));
    if 1
        plotfig=0;
        if plotfig
            figure()
            plot(t(ad>=2),y(ad>=2)./x(ad>=2),'r.');
            hold('on');
            plot(t(a1),y(a1)./x(a1),'.');
            plot(t,nanmedian(y(ad>=2)./x(ad>=2))+0*t);
        end
        temp = t(ad>=2);
        data = y(ad>=2)./x(ad>=2);
        [temp,ia,ic] = unique(temp);
        data = data(ia);
        data(data>prctile(data,99))=NaN;
        data(data<prctile(data,1))=NaN;
        temp(isnan(data))=[]; data(isnan(data))=[];
        maxt = 1800;
        Gmovingave=NaN*temp;
        for j=1:length(temp)
            temp2 = abs(temp(j)-temp);
            weights = maxt-temp2; weights(weights<0)=0; weights = weights/sum(weights);
            Gmovingave(j)=nansum(data(weights>0).*weights(weights>0));
        end
        temp3=sortrows([temp,Gmovingave]);
        Gtime{n}=temp3;
        if plotfig
            plot(Gtime{n}(:,1),Gtime{n}(:,2),'g-');
        end
    end
    
       g = G(n);
        
       %Continuous time drive calibration
            Gtemp=interp1(Gtime{n}(:,1),Gtime{n}(:,2),t);
            if isnan(Gtemp(1))
                I=find(~isnan(Gtemp),1,'first');
                Gtemp(1:I-1)=Gtemp(I);
            end
            if isnan(Gtemp(end))
                I=find(~isnan(flipud(Gtemp)),1,'first');
                Gtemp((end-I+2):end)=Gtemp(end-I+1);
            end
            
            x_lps=x.*Gtemp;
            
%             try
%             y_=nanmedian(veup);
%             x_ = y_/g;
%             figure(9)
%             [Vpassive(n),Vactive(n),~,~,~,~]=VEVdriveAnalysis(x,y,arthres(n),1,x_,10);
%             catch me
%             end
   
       ArOn=[NaN;diff(a)];
       I = find(ArOn==1);
   
       I((I-gobackN)<1)=[]; %remove if too close to start
       I(win(I-gobackN)~=win(I))=[]; %remove if N breaths back is in prior window
       for i=length(I):-1:1
            if sum(a((I(i)-gobackN2):I(i)-1)==1)~=0 %if not enough sleep prior to arousal
                I(i)=[];
            end
       end
       EdiPrearArray = zeros(length(I),gobackN+1);
       PrearArray_Pbeta = zeros(length(I),gobackN+1);
       PrearArray_Palpha = zeros(length(I),gobackN+1);
       PrearArray_Pdelta = zeros(length(I),gobackN+1);
       PrearArray_Ptheta = zeros(length(I),gobackN+1);
       PrearArray_SpO2 = zeros(length(I),gobackN+1);
       PrearArray_VE = zeros(length(I),gobackN+1);
       for i=1:length(I)
           irange = (I(i)-gobackN):I(i);
           PrearArray_Pbeta(i,:) = Pbeta(irange);
           PrearArray_Palpha(i,:) = Palpha(irange);
           PrearArray_Pdelta(i,:) = Pdelta(irange);
           PrearArray_Ptheta(i,:) = Ptheta(irange);
           PrearArray_SpO2(i,:) = s(irange);
           EdiPrearArray(i,:) = x_lps(irange);
           PrearArray_VE(i,:) = y(irange)./veup(irange);
       end
       EdiPrearVeup = veup(I);     
       EdiPrearSubj = n+0*EdiPrearVeup;
       
       EdiPrearSubj_N = [EdiPrearSubj_N;EdiPrearSubj];
       EdiPrearArray_N = [EdiPrearArray_N;EdiPrearArray/median(EdiPrearVeup)];     
       PrearArray_Pbeta_N = [PrearArray_Pbeta_N;PrearArray_Pbeta];     
       PrearArray_Palpha_N = [PrearArray_Palpha_N;PrearArray_Palpha];     
       PrearArray_Pdelta_N = [PrearArray_Pdelta_N;PrearArray_Pdelta];     
       PrearArray_Ptheta_N = [PrearArray_Ptheta_N;PrearArray_Ptheta];
       PrearArray_SpO2_N = [PrearArray_SpO2_N;PrearArray_SpO2];
       PrearArray_VE_N = [PrearArray_VE_N;PrearArray_VE];
       
       % matrix of VE during sleep. todo: add in sleep state info
       I = find(nota==1);   
       I((I-gobackN)<1)=[]; %remove if too close to start
       I(win(I-gobackN)~=win(I))=[]; %remove if N breaths back is in prior window
       PreVEsleepArray_VE = zeros(length(I),gobackN+1);
       PreVEsleepArray_Edi = zeros(length(I),gobackN+1); 
       for i=1:length(I)
           irange = (I(i)-gobackN):I(i);
           PreVEsleepArray_Edi(i,:) = x_lps(irange)./veup(irange);
           PreVEsleepArray_VE(i,:) = y(irange)./veup(irange);
       end
       PreVEsleepArray_VE_N = [PreVEsleepArray_VE_N;PreVEsleepArray_VE];
       PreVEsleepArray_Edi_N = [PreVEsleepArray_Edi_N;PreVEsleepArray_Edi];
       PreVEsleepArray_Subj_N = [PreVEsleepArray_Subj_N;n+0*I];
end    

%remove rows with NaN
IremoveNaN=isnan(sum(EdiPrearArray_N'));
EdiPrearArray_N(IremoveNaN,:)=[];
EdiPrearSubj_N(IremoveNaN,:)=[];
PrearArray_Pbeta_N(IremoveNaN,:)=[];
PrearArray_Palpha_N(IremoveNaN,:)=[];
PrearArray_Pdelta_N(IremoveNaN,:)=[];
PrearArray_Ptheta_N(IremoveNaN,:)=[];
PrearArray_SpO2_N(IremoveNaN,:)=[];
PrearArray_VE_N(IremoveNaN,:)=[];

IremoveNaN=isnan(sum((PreVEsleepArray_Edi_N+PreVEsleepArray_VE_N)'));
PreVEsleepArray_VE_N(IremoveNaN,:)=[];
PreVEsleepArray_Edi_N(IremoveNaN,:)=[];
PreVEsleepArray_Subj_N(IremoveNaN,:)=[];
% 
% predictI=8;
% predictstart=6;
% X5=[EdiPrearArray_N(:,predictstart:predictI-1)];
% y5=EdiPrearArray_N(:,predictI);
% 
% predictI=7;
% predictstart=5;
% X4=[EdiPrearArray_N(:,predictstart:predictI-1)];
% y4=EdiPrearArray_N(:,predictI);
% 
% predictI=6;
% predictstart=4;
% X3=[EdiPrearArray_N(:,predictstart:predictI-1)];
% y3=EdiPrearArray_N(:,predictI);
%%

%%
% predictI=6;
% predictstart=3;
% X2=[EdiPrearArray_N(:,predictstart:predictI-1)];
% y2=EdiPrearArray_N(:,predictI);
if 1
predictI=5+gobackN-gobackN2;
predictstart=2+gobackN-gobackN2;
X1=[EdiPrearArray_N(:,predictstart:predictI-1)];
y1=EdiPrearArray_N(:,predictI);
% 
predictI=4+gobackN-gobackN2;
predictstart=1+gobackN-gobackN2;
X0=[EdiPrearArray_N(:,predictstart:predictI-1)];
y0=EdiPrearArray_N(:,predictI);

X=[X0;X1];
%X = [X ones(size(X,1),1)];
%X = [X(:,2)-X(:,1) X(:,2)];
%X = [X(:,end)-X(:,end-1) X(:,end)];
%slope = 1/6*(X(:,end)-X(:,end-2)) + 1/3*(X(:,end)-X(:,end-1)) + 1/3*(X(:,end-1)-X(:,end-2));
%slope = [(X(:,end)-X(:,end-1)) X(:,end)-X(:,end-2)];
%slope = [(X(:,end)-X(:,end-1))];

%X=X(:,end);
%X = [X(:,end) slope];
y=[y0;y1];

A = X\y;
[b,bint,r,rint,stats] = regress(y,X);
stats(1)
else

b=[ 0.1306 ...
    0.2446 ...
    0.6737]'
end
%%
m = length(b);
X = EdiPrearArray_N(:,end-m:end-1);
% slope = 1/6*(X(:,end)-X(:,end-2)) + 1/3*(X(:,end)-X(:,end-1)) + 1/3*(X(:,end-1)-X(:,end-2));
% slope = 0.5*(X(:,end)-X(:,end-2));
% slope = [(X(:,end)-X(:,end-1)) X(:,end)-X(:,end-2)];
% slope = [(X(:,end)-X(:,end-1))];
% X = [X(:,end) slope];
%X=X(:,end);
%b = b/sum(b);
if m>1
Predict = sum((X.*b')')';
else
Predict = X.*b;
end
Delta = EdiPrearArray_N(:,end)-Predict;

for i=1:M
    VRAmedian(i) = 100*nanmedian(Delta(EdiPrearSubj_N==i));
    VRAmean(i) = 100*nanmean(Delta(EdiPrearSubj_N==i));
    VRAstd(i) = 100*nanstd(Delta(EdiPrearSubj_N==i));
    VRAn(i) = nansum((EdiPrearSubj_N==i));
end
VRAlowerCI = VRAmean-1.96*(VRAstd./(VRAn.^.5));
VRAupperCI = VRAmean+1.96*(VRAstd./(VRAn.^.5));

VRAdata = [VRAmean' VRAlowerCI' VRAupperCI']


%% %Setup mixed effects models to include effects of multiple predictors

% predictors can be within subject variables

Tbl=table(EdiPrearSubj_N,Delta*100, ...
    'VariableNames',{'Subjects','VRA'});
Tbl.Subjects = nominal(Tbl.Subjects);
%         Tbl.Edi = ordinal(Tbl.Edi);
%         Tbl.Vol = ordinal(Tbl.Vol);
%         Tbl.Pes = ordinal(Tbl.Pes);
lme = fitlme(Tbl,['VRA ~ 1 + (1|Subjects) '])
%lme = fitlme(Tbl,['VRA ~ -1 + (1|Subjects) '])
[beta] = fixedEffects(lme);
[~,~,STATS] = randomEffects(lme); % Compute the random-effects statistics (STATS)
STATS.Level = nominal(STATS.Level);
% Q: are values different between subjects? Y

%lme = fitlme(Tbl,['Edi ~ 1 + Pes + Vol + (Vol | Subjects)'])
%lme = fitlme(Tbl,['Edi ~ 1 + Pes + Vol + (Pes | Subjects) + (Vol | Subjects) '])
%lme = fitlme(Tbl,['Vol ~ 1 + Edi + Pes + (Edi| Subjects)'])

lme.Rsquared

temp = STATS.Estimate+beta;


%%


Y=(PreVEsleepArray_Edi_N(PreVEsleepArray_Subj_N==i,end)-1);

M_ = sum(~Exclude);
J = 10;
p1=PreVEsleepArray_VE_N(:,end-1)-1;
p2=PreVEsleepArray_VE_N(:,end-2)-1;
p3=PreVEsleepArray_VE_N(:,end-3)-1;
p4=PreVEsleepArray_VE_N(:,end-4)-1;
p5=PreVEsleepArray_VE_N(:,end-5)-1;
p6=PreVEsleepArray_VE_N(:,end-6)-1;
p7=PreVEsleepArray_VE_N(:,end-7)-1;
p8=PreVEsleepArray_VE_N(:,end-8)-1;
p9=PreVEsleepArray_VE_N(:,end-9)-1;
p10=PreVEsleepArray_VE_N(:,end-10)-1;


%% MLR using VE and Edi history to predict VRA on breath 1, 2, 3 etc.

clear A Rsq
for i=1:4
    predforward=i;
    y = (PreVEsleepArray_Edi_N(:,end));
    X = (PreVEsleepArray_VE_N(:,[end-5-predforward:end-predforward]));
    temp = PreVEsleepArray_VE_N(:,end);
    if 1 %0: no intercept
        X = [X 1+0*y]; %add intercept, column of ones
    end
    if 1 %incl Edi
        X = [X (PreVEsleepArray_Edi_N(:,[end-5-predforward:end-predforward]))]; %add intercept, column of ones
    end  
    if 1 %only predict drive in breaths with VE>eupnea: better Rsq
        X(temp<1,:)=[];
        y(temp<1)=[];
    end
A = X\y;
[b,bint,r,rint,stats] = regress(y,X);
Rsq(i)=stats(1);
end


clear A Rsq
for i=1:M_
    predforward=1;
    y = PreVEsleepArray_Edi_N(PreVEsleepArray_Subj_N==i,end);
    X = PreVEsleepArray_VE_N(PreVEsleepArray_Subj_N==i,[end-5-predforward:end-predforward]);
    temp = PreVEsleepArray_VE_N(PreVEsleepArray_Subj_N==i,end);
    if 1 %0: no intercept
        X = [X 1+0*y]; %add intercept, column of ones
    end
    if 1 %incl Edi
        X = [X (PreVEsleepArray_Edi_N(PreVEsleepArray_Subj_N==i,[end-5-predforward:end-predforward]))]; %add intercept, column of ones
    end
    if 1 %only predict drive in breaths with VE>eupnea: better Rsq
        X(temp<1,:)=[];
        y(temp<1)=[];
    end
A(:,i) = X\y;
[b,bint,r,rint,stats] = regress(y,X);
Rsq(i)=stats(1);
end
nansum(Rsq)/sum(~isnan(Rsq))


%% Setup mixed effects models to predict drive during sleep [VE->Edi]

Y=(PreVEsleepArray_Edi_N(:,end)-1);

M_ = sum(~Exclude);
J = 10;
p1=PreVEsleepArray_VE_N(:,end-1)-1;
p2=PreVEsleepArray_VE_N(:,end-2)-1;
p3=PreVEsleepArray_VE_N(:,end-3)-1;
p4=PreVEsleepArray_VE_N(:,end-4)-1;
p5=PreVEsleepArray_VE_N(:,end-5)-1;
p6=PreVEsleepArray_VE_N(:,end-6)-1;
p7=PreVEsleepArray_VE_N(:,end-7)-1;
p8=PreVEsleepArray_VE_N(:,end-8)-1;
p9=PreVEsleepArray_VE_N(:,end-9)-1;
p10=PreVEsleepArray_VE_N(:,end-10)-1;

%p9=1-PreVEsleepArray_VE_N(:,end);

%Y=100*Delta./p5; %dramatically lowers goodness of fit if removing this variance: success
% predictors can be within subject variables

Tbl=table(PreVEsleepArray_Subj_N,Y, ...
    'VariableNames',{'Subjects','Vdrive'});
Tbl.Subjects = nominal(Tbl.Subjects);

p=[];
for i=1:J
    p=eval(['[p p' num2str(i) ']']);
end
if 0
    i=J+1;
    p=eval(['[p p' num2str(i) ']']);
end

B = array2table(p);
Tbl = horzcat(Tbl,B);

fitstr = 'Vdrive ~ -1 ';
%fitstr = 'Vdrive ~ 1 ';
if 0
    if J==1
        fitstr = [fitstr ' + p'];
    else
    for i=1:J
        fitstr = [fitstr ' + p' num2str(i)];
    end
    end
elseif 0
    i=4
    fitstr = [fitstr ' + p' num2str(i)];    
end
if 0
    i=J+1;
    fitstr = [fitstr ' + p' num2str(i)];
end
%random effects must come after
if 0
fitstr = [fitstr ' + (1|Subjects)']; %subject-based intercept
end
if 1
    if J==1
        fitstr = [fitstr ' + (p-1|Subjects)'];
    else
    for i=1:J
        fitstr = [fitstr ' + (p' num2str(i) '-1|Subjects)'];
    end
    end
elseif 0
    i=4
    fitstr = [fitstr ' + (p' num2str(i) '-1|Subjects)'];
end
if 0
    i=J+1;
    fitstr = [fitstr ' + (p' num2str(i) '-1|Subjects)'];
end

lme = fitlme(Tbl,fitstr) %(p1-1|Subjects) + (p2-1|Subjects)

[beta] = fixedEffects(lme);
[~,~,STATS] = randomEffects(lme); % Compute the random-effects statistics (STATS)
STATS.Level = nominal(STATS.Level);

randomeffectp = STATS.pValue;
K = length(randomeffectp)/M_;
Fsig=[];
for i=1:K
    Fsig(i)=sum(randomeffectp(1+(i-1)*M_:i*M_)<0.05)/M_;
end
Fsig

lme.Rsquared

%temp = STATS.Estimate+beta;

%% %Setup mixed effects models to include effects of multiple predictors

Y=Delta*100;

M_ = sum(~Exclude)
J = 12;
p1 = PrearArray_Pbeta_N(:,end) - PrearArray_Pbeta_N(:,end-1);
p1(p1<0)=0;

%pred1(pred1>2)=2;
%pred2 = PrearArray_Pbeta_N(:,end) - PrearArray_Pbeta_N(:,end-2);
p2 = PrearArray_Palpha_N(:,end) - PrearArray_Palpha_N(:,end-1);
%p3 = PrearArray_Pdelta_N(:,end) - PrearArray_Pdelta_N(:,end-1);
p3 = PrearArray_Ptheta_N(:,end) - PrearArray_Ptheta_N(:,end-1);

p4 = 100-PrearArray_SpO2_N(:,end);

p13=144.3*p1+36.974*p2-27.815*p3; %simple model with (100-SpO2 included as a covariate)

p5=PrearArray_VE_N(:,end-1);
p6=PrearArray_VE_N(:,end-2);
p7=PrearArray_VE_N(:,end-3);
p8=PrearArray_VE_N(:,end-4);
p9=PrearArray_VE_N(:,end-5);
p10=PrearArray_VE_N(:,end-6);
p11=PrearArray_VE_N(:,end-7);
p12=PrearArray_VE_N(:,end-8);

if 1
    for i=[1:3 13]
        eval(['p' num2str(i) '=p' num2str(i) '/median(p' num2str(i) ');']);
    end
end
% pred1 = PrearArray_Pbeta_N(:,end) - PrearArray_Pbeta_N(:,end-1);

%p4(p4<0)=0;

p13(p13<0.01)=0.01;

%Y=100*Delta./p5; %dramatically lowers goodness of fit if removing this variance: success
% predictors can be within subject variables

Tbl=table(EdiPrearSubj_N,Y, ...
    'VariableNames',{'Subjects','VRA'});
Tbl.Subjects = nominal(Tbl.Subjects);

p=[];
for i=1:J
    p=eval(['[p p' num2str(i) ']']);
end
if 0
    i=J+1;
    p=eval(['[p p' num2str(i) ']']);
end

B = array2table(p);
Tbl = horzcat(Tbl,B);

fitstr = 'VRA ~ -1 ';
if 1
    for i=1:J
        fitstr = [fitstr ' + p' num2str(i)];
    end
elseif 1
    i=4
    fitstr = [fitstr ' + p' num2str(i)];    
end
if 0
    i=J+1;
    fitstr = [fitstr ' + p' num2str(i)];
end
%random effects must come after
if 0
fitstr = [fitstr ' + (1|Subjects)']; %subject-based intercept
end
if 0
    for i=1:J
        fitstr = [fitstr ' + (p' num2str(i) '-1|Subjects)'];
    end
elseif 1
    i=4
    fitstr = [fitstr ' + (p' num2str(i) '-1|Subjects)'];
end
if 0
    i=J+1;
    fitstr = [fitstr ' + (p' num2str(i) '-1|Subjects)'];
end

lme = fitlme(Tbl,fitstr) %(p1-1|Subjects) + (p2-1|Subjects)

[beta] = fixedEffects(lme);
[~,~,STATS] = randomEffects(lme); % Compute the random-effects statistics (STATS)
STATS.Level = nominal(STATS.Level);

randomeffectp = STATS.pValue;
K = length(randomeffectp)/M_;
Fsig=[];
for i=1:K
    Fsig(i)=sum(randomeffectp(1+(i-1)*M_:i*M_)<0.05)/M_;
end
Fsig

lme.Rsquared

%temp = STATS.Estimate+beta;
%%
for n=1:M
    try
    figure(1)
    plotregressionwithSEM(Delta(EdiPrearSubj_N==n),pred2(EdiPrearSubj_N==n))
    hold('off')
    pause(0.5)
    catch me
    end
end

%% Test validation that model VRA represents actual VRA
figure()
I=find(~isnan(VRAdata(:,1)));
plotregressionwithSEM(temp,VRA1(I))

%%

if 1
%made up perfect linear within subject predictor
predictor = 0*EdiPrearSubj_N;
k=rand(M,1);
k2=rand(M,1);
for i=1:M
    rangei = find(EdiPrearSubj_N==i);
    temp_ = 0.0*rand(length(rangei),1);
    temp_1 = 1/k(i)*(Delta(rangei)-k2(i));
    temp = temp_ + temp_1;
    predictor(rangei) = temp;
end

Tbl=table(EdiPrearSubj_N,Delta,predictor, ...
    'VariableNames',{'Subjects','VRA','predictor'});
Tbl.Subjects = nominal(Tbl.Subjects);
%lme = fitlme(Tbl,['VRA ~ -1 + (predictor-1|Subjects) ']) %y=mx
lme = fitlme(Tbl,['VRA ~ -1 + (predictor|Subjects)']); %y=mx+c
lme.Rsquared.Ordinary

REM = 0*Delta; REM(1:2:end)=1;
Delta2 = Delta.*(1+0.2*REM);
guesskREM = [0:0.1:1];
clear F
for i=1:length(guesskREM)
    
Delta3 = Delta2./(1+guesskREM(i).*REM);
Tbl=table(EdiPrearSubj_N,Delta3,predictor,REM, ...
    'VariableNames',{'Subjects','VRA','predictor','REM'});
Tbl.Subjects = nominal(Tbl.Subjects);
%lme = fitlme(Tbl,['VRA ~ -1 + (predictor-1|Subjects) ']) %y=mx
lme = fitlme(Tbl,['VRA ~ -1 + (predictor|Subjects)']); %y=mx+c
%same as
%lme = fitlme(Tbl,['VRA ~ -1 + (1|Subjects) + (predictor-1|Subjects) ']) %y=mx+c

F(i)=lme.Rsquared.Ordinary;
end


[beta] = fixedEffects(lme);
[~,~,STATS] = randomEffects(lme); % Compute the random-effects statistics (STATS)
STATS.Level = nominal(STATS.Level);
    %k2-beta
end
    
    
    %%
    
    
    
    
    

