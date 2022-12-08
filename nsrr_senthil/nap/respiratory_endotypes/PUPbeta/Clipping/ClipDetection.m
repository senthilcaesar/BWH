function [FclippedIB,FclippedEB,FclippedI,FclippedE,Flowclipping] = ClipDetection(ClipdetectionFlow,BB_i_start,BB_i_mid,BB_i_end,Threshold,NormalFclipped99,ploton)

if ~(exist('Threshold')==1) || isempty(Threshold)
    Threshold=0.90;
end

if ~(exist('NormalFclipped99')==1) || isempty(NormalFclipped99)
    NormalFclipped99=0.002;
end

if ~(exist('ploton')==1) || isempty(ploton)
    ploton=1;
end


maxWindow = max(ClipdetectionFlow);
minWindow = min(ClipdetectionFlow);
ClipdetectionFlow2 = ClipdetectionFlow-(0.5*(maxWindow-minWindow)+minWindow); %make flow symmetrical on range (clipping insp at same point as clipping exp)

maxWindow2 = max(ClipdetectionFlow2);

ThresU = Threshold*maxWindow2;
ThresL = -Threshold*maxWindow2;
FclippedInsp = sum(ClipdetectionFlow2>=ThresU)/length(ClipdetectionFlow2);
FclippedUse = FclippedInsp - NormalFclipped99;
if FclippedUse<0
    FclippedUse=0;
    ThresU = max(ClipdetectionFlow2)*1.01;
else
    ThresU = prctile(ClipdetectionFlow2,100*(1-FclippedUse));
end

FclippedExp = sum(ClipdetectionFlow2<=ThresL)/length(ClipdetectionFlow2);
FclippedUse = FclippedExp - NormalFclipped99;
if FclippedUse<0
    FclippedUse=0;
    ThresL = min(ClipdetectionFlow2)*1.01;
else
    ThresL = prctile(ClipdetectionFlow2,100*(FclippedUse));
end

FlowclippingInsp = zeros(length(ClipdetectionFlow2),1); %0-1 function 0 if no clipping, 1 if clipping
FlowclippingExp = zeros(length(ClipdetectionFlow2),1); %0-1 function 0 if no clipping, -1 if clipping
FlowclippingInsp(ClipdetectionFlow2>=ThresU) = 1;
FlowclippingExp(ClipdetectionFlow2<=ThresL) = -1;

FclippedI = sum(ClipdetectionFlow2>=ThresU)/length(ClipdetectionFlow2);
FclippedE = sum(ClipdetectionFlow2<=ThresL)/length(ClipdetectionFlow2);

%%
if ploton
figure(999); clf(999);
set(gcf,'color',[1 1 1]);
Time = 1:length(ClipdetectionFlow);
plot(Time,ClipdetectionFlow2);
hold on
% plot(Time,[Threshold*maxWindow2 ThresU] +0*ClipdetectionFlow2);
% plot(Time,[-Threshold*maxWindow2 ThresL] +0*ClipdetectionFlow2);
plot(Time,[ThresU] +0*ClipdetectionFlow2,'color',[0.7 0.4 0.1]);
plot(Time,[ThresL] +0*ClipdetectionFlow2,'color',[0.7 0.4 0.1]);

ClipdetectionFlow3 = ClipdetectionFlow2;
ClipdetectionFlow3(FlowclippingInsp==0 & FlowclippingExp==0)=NaN;
sum(~isnan(ClipdetectionFlow3));
plot(Time,ClipdetectionFlow3,'r','linewidth',3);
set(gca,'box','off');
end

%% [FclippedIB,FclippedEB,FclippedI,FclippedE,Flowclipping]
Flowclipping = FlowclippingInsp + FlowclippingExp;

%% per breath 
FclippedIB = zeros(length(BB_i_start),1);
FclippedEB = zeros(length(BB_i_start),1);
for i = 1:length(BB_i_start)
    rangeiI = BB_i_start(i):(BB_i_mid(i)-1);
    rangeiE = BB_i_mid(i):(BB_i_end(i)-1);
    FclippedIB(i)=sum(FlowclippingInsp(rangeiI))/length(rangeiI);
    FclippedEB(i)=sum(-FlowclippingExp(rangeiE))/length(rangeiE);
end







































