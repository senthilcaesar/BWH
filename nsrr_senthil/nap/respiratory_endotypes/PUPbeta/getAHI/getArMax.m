function [WPrMax,ArPrMax,WSBalanceMax,ArBalanceMax,WPrMax3,ArPrMax3,WSBalanceMax3,ArBalanceMax3]=getArMax(ArT,WPr,ArPr,dt)

logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p));


ArPr=max([WPr ArPr]')';

clear WPrMax ArPrMax WPrMax3 ArPrMax3
for i=1:length(ArT.starttimesi)
    li=ArT.starttimesi(i);
    ri=ArT.endtimesi(i)-1;
    ri((ri-li)*dt>30)=li + round(30/dt);
    WPrMax(i,1) = max(WPr(li:ri));
    ArPrMax(i,1) = max(ArPr(li:ri));
    dF = 3/ArT.EventDuration(i);
    dF(dF>1)=1;
    WPrMax3(i,1) = prctile(WPr(li:ri),(1-dF)*100);
    ArPrMax3(i,1) = prctile(ArPr(li:ri),(1-dF)*100);
end
WPrMax=WPrMax;
ArPrMax=ArPrMax;
WSBalanceMax=logit(WPrMax);
ArBalanceMax=logit(ArPrMax);
WPrMax3=WPrMax3;
ArPrMax3=ArPrMax3;
WSBalanceMax3=logit(WPrMax3);
ArBalanceMax3=logit(ArPrMax3);

