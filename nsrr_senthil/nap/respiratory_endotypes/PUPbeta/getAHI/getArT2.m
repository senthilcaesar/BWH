function [ArIInfo,ArT]=getArT2(ArT,Time,ArIntOr,ArInt,ArIInfo)

dt = Time(2)-Time(1);
maxI = round(15/dt);

clear Tnew
for i=1:height(ArT)
    I = find(Time>=ArT.EventStart(i) & Time<=(ArT.EventEnd(i)));
    I(maxI+1:end)=[];
    Tnew.ArIntOr(i,1) = nanmax(ArIntOr(I));
    Tnew.ArInt(i,1) = nanmax(ArInt(I));
%     Tnew.WPr(i,1) = nanmax(WPr(I));
%     Tnew.ArPr(i,1) = nanmax(ArPr(I));
%     Tnew.WArPr(i,1) = nanmax(WArPr(I));   
end
Tnew = struct2table(Tnew);
ArT = [ArT Tnew];

% Evts.ArT.WPrScore = logit(Evts.ArT.WPr);
%     Evts.ArT.WPrScore(Evts.ArT.WPrScore>10)=10;
%     Evts.ArT.WPrScore(Evts.ArT.WPrScore<-10)=-10;
% Evts.ArT.ArPrScore = logit(Evts.ArT.ArPr);
%     Evts.ArT.ArPrScore(Evts.ArT.ArPrScore>10)=10;
%     Evts.ArT.ArPrScore(Evts.ArT.ArPrScore<-10)=-10;
% Evts.ArT.WArPrScore = logit(Evts.ArT.WArPr);
%     Evts.ArT.WArPrScore(Evts.ArT.WArPrScore>10)=10;
%     Evts.ArT.WArPrScore(Evts.ArT.WArPrScore<-10)=-10;

    
    
%Keep (More) Things
Icrit = ArT.AASMarousal==1;
    
varlist = ArT.Properties.VariableNames(10:end);
Tmedians = array2table(nanmedian(ArT{Icrit,varlist}));
Tmedians.Properties.VariableNames = varlist;
for i=1:width(Tmedians)
    Tmedians.Properties.VariableNames{i} = ['ArMed_' Tmedians.Properties.VariableNames{i}];
end
Tmedians = table2struct(Tmedians);

temp = ArT.EventDuration(Icrit);
    temp(temp>30)=30;
    Tmedians.ArMed_Duration = nanmedian(temp);
    
 f = fieldnames(Tmedians);
 for i = 1:length(f)
   ArIInfo.(f{i}) = Tmedians.(f{i});
 end
 
% disp(ArIInfo);