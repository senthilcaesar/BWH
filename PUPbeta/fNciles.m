function [Tout,Tset,Xtile]=fNciles(Tdata,Nciles)

%Tdata =  array2table([x_ y_]);

if ~exist('Nciles')
    Nciles=10;
end

x_ = Tdata{:,1};
Xtile = zeros(height(Tdata),1);

for j=1:Nciles
    lower = (j-1)*(100/Nciles);
    upper = j*(100/Nciles);
  
    if lower<0, lower=0; end
    if lower>100, lower=100; end
    if upper>100, upper=100; end
    I = (x_>=prctile(x_,lower)&x_<=prctile(x_,upper));
    Xtile(I)=j;
end
Xtile(Xtile==0)=NaN;

clear Tset
Tmedian = [];
for j=1:Nciles
    Tmedian = [Tmedian;nanmedian(Tdata{Xtile==[j],:})];
end
Tmedian = array2table(Tmedian);
Tmedian.Properties.VariableNames = Tdata.Properties.VariableNames;
Tset.Tmedian=Tmedian;
Tmedian.Properties.VariableNames = strcat('Median_',Tmedian.Properties.VariableNames);

Tp25 = [];
for j=1:Nciles
    Tp25 = [Tp25;prctile(Tdata{Xtile==[j],:},25)];
end
Tp25 = array2table(Tp25);
Tp25.Properties.VariableNames = Tdata.Properties.VariableNames;
Tset.Tp25=Tp25;
Tp25.Properties.VariableNames = strcat('p25_',Tp25.Properties.VariableNames);

Tp75 = [];
for j=1:Nciles
    Tp75 = [Tp75;prctile(Tdata{Xtile==[j],:},75)];
end
Tp75 = array2table(Tp75);
Tp75.Properties.VariableNames = Tdata.Properties.VariableNames;
Tset.Tp75=Tp75;
Tp75.Properties.VariableNames = strcat('p75_',Tp75.Properties.VariableNames);

Tmean = [];
for j=1:Nciles
    Tmean = [Tmean;nanmean(Tdata{Xtile==[j],:})];
end
Tmean = array2table(Tmean);
Tmean.Properties.VariableNames = Tdata.Properties.VariableNames;
Tset.Tmean=Tmean;
Tmean.Properties.VariableNames = strcat('Mean_',Tmean.Properties.VariableNames);

Tout = [Tmedian,Tp25,Tp75,Tmean];
