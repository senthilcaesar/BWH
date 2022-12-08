function DotAcrossPlot(X)
%%

h=plot(1,1);
clear h
colors = [0 0 0];

J = nanstd(X);
N = size(X,1);
trans = 0.7 - 0.2* log10(N/100);
    trans(trans>0.5)=0.5;
    trans(trans<0.1)=0.1;

Y = 0.*X + [1:size(X,2)];
for j=1:size(X,2)
    if J(j)>0
       I=~isnan(X(:,j));
       temp = BeeSwarm(X(:,j)); 
       Y(I,j) = Y(I,j)+ temp(:,1) -1;
    end
end

for i=1:size(X,1)
    data = X(i,:);
    scatter(Y(i,:),data,80,colors,'filled','markerfacealpha',trans)
    hold on
    plot(Y(i,:),data,'k-');
end
xticks([]);

data = nanmedian(X);
for j=1:size(X,2)
    data1 = data(j);
    if data1~=0
        plot(j+[-0.4 0.4],data1+[0 0],'k-','linewidth',4);
    end
end
xlim([1-0.4 size(X,2)+0.4]);


