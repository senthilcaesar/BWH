function [Xbin,Ybin] = QuantileBin(X,Y,N,ploton)

A = [X(:) Y(:)];
A = sortrows(A,1);

if ~exist('N') || isempty(N)
    N=10;
end
if ~exist('ploton')
    ploton=0
end

centilestep = 100/N;
thres = 0:centilestep:100;
thres(end)=[];
Xs = prctile(A(:,1),thres);
Bin = sum(A(:,1)>Xs,2);

clear Xbin Ybin Yu Yl
for i=1:N %will be a faster way to do this but it works ok
    I = Bin==i;
    Xbin(i,1)=nanmedian(A(I,1));
    Ybin(i,1)=nanmedian(A(I,2));
    Yu(i,1)=prctile(A(I,2),75);
    Yl(i,1)=prctile(A(I,2),25);
end

if ploton
    hold on
    h=fill([Xbin;flipud(Xbin)],[Yu;flipud(Yl)],[0.7 0.6 0.2],'EdgeColor','None','Facealpha',0.4);
    plot(Xbin,Ybin,'k-');
end




