function h = plotellipsefill(x,y,xd,yd,C,alpha)
% 
% x = 3;
% y = 4;
% xd = 1;
% yd = 0.5;

X = x-xd:xd/100:x+xd;

Yupper = ((1 - ((X - x).^2/(xd^2))).*(yd^2)).^0.5 + y
Ylower = -((1 - ((X - x).^2/(xd^2))).*(yd^2)).^0.5 + y

h = fill([X fliplr(X)],[Yupper fliplr(Ylower)],C,'EdgeAlpha',0,'facealpha',alpha);
