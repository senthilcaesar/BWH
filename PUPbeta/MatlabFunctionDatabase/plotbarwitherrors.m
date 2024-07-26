function plotbarwitherrors(Y,X,Errors,BarColor,baseline,barwidth,hatchwidth)

%written by Scott Sands 2016-11-04, inspired by core function errorbars.m

if 0 %example
Y = [5;10];
X = [1;2];
Errors = [1;2];
BarColor=[1 0 0;0 0.5 0];
baseline=0;
barwidth=0.8;
hatchwidth = 0.1;
end

%figure();
set(gcf,'color',[1 1 1]);
for i = 1:length(Y)
h(i) = bar(X(i),Y(i));
set(h(i),'BarWidth',barwidth,'EdgeColor',BarColor(i,:),'FaceColor',BarColor(i,:),'showbaseline','off','basevalue',baseline);    % 1= touching
hold on;
end

uppernotlower=Y>baseline;
ErrorsU = Errors.*uppernotlower;
ErrorsL = Errors.*(1-uppernotlower);
for i = 1:length(Y)
ytop = Y(i) + ErrorsU(i);
ybot = Y(i) - ErrorsL(i);
xb = [X(i) X(i) NaN X(i)-hatchwidth X(i)+hatchwidth NaN X(i)-hatchwidth X(i)+hatchwidth NaN];
yb = [ytop ybot NaN ytop ytop NaN ybot ybot NaN];
if uppernotlower(i)
    xb(7:9)=NaN;
else
    xb(4:6)=NaN;
end
herror(i) = plot(xb,yb,'k'); 
end
plot([min(X)-barwidth*0.75 max(X)+barwidth*0.75],baseline*[1 1],'k')
set(gca,'box','off','tickdir','out','fontname','arial narrow');
