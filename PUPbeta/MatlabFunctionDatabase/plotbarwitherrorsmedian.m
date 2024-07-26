function plotbarwitherrorsmedian(Y,X,ErrorU,ErrorL,BarColor,baseline,barwidth,linewidth)

%written by Scott Sands 2016-11-04, inspired by core function errorbars.m

if 0 %example
Y = [5;8];
X = [1;2];
ErrorU = [10;12];
ErrorL = [-3 -5];
BarColor=[1 0 0;0 0.5 0];
baseline=0;
barwidth=0.8;
linewidth=0.1;
end

%figure();
set(gcf,'color',[1 1 1]);
for i = 1:length(Y)
    fill([X(i)-barwidth/2 X(i)+barwidth/2 X(i)+barwidth/2 X(i)-barwidth/2 X(i)-barwidth/2],...
        [ErrorL(i) ErrorL(i) Y(i)-linewidth/2 Y(i)-linewidth/2 ErrorL(i)],...
        BarColor(i,:),'EdgeColor','none');
    hold on;
    fill([X(i)-barwidth/2 X(i)+barwidth/2 X(i)+barwidth/2 X(i)-barwidth/2 X(i)-barwidth/2],...
        [Y(i)+linewidth/2 Y(i)+linewidth/2 ErrorU(i) ErrorU(i) Y(i)+linewidth/2],...
        BarColor(i,:),'EdgeColor','none');
%line([X(i)-barwidth/2 X(i)+barwidth/1.92],[Y(i) Y(i)],'linewidth',linewidth,'color',[1 1 1]);
hold on;
end
if ~isnan(baseline)
line([min(X)-barwidth/2*1.2 max(X)+barwidth/2*1.2],[baseline baseline],'linewidth',1,'color',[0 0 0]);
end

% uppernotlower=Y>baseline;
% ErrorsU = Errors.*uppernotlower;
% ErrorsL = Errors.*(1-uppernotlower);
% for i = 1:length(Y)
% ytop = Y(i) + ErrorsU(i);
% ybot = Y(i) - ErrorsL(i);
% xb = [X(i) X(i) NaN X(i)-hatchwidth X(i)+hatchwidth NaN X(i)-hatchwidth X(i)+hatchwidth NaN];
% yb = [ytop ybot NaN ytop ytop NaN ybot ybot NaN];
% if uppernotlower(i)
%     xb(7:9)=NaN;
% else
%     xb(4:6)=NaN;
% end
% herror(i) = plot(xb,yb,'k'); 
% end
% plot([min(X)-barwidth*0.6 max(X)+barwidth*0.8],baseline*[1 1],'k')
set(gca,'box','off','tickdir','out','fontname','arial narrow');
