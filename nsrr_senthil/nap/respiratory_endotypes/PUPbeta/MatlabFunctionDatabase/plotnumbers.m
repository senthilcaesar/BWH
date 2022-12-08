function plotnumbers(x,y,num)
xlims=get(gca,'xlim'); dx = diff(xlims)/20;
ylims=get(gca,'ylim'); dy = diff(ylims)/20;
if isempty(num)
    for i=1:length(x)
        h=text(x(i)+dx,y(i)+dy,num2str(i),'fontname','arial narrow','fontsize',9);
    end
else
    for i=1:length(x)
        h=text(x(i)+dx/5,y(i)+dy/5,num2str(num(i)),'fontname','arial narrow','fontsize',8);
    end
end