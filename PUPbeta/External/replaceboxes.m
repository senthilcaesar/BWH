function replaceboxes
% from https://au.mathworks.com/matlabcentral/answers/5500-nice-boxes-for-great-linewidths
h = findobj('tag','Box');
n = length(h);
for k = 1:n
    x = get(h(k),'XData');
    y = get(h(k),'YData');
    c = get(h(k),'Color');
    l = get(h(k),'LineWidth');
    ht = y(2)-y(1);
    wd = x(3)-x(1);
    rectangle('position',[x(1),y(1),wd,ht],'EdgeColor',c,'LineWidth',l)
end
delete(h);