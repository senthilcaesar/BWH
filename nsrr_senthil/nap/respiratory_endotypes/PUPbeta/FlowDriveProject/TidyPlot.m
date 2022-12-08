function TidyPlot()
% nothing special to see here
fig = gcf;
fig.Color = [1 1 1]; % set background colour to white
fig.Units = 'inches';
fig.Position = [10 1 7 6];
ax = gca;
ax.TickDir = 'out';
box off

if 0
    
    savefig(fig, 'c:\Temp\AWP_addingFtrs.fig');
    title('');
    suptitle('');
    legend on
    legend off
    box on
    box off
end
end





