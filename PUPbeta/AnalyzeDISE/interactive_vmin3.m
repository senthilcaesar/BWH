function interactive_vmin3(time, vol, flow, snore, snoreDB, iv1,iv2,n,l,plotSnore)
scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[10 70 scrsz(3)-20 scrsz(4)-200]);
nfig=0;
totalfig=2;
iv1 = iv1(~isnan(iv1));
iv2 = iv2(~isnan(iv2));

%% end condition
if n+l > length(vol)
    l = length(vol)-n;
end

%%
if plotSnore
    ax1=subplot(4,1,1);
    plot(time, snore), hold on
    plot([time(iv1) time(iv1)],[min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1],'r')
    plot([time(iv2) time(iv2)],[min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
    ylabel('Snore')
    xlim([time(n) time(n+l)])
    ylim([min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1])

    ax2=subplot(4,1,2);
    plot(time, snoreDB), hold on
    plot([time(iv1) time(iv1)],[min(snoreDB(n:n+l))-0.1 max(snoreDB(n:n+l))+.1],'r')
    plot([time(iv2) time(iv2)],[min(snoreDB(n:n+l))-0.1 max(snoreDB(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
    ylabel('SnoreDB')
    xlim([time(n) time(n+l)])
    ylim([min(snoreDB(n:n+l)) max(snoreDB(n:n+l))+2])

    nfig=2;
    totalfig=4;
end

ax3=subplot(totalfig,1,nfig+1);
plot(time, flow), hold on
plot([time(iv1) time(iv1)],[min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1],'r')
plot([time(iv2) time(iv2)],[min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
ylabel('Flow')
xlim([time(n) time(n+l)])
ylim([min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1])

ax4=subplot(totalfig,1,nfig+2);
plot(time, vol), hold on
plot([time(iv1) time(iv1)],[min(vol(n:n+l))-0.1 max(vol(n:n+l))+.1],'r')
plot([time(iv2) time(iv2)],[min(vol(n:n+l))-0.1 max(vol(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
ylabel('Vol')
xlabel('Time')
xlim([time(n) time(n+l)])
ylim([min(vol(n:n+l))-0.5 max(vol(n:n+l))+.5])

if plotSnore
    linkaxes([ax1 ax2 ax3 ax4],'x');
else
    linkaxes([ax3 ax4],'x');
end
