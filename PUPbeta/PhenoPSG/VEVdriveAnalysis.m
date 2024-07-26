function [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,VEdecilesUpper,VEdecilesLower, VEdecilesMean] = VEVdriveAnalysis(x,y,arthres,ploton,xeup,Nciles,Nbootstrap,NcileLimits,filt121)
global SAfontsize
I = isnan(y)|isnan(x);
y(I)=[];
x(I)=[];
if ~exist('filt121')
    filt121=1;
end
%x is drive, binned data by values of x, find y

% Create summary plot of VE versus Vdrive in deciles based on all sleep data
%Nbootstrap=0; %larger is more accurate but slower, previously used 50, 0 for exclude
if ~exist('Nbootstrap','var')
    Nbootstrap=0;
end
    
disp(['ploton: ' num2str(ploton)]);

%line of identity or gridline
if ploton&&xeup==1
    plot([1 1],[0 1.6],':','color',[0.7 0.7 0.7]);
    plot([0 3],[1 1],':','color',[0.7 0.7 0.7]);
end

x_=x;
y_=y;

% Ncile choices   :  111   121   131   141   151   161   171
% Nciles smoothing:   10    15    20    25    30    35    40
%instead, set nciles and limits in start here, define smoothing form there
if Nciles>100 %131
    smoothNciles=1;
    if ~exist('NcileLimits')
        NcileLimits=[5 95]; %[]
    end
    NcileW = ((Nciles-1) - NcileLimits(2) + NcileLimits(1))/2;
else
    smoothNciles=0;
end

clear ybindata
for j=1:Nciles
    if smoothNciles
        mid = (j-1)-NcileW*((NcileW-NcileLimits(1))/NcileW);
        lower = mid - NcileW;
        upper = mid + NcileW;
    else
        lower = (j-1)*(100/Nciles);
        upper = j*(100/Nciles);
    end
    if lower<0, lower=0; end
    if lower>100, lower=100; end
    if upper>100, upper=100; end
    I = (x_>=prctile(x_,lower)&x_<=prctile(x_,upper));
    if Nbootstrap>0
        Yci=bootci(Nbootstrap,@median,y_(I)); 
    else
        Yci=[NaN;NaN];
    end
    VEdecilesUpper = prctile(y_(I),75);
    VEdecilesLower = prctile(y_(I),25);
    ybindata(j,:)=[median(y_(I)),VEdecilesLower,VEdecilesUpper,sum(I),nanmedian(x_(I)),Yci',nanmean(y_(I))];
   
end

if filt121>0
    ybindata(:,1) = filter121(ybindata(:,1),filt121);
    ybindata(:,2) = filter121(ybindata(:,2),filt121);
    ybindata(:,3) = filter121(ybindata(:,3),filt121);
end

VEdecilesUpper = ybindata(:,3);
VEdecilesLower = ybindata(:,2);

%Vpassive
if xeup<ybindata(1,5)||xeup>ybindata(end,5)
    [~,ii] = unique(ybindata(:,5));
    Vpassive=interp1(ybindata(ii,5),ybindata(ii,1),xeup,'nearest','extrap');
    VpassiveCI=[interp1(ybindata(ii,5),ybindata(ii,6),xeup,'nearest','extrap') ...
                 interp1(ybindata(ii,5),ybindata(ii,7),xeup,'nearest','extrap')];
else
    [~,ii] = unique(ybindata(:,5));
    Vpassive=interp1(ybindata(ii,5),ybindata(ii,1),xeup);
    VpassiveCI=[interp1(ybindata(ii,5),ybindata(ii,6),xeup) ...
                 interp1(ybindata(ii,5),ybindata(ii,7),xeup)];
end

%ci=bootci(length(data)*20,@median,data);
%deltaci=abs(median(data)-ci);

%Vactive

Xtemp = arthres;

if Xtemp<ybindata(1,5)||Xtemp>ybindata(end,5)
    Vactive=interp1(ybindata(ii,5),ybindata(ii,1),Xtemp,'nearest','extrap');
    VactiveCI=[interp1(ybindata(ii,5),ybindata(ii,6),Xtemp,'nearest','extrap') ...
                 interp1(ybindata(ii,5),ybindata(ii,7),Xtemp,'nearest','extrap')];
else
    Vactive=interp1(ybindata(ii,5),ybindata(ii,1),Xtemp);
    VactiveCI=[interp1(ybindata(ii,5),ybindata(ii,6),Xtemp) ...
                 interp1(ybindata(ii,5),ybindata(ii,7),Xtemp)];
end

if ploton
    hold('on');
    fill([ybindata(:,5);flipud(ybindata(:,5))],[ybindata(:,2);flipud(ybindata(:,3))],[0.9500 0.2000  0.0500],'edgecolor','none','facealpha',0.5);
    
    plot(ybindata(:,5),ybindata(:,1),'k');
    
    plot(ybindata(:,5),ybindata(:,1),'k-','linewidth',2);
    
    try
        ylim([0 max(ybindata(:,3))]);
    catch me
    end
    
    %Plot Traits
%     if xeup==1
%         plot(xeup,1,'b.','markersize',25); %Eupnea
%     end
    xplotpoint = 100;
    if ybindata(1,5)>xeup 
        xplotpoint = ybindata(1,5);
    end
    plot(xplotpoint,Vpassive,'.','markersize',25,'color',[0.6 0.3 0.9]); %Vpassive
    plot([arthres arthres],[0 max(arthres,Vactive)],'g'); %Arthres line
    plot(arthres,Vactive,'r.','markersize',25);  %Vactive
    
    set(gca,'tickdir','out','fontsize',SAfontsize);
    axis([0 3.2 0 1.7]);
    box('off');
end
VEdeciles = ybindata(:,1);
Vdrivedeciles = ybindata(:,5);
VEdecilesMean = ybindata(:,8);