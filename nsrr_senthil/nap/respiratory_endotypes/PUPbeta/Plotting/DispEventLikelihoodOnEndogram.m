function DispEventLikelihoodOnEndogram(Drive,VE,PrE,plotcolorbar)
hh = gcf;
c = hot(101);
% c = hot(116); c(102:end,:)=[];
cPr = (0:0.01:1)';
colormap(c);

if ~exist('Drive') || isempty(Drive)
    Drive = [20:20:200];
    VE = [10:10:100];
    PrE = 1-[10:9:91]/100;
end

% plot(Drive,VE,'.-');
if exist('plotcolorbar') && plotcolorbar==0
    
else
    h = colorbar();
    if 1
        set(h,'Position',[0.8216 0.2127 0.0549 0.3905]);
        yticks = get(h,'ytick')
        for i=1:length(yticks)
            ystr{i}=num2str(yticks(i)*100);
        end
        set(h,'yticklabel',ystr);
    end
end

hold on

dDrive = diff(Drive);
dVE = diff(VE);
for i=1:length(Drive)
    [~,j] = min(abs(PrE(i)-cPr));
    colortemp = c(j,:);
    if i~=length(Drive)
        plot(Drive(i)+[0 dDrive(i)/2],VE(i)+[0 dVE(i)/2],'-','color',colortemp,'linewidth',2);
    end
    if i~=1
        plot(Drive(i)-[0 dDrive(i-1)/2],VE(i)-[0 dVE(i-1)/2],'-','color',colortemp,'linewidth',2);
    end
    hold on
end
%hold off
