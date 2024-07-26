function dataOut = BeeSwarm(X)

data = X;

N=size(data,1);

bw=0.05/(2^log10(N/100));
trans = 0.7 - 0.2* log10(N/100);
    trans(trans>1)=1;
    trans(trans<0.1)=0.1;
mn = 1-(1./(N.^0.5));

catIdx = ones(length(data),1);
figure(901010);
[~,dataOut,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
close(901010);
% data = Tout.FlowDrive(Tout.Hypopnea==0);
% catIdx = ones(length(data),1);
% [~,dataOut2,~] = plotSpread(data,'categoryIdx',catIdx,'categoryMarkers',{'.'},'categoryColors',{'r'},'binWidth',bw,'magicNumber',mn);
% close(1)
% 
% figure(1); clf(1); set(gcf,'color',[1 1 1]);
% scatter(dataOut(:,1),dataOut(:,2),10,[0.9 0.2 0.05],'filled','markerfacealpha',trans);
% hold on
% scatter(dataOut2(:,1)+1,dataOut2(:,2),10,fliplr([0.9 0.2 0.05]),'filled','markerfacealpha',trans);