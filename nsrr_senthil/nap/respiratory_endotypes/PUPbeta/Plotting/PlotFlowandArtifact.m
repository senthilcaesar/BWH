% values for pt 929
FlowArt = [0	8091;
8897	8936;
9078	9200;
13219	14134;
20935	21074;	
25211	25623;
32201	99999];
PnasalArt = [0	8089;
9081	9196;
13269	13978;
14912	21073;
25186	25586;
32239	33410];

% index to time conversion, etc.
FlowArt = FlowArt.*125;
FlowArt(1:1) = 1;
FlowArt(end:end) = length(Flow);
PnasalArt = PnasalArt.*125;
PnasalArt(1:1) = 1;
PnasalArt(end:end) = length(Flow);

% Time flow and pnasal for current pt
Time = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1));
Flow = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1));
Pnasal = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Pnasal')==1));

% do plot
figure(10); clf(figure(10));
ax(1) = subplot(6,1,1);
plot([Time(FlowArt(:,1)), Time(FlowArt(:,2))],  [0, 0], 'k-'); 
ylabel('F-art');
ax(2) = subplot(6,1,[2:3]);
plot(Time, Flow);
ylabel('Flow');
ax(3) = subplot(6,1,4);
plot([Time(PnasalArt(:,1)), Time(PnasalArt(:,2))],  [0, 0], 'k-'); 
ylabel('P-art');
ax(4) = subplot(6,1,[5:6]);
plot(Time, Pnasal);
ylabel('Pnasal');
linkaxes(ax);
