function [percentageclippedInsp, percentageclippedExp, NumberClipped, Flowclipping,RobustFlowclipping, ClipFractionBreath] = ClipDetection(ClipdetectionFlow, BB_i_start, BB_i_end,Threshold)
% Calculate range
ClipdetectionFlow1 = ClipdetectionFlow-mean(ClipdetectionFlow);
range = zeros(length(BB_i_start),1);
for j = 1:length(BB_i_start)
    range(j) = abs(max(ClipdetectionFlow1(BB_i_start(j):BB_i_end(j))))+abs(min(ClipdetectionFlow1(BB_i_start(j):BB_i_end(j))))-abs(mean(ClipdetectionFlow(BB_i_start(j):BB_i_end(j))));
end

% [a_insp,b_insp] = findpeaks((ClipdetectionFlow-mean(ClipdetectionFlow)),'MinPeakDistance',150,'MinPeakProminence',0.20,'Annotate','extents');
% [a_exp,b_exp] = findpeaks(-(ClipdetectionFlow-mean(ClipdetectionFlow)),'MinPeakDistance',150,'MinPeakProminence',0.20,'Annotate','extents');
% 
% b_tot = b_insp + b_exp;
% a_tot = abs(a_insp) + abs(a_exp);
% max(a_tot)

% Calculate average peak flow
meanPeakflowrange = mean(range);
maxPeakflowrange = max(range);
maxWindow = max(ClipdetectionFlow);
minWindow = min(ClipdetectionFlow);

ClipdetectionFlow2 = ClipdetectionFlow+(0.5*(abs(maxWindow)-minWindow)-maxWindow); %make flow symmetrical on range (clipping insp at same point as clipping exp)
maxWindow2 = max(ClipdetectionFlow2);

NumberClippedInsp = 0;
NumberClippedExp = 0;
FlowclippingInsp = zeros(length(ClipdetectionFlow2),1); %0-1 function 0 if no clipping, 1 if clipping
FlowclippingExp = zeros(length(ClipdetectionFlow2),1); %0-1 function 0 if no clipping, -1 if clipping
% For entire window

FlowclippingInsp(ClipdetectionFlow2>Threshold*maxWindow2) = 1;
NumberClippedInsp = sum(FlowclippingInsp);
FlowclippingExp(ClipdetectionFlow2<-Threshold*maxWindow2) = -1;
NumberClippedExp = sum(abs(FlowclippingExp));


NumberClipped = [NumberClippedInsp; NumberClippedExp];
Flowclipping = FlowclippingInsp + FlowclippingExp;

percentageclippedInsp = NumberClippedInsp/length(ClipdetectionFlow2>=0)*100;%bekijken wss niet 100% juist om te delen door length(ClipdetectionFlow)
percentageclippedExp = NumberClippedExp/length(ClipdetectionFlow2<=0)*100;

% Robuster version -> remove clipping if abs(flow) = min
[minClipping,~]= min(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)));
[maxClipping,~] = max(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)));
Thresholdremoval = minClipping + 0.20*(maxClipping-minClipping);
RobustFlowclipping = Flowclipping;
RobustFlowclipping(abs(ClipdetectionFlow2)<Thresholdremoval,:)=0;

if 1
    figure(999); plot(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)))
    hold on; figure; plot(abs(ClipdetectionFlow2(abs(RobustFlowclipping)==1,:)))
    min(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)))
    max(ClipdetectionFlow2(abs(Flowclipping)==1,:))
    min(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)))
end

% % Very robust version -> only retain clipping if abs(flow) = max
% [a,~]=max(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)));
% RobustFlowclipping = Flowclipping;
% RobustFlowclipping(abs(ClipdetectionFlow2)~=a,:)=0;
% if 0
%     figure; plot(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)))
%     min(ClipdetectionFlow2(abs(Flowclipping)==1,:))
%     max(ClipdetectionFlow2(abs(Flowclipping)==1,:))
%     min(abs(ClipdetectionFlow2(abs(Flowclipping)==1,:)))
% end

%% per breath (very robust version)
FlowclippingInspBreath = zeros(length(BB_i_start),1);
FlowclippingExpBreath = zeros(length(BB_i_start),1);
ClipFractionBreath = zeros(length(range),2);
for k = 1:length(range)
    NumberClippedInspBreath = 0;
    NumberClippedExpBreath = 0;
    for l = BB_i_start(k):BB_i_end(k)
        if RobustFlowclipping(l) == 1
            NumberClippedInspBreath = NumberClippedInspBreath + 1; %clipped datapoints per breath
            FlowclippingInspBreath(k) = 1; % 1 if breath clips during insp
        elseif RobustFlowclipping(l) == -1
            NumberClippedExpBreath = NumberClippedExpBreath+1;
            FlowclippingExpBreath(k) = -1; %-1 if breath clips during exp
        end
%         if ClipdetectionFlow2(l) > Threshold*maxWindow2
%             NumberClippedInspBreath = NumberClippedInspBreath+1;
%             FlowclippingInspBreath(l) = 1;
%         elseif ClipdetectionFlow2(l) < -Threshold*maxWindow2
%             NumberClippedExpBreath = NumberClippedExpBreath+1;
%             FlowclippingExpBreath(l) = -1;
%         end
        ClipFractionBreath(k,:)= [NumberClippedInspBreath NumberClippedExpBreath]./(BB_i_end(k)-BB_i_start(k))*100;
    end
end

%% per breath (normal version)
FlowclippingInspBreath = zeros(length(BB_i_start),1);
FlowclippingExpBreath = zeros(length(BB_i_start),1);
ClipFractionBreath = zeros(length(range),2);
for k = 1:length(range)
    NumberClippedInspBreath = 0;
    NumberClippedExpBreath = 0;
    for l = BB_i_start(k):BB_i_end(k)
        if Flowclipping(l) == 1
            NumberClippedInspBreath = NumberClippedInspBreath + 1; %clipped datapoints per breath
            FlowclippingInspBreath(k) = 1; % 1 if breath clips during insp
        elseif Flowclipping(l) == -1
            NumberClippedExpBreath = NumberClippedExpBreath+1;
            FlowclippingExpBreath(k) = -1; %-1 if breath clips during exp
        end
%         if ClipdetectionFlow2(l) > Threshold*maxWindow2
%             NumberClippedInspBreath = NumberClippedInspBreath+1;
%             FlowclippingInspBreath(l) = 1;
%         elseif ClipdetectionFlow2(l) < -Threshold*maxWindow2
%             NumberClippedExpBreath = NumberClippedExpBreath+1;
%             FlowclippingExpBreath(l) = -1;
%         end
        ClipFractionBreath(k,:)= [NumberClippedInspBreath NumberClippedExpBreath]./(BB_i_end(k)-BB_i_start(k))*100;
    end
end


end