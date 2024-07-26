%% curve fits
function addCurveFits(InspTime, InspFlow, InspCurveColor, ExpTime, ExpFlow, ExpCurveColor)
Ti25 = round(0.25*(length(InspTime)),0); Ti75 = round(0.75*(length(InspTime)),0);   

% sine fit to expiration
threeX = [ExpTime(1); ExpTime(round(length(ExpTime)/2,0)); ExpTime(end)];
threeY = [ExpFlow(1); min(ExpFlow); ExpFlow(end)];
[pp,~,mu] = polyfit(threeX, threeY, 2);
y1 = polyval(pp, ExpTime,[],mu);

% sine fit to middle 50% of inspiration
threeX = [InspTime(1); InspTime(round(length(InspTime)/2,0)); InspTime(end)];
threeY = [InspFlow(1); max(InspFlow); InspFlow(end)];
[pp,~,mu] = polyfit(threeX, threeY, 2);
y2 = polyval(pp, InspTime,[],mu);

% overlay on plot
plot(ExpTime, y1, 'k--','linewidth', 1.5);
%jbfill(ExpTime,  ExpFlow', y1, ExpCurveColor, ExpCurveColor, 0, 0.1);

%plot(InspTime(Ti25:Ti75), y2(Ti25:Ti75), 'k--','linewidth', 1.5);
%jbfill(InspTime(Ti25:Ti75), InspFlow(Ti25:Ti75)', y2(Ti25:Ti75), InspCurveColor, InspCurveColor, 0, 0.5);

% % parameters for next curve fits
% c=max(InspFlow)*1;
% N=length(InspFlow);
% x=(1:N)-round(N/2);   %x=(1:N);
% % inverted parabola
% y3=c*((1-(x.^2)/(N/2)^2));
% % ellipseI
% y4=c*(((1-(x.^2)/(N/2)^2)).^(1/2));
% % hypcosI, hyperbolic cosine arch 0<b<Inf
% b=10;
% y5=c*((1-cosh((x/(N/2))*log(1+b+sqrt(2+b)*sqrt(b)))+b)/b);

% %plot(InspTime, y3, 'c');
% %plot(InspTime, y4, 'm');
% plot(InspTime, y5, 'g--','linewidth', 1.5);
end
