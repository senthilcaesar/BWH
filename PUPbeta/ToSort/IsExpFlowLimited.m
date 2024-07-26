function [EndVolInspExpRatio]=IsExpFlowLimited(Flow2,BrStart,BrEnd,Fs)
Flow=Flow2(BrStart:BrEnd);
t=(0:length(Flow)-1)'/Fs;
ZeroFlowShift=[-0.1 -0.05 0 0.05 0.1];
endVolRatio=NaN(size(ZeroFlowShift));
for kk=1:length(ZeroFlowShift)
    p(4)=Flow(1)-ZeroFlowShift(kk);
    Urange=(prctile(Flow,95)-p(4));
    Lrange=(p(4)-prctile(Flow,5));
    if Urange<=0 || Lrange <=0
        endVolRatio(kk)=NaN;
    else  
        p(1)=p(4)-0.25*Lrange;
        p(2)=p(4)-0.5*Lrange;
        p(3)=p(4)-0.75*Lrange;
        p(5)=p(4)+0.25*Urange;
        p(6)=p(4)+0.5*Urange;
        p(7)=p(4)+0.75*Urange;
        for ii=1:length(p)
            vol(ii,:) = cumtrapz(t,Flow-p(ii));
            endpoint(kk,ii)=vol(ii,end);
        end
        endVolRatio(kk)=sum(endpoint(kk,4)-endpoint(kk,5:7))/sum(endpoint(kk,1:3)-endpoint(kk,4));
        % endVolRatio(kk)=(endpoint(kk,4)-endpoint(kk,7))/(endpoint(kk,1)-endpoint(kk,4));
    end 
end
%     Coeff1 = polyfit(ZeroFlowShift(endVolRatio>=0.2),endVolRatio(endVolRatio>=0.2),1);
EndVolInspExpRatio=endVolRatio(3);
%     ExpLimitedSlope=Coeff1(1);
end