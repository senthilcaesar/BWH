function [h,x,t,fs,label] = plotFromEDF(EDFfilenamedir,Ch)
if 0
    Ch=T.ChannelNumberA(n);
    
end
[x,fs,~,~,label,~,~,~,~] = readedfrev3(EDFfilenamedir,Ch-1,0,Inf);
         t=(0:(1/fs):0+(length(x)-1)*(1/fs))'; % This is the time vector associated with the _XHz Flow data.
         
         h=plot(t,x);
         plotlabel = [label '|Ch=' num2str(Ch)];
         ylabel(plotlabel);
         box('off');
         set(gca,'tickdir','out')
