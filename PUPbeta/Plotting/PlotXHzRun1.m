
logit = @(p) log(p./(1-p));
logitinverse = @(x) 1./(1+exp(-x));

plotchanneltext = { ...
        {'Epochs','EventsResp'},{'Flow'},{'SpO2'},{'EKG'},{'HR'},{'EventsAr'},{'EventsArWS'},{'WPr','ArPr'},{'WPr','ArPr'},{'ArIntensity','ArIntensityOrig'},{'EEG'}...
        };

    plotTransformLogitList = [8];
%plotchanneltext = { ...
%        {'Epochs','EventsResp'},{'Flow'},{'SpO2'},{'HR'},{'EventsAr'},{'EventsArWSB','WPrB','ArPrB'}...
%        };

plotTransformLogit = zeros(length(plotchanneltext),1);
    plotTransformLogit(plotTransformLogitList)=1;
    
run PlotXHzData

ax2 = ax;

global xvalues yvalues range ax2 xvalues_rm yvalues_rm
xvalues=[];
yvalues=[];
tminmax = [DataEventHypnog_Mat(1,1) DataEventHypnog_Mat(end,1)];
range = diff(tminmax)*0.5;
plotwithslider(tminmax)





