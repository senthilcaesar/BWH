function [Rvalue,Pvalue,Pvalue2,Xmodel]=plotregressionwithSEM_powerlaw(xdata,ydata)

removenan=isnan(xdata)|isnan(ydata); xdata(removenan)=[]; ydata(removenan)=[]; 

%can not have negative values in power law function
xdata(xdata<0)=0;

ploton=1;
newfigure=0;
plotfill=1;
plotlines=0;
%xdata=Vpassive;
%ydata=VpassivePSG;
%regression with confidence intervals

lsqoptions=optimset('display','off');
fun = @(x,xdata)(x(1)*xdata.^x(3)+x(2));
lower=[-Inf -Inf 0.5];
upper=[Inf Inf 2];
x0=[0 0 1];
[Xmodel,Fres,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper,lsqoptions);

Ftot=sum((ydata-mean(ydata)).^2);
Rvalue=(1-Fres/Ftot)^0.5;

%get pvalue (adapted code from matlab's "NonLinearModel" function)
ssr = max(Ftot - Fres,0);
dfr = length(x0) - 1;
dfe = length(xdata) - 1 - dfr;
f = (ssr./dfr)/(Fres/dfe);
Pvalue = fcdf(1./f,dfe,dfr); % upper tail

%get pvalue versus linear model
fun2 = @(x,xdata)(x(1)*xdata+x(2));
lower=[-Inf -Inf];
upper=[Inf Inf];
x02=polyfit(xdata,ydata,1);
[Xmodel2,Fres2,~,~,~,~,~]=lsqcurvefit(fun2,x02,xdata,ydata,lower,upper,lsqoptions);
ssr = max(Fres2 - Fres,0);
dfr = length(x0) - length(x02);
dfe = length(xdata) - length(x0);
f = (ssr./dfr)/(Fres/dfe);
Pvalue2 = fcdf(1./f,dfe,dfr); % upper tail

Nx=101;
dX=(max(xdata)-min(xdata))/Nx;
xdataline = (min(xdata):dX:max(xdata))';

predopt1='observation'; %'observation' 'curve'
[ydatamodel,delta] = nlpredci(@(Xmodel,xdataline) fun(Xmodel,xdataline),xdataline,Xmodel,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);
upperline=ydatamodel+delta/1.96; %SEM or SD
lowerline=ydatamodel-delta/1.96; %SEM or SD
predopt1='curve'; %'observation' 'curve'
[~,delta] = nlpredci(@(Xmodel,xdataline) fun(Xmodel,xdataline),xdataline,Xmodel,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);
upperline2=ydatamodel+delta/1.96; %SEM or SD
lowerline2=ydatamodel-delta/1.96; %SEM or SD



if ploton
    if newfigure
        figure(88); 
    end
    plot(xdataline,ydatamodel,'color',[1 1 1]); %dummy plot
    hold('on');
    if plotfill
        fill([xdataline;flipud(xdataline)],[lowerline;flipud(upperline)],[0.95 0.95 0.9875],'linestyle','none');
        fill([xdataline;flipud(xdataline)],[lowerline2;flipud(upperline2)],[0.8 0.8 0.95],'linestyle','none');
    end
    if plotlines
        plot(xdataline,upperline,'-','color',[0.5 0.5 0.5]); 
        plot(xdataline,lowerline,'-','color',[0.5 0.5 0.5]); 
    end
    plot(xdataline,ydatamodel,'k-');
    plot(xdata,ydata,'k.','markersize',16); box('off'); 
    
end