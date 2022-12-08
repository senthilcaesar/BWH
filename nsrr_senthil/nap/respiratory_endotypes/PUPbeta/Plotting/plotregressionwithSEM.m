function [Rvalue,Pvalue,Slope,Intercept,ModelSEM]=plotregressionwithSEM(xdata,ydata,plotSEMoff,plotfill,mkcol,filcol,mksize,mkalpha)

if ~exist('plotSEMoff') || isempty(plotSEMoff)
    plotSEMoff=0;
end
if ~exist('plotfill') || isempty(plotfill)
    plotfill=2;
end
if ~exist('mkcol') || isempty(mkcol)
    mkcol = [0 0 0];
end
if ~exist('filcol') || isempty(filcol)
    filcol = [0.8 0.8 0.95];
end
if ~exist('mksize') || isempty(mksize)
    mksize = 15;
end
if ~exist('mkalpha') || isempty(mkalpha)
    mkalpha = 0.5;
end

filcol2 = 1-((1-filcol)./4);

removenan=isnan(xdata)|isnan(ydata); xdata(removenan)=[]; ydata(removenan)=[];

ploton=1;
newfigure=0;
plotlines=0;


%xdata=Vpassive;
%ydata=VpassivePSG;
%regression with confidence intervals

lsqoptions=optimset('display','off');
fun = @(x,xdata)(x(1)*xdata+x(2));
lower=[-Inf -Inf];
upper=[Inf Inf];
x0=polyfit(xdata,ydata,1);
[Xmodel,Fres,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper,lsqoptions);
Slope = Xmodel(1); Intercept = Xmodel(2);
Ftot=sum((ydata-mean(ydata)).^2);
Rvalue=(1-Fres/Ftot)^0.5;
if Slope<0
    Rvalue=-1*Rvalue;
end

%get pvalue (adapted code from matlab's "NonLinearModel" function)
ssr = max(Ftot - Fres,0);
dfr = length(x0) - 1;
dfe = length(xdata) - 1 - dfr;
f = (ssr./dfr)/(Fres/dfe);
Pvalue = fcdf(1./f,dfe,dfr); % upper tail

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

Xmodel_ci = nlparci(Xmodel,RESIDUAL,'jacobian',JACOBIAN);
ModelSEM = (Xmodel_ci(:,2)'-Xmodel)/1.96;

if ploton
    if newfigure
        figure(88);
    end
    if ~plotSEMoff
        plot(xdataline,ydatamodel,'color',[1 1 1]); %dummy plot
        hold('on');
        if plotfill>0
            if plotfill==2
                fill([xdataline;flipud(xdataline)],[lowerline;flipud(upperline)],filcol2,'linestyle','none');
            end
            fill([xdataline;flipud(xdataline)],[lowerline2;flipud(upperline2)],filcol,'linestyle','none');
            
        end
        if plotlines
            plot(xdataline,upperline,'-','color',[0.5 0.5 0.5]);
            plot(xdataline,lowerline,'-','color',[0.5 0.5 0.5]);
        end
        plot(xdataline,ydatamodel,'k-');
    end
    %plot(xdata,ydata,'k.', 'markersize',mksize, 'color',mkcol); box('off');
    scatter(xdata,ydata,mksize*2.5,mkcol,'filled','markerfacealpha',mkalpha); box('off');
end

set(gca,'tickdir','out')