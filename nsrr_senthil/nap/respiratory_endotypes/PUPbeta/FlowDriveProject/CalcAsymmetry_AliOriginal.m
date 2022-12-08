function [AsymIndex]=CalcAsymmetry(Flow,Br,Fs_in)

InspFlow=Flow(1:Br.BrMid);
ExpFlow=Flow(Br.BrMid:end);

Fs=250;
FlowA=resample(Flow,Fs,Fs_in);
InspFlow=resample(InspFlow,Fs,Fs_in);
ExpFlow=resample(ExpFlow,Fs,Fs_in);

% InspFlow=InspFlow/mean(InspFlow);
% ExpFlow=ExpFlow/mean(-ExpFlow);
% ExpFlow=ExpFlow(ExpFlow>=1);

% if length(ExpFlow)>3*Fs
%     ExpFlow=ExpFlow(1:3*Fs);
% end

% InspFlowhat=InspModel(InspFlow,ExpFlow);
% [~,~,ExpFlowhat]=findExpTau(ExpFlow,Fs);
% if length(ExpFlow)>length(InspFlowhat) 
%     sighat=[InspFlowhat;ExpFlowhat(1:length(InspFlowhat))];
%     origSig=[InspFlow;ExpFlow(1:length(InspFlow))];
% else
%     sighat=[InspFlowhat;ExpFlowhat];
%     origSig=[InspFlow;ExpFlow];
% end
% % tau_neural=0.3890;
% % tau_mech=0.3357;
% %     
% % [~,ExpFlowHat,~]=SingleBreathModelExp([tau_neural;tau_mech], ExpFlow ,Fs);
% % ExpFlow_t=-fliplr(ExpFlowHat)';
% % 
% % sighat=[ExpFlow_t; 0 ;ExpFlow];
% % 
% % tt=(0:length(sighat)-1)/Fs;
% % t_mid=tt(length(ExpFlow)+1);
% % 
% % [AsymIndex]=findAssymIdx(sighat,Fs,t_mid);

[AsymIndex]=1-findAssymIdx(ExpFlow,Fs,FlowA);
% [AsymIndexInsp]=findAssymIdx(InspFlow,Fs);
% AsymIndex=AsymIndexExp/(AsymIndexInsp+AsymIndexExp);
end

function [AssymIdx]=findAssymIdx(x,Fs,x_t,t_mid)
    x=x';
    t2=(0:length(x_t)-1)/Fs; %vector time
    t=(0:length(x)-1)/Fs; %vector time
    if nargin<=3
        if rem(length(x),2)==0       
            x=[x x(end)];
        end
        t_mid=t((length(x)+1)/2); 
    end
    t=(0:length(x)-1)/Fs; %vector time    
        
    xmt=[fliplr(x(t>=t_mid)) fliplr(x(t<t_mid))];
    xe=0.5*(xmt+x);
    xo=0.5*(x-xmt);
    AssymIdxE=sum(xe.^2)/(length(xe)*(max(abs(xe)))^2);

    AssymIdx=sum(xo.^2)/(length(xo)*(max(abs(xo)))^2);
    AssymIdxT=sum(x.^2)/(length(x)*(max(abs(x)))^2);

    mn=min(x_t);
    mx=max(x_t);
    
    mn2=min(x);
    mx2=max(x);
    
    mn3=min(xo);
    mx3=max(xo);
    
    mn4=min(xe);
    mx4=max(xe);
        
    if 0
        figure
        clf
%         ax1=subplot_tight(2,2,1,0.04);
%         plot(x_t,'k','LineWidth',2);axis tight;ylim([mn-0.1 mx+0.1])
%         
%         set(ax1,'XTick',[]);
%         set(ax1,'XTickLabel',[]);
% %         set(ax1,'YTick',[]);
% %         set(ax1,'YTickLabel',[]);
        ax1=subplot_tight(1,3,1,0.05);
        plot(x,'k','LineWidth',2);axis tight;ylim([mn2-0.05 mx3+0.05])
        hold on;plot(round(t_mid*Fs)*ones(2,1),[mn2-0.05 mx3+0.05],'--k','LineWidth',1)
        hold off
        set(ax1,'XTick',[]);
        set(ax1,'XTickLabel',[]);
        set(ax1,'YTick',[mn2 0]);
        set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',mn2)) 0 ]);
%         title(['EFL=' num2str(AssymIdxT)])
        text(round(length(x)/5),0.1,['EFL=' num2str(AssymIdxT)])
        ax1=subplot_tight(1,3,2,0.05);
        plot(xo,'k','LineWidth',2);axis tight;ylim([mn2-0.05 mx3+0.05])
        hold on;plot(round(t_mid*Fs)*ones(2,1),[mn2-0.05 mx3+0.05],'--k','LineWidth',1)
        hold off
        text(round(length(x)/5),0.1,['EFL=' num2str(AssymIdx)])
%         title(['EFL=' num2str(AssymIdx)])
        set(ax1,'XTick',[]);
        set(ax1,'XTickLabel',[]);
        set(ax1,'YTick',[mn3 0 mx3]);
        set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',mn3)) 0 str2num(sprintf('%0.2f',mx3))]);
        ax1=subplot_tight(1,3,3,0.05);
        plot(xe,'k','LineWidth',2);axis tight;ylim([mn2-0.05 mx3+0.05])
        hold on;plot(round(t_mid*Fs)*ones(2,1),[mn2-0.05 mx3+0.05],'--k','LineWidth',1)
        hold off
%         title(['EFL=' num2str(AssymIdxE)])
        text(round(length(x)/5),0.1,['EFL=' num2str(AssymIdxE)])
        
        set(ax1,'XTick',[]);
        set(ax1,'XTickLabel',[]);
        set(ax1,'YTick',[mn3 0 mx3]);
        set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',mn3)) 0 str2num(sprintf('%0.2f',mx3))]);
    end
    %     AssymIdx=AssymIdxE*AssymIdx;
%     title(['Odd part' num2str(AssymIdx)])
%     figure,
%     ax1=subplot_tight(2,4,1,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.62 0]);set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',-0.62)) 0 ]);
%     ax1=subplot_tight(2,4,2,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.62 0]);set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',-0.62)) 0 ]);
%     ax1=subplot_tight(2,4,3,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.62 0]);set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',-0.62)) 0 ]);
%     ax1=subplot_tight(2,4,4,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.62 0]);set(ax1,'YTickLabel',[str2num(sprintf('%0.2f',-0.62)) 0 ]);
%     ax1=subplot_tight(2,4,5,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.3 0 0.3]);set(ax1,'YTickLabel',[str2num(sprintf('%0.1f',-0.3)) 0 str2num(sprintf('%0.1f',0.3))]);
%     ax1=subplot_tight(2,4,6,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.3 0 0.3]);set(ax1,'YTickLabel',[str2num(sprintf('%0.1f',-0.3)) 0 str2num(sprintf('%0.1f',0.3))]);
%     ax1=subplot_tight(2,4,7,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.3 0 0.3]);set(ax1,'YTickLabel',[str2num(sprintf('%0.1f',-0.3)) 0 str2num(sprintf('%0.1f',0.3))]);
%     ax1=subplot_tight(2,4,8,0.03);set(ax1,'XTick',[]);set(ax1,'XTickLabel',[]);
%     set(ax1,'YTick',[-0.3 0 0.3]);set(ax1,'YTickLabel',[str2num(sprintf('%0.1f',-0.3)) 0 str2num(sprintf('%0.1f',0.3))]);
    
end

function [y_insp]=InspModel(InspFlow,ExpFlow)

    FlowData=InspFlow;
    
    
%     start=[max(FlowData)*2]; % A reasonable estimate of a starting physiological value for each parameter
%     lower=[0]; % Minimum conveivable physiological value for each parameter
%     upper=[max(FlowData)*20]; % Maximum conveivable physiological value for each parameter
% 
%     OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',4000,'MaxFunEvals',150,'Algorithm','interior-point');
%     Parameters=start;
%     [c,~,~,~] = fmincon(@(Parameters) TheInvertedParabola(Parameters,FlowData),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]

    c=abs(prctile(ExpFlow,0.5));
    [~,y_insp]=TheInvertedParabola(c,FlowData);
    y_insp=y_insp';
%     figure(1);subplot(211);plot(FlowData);hold on;plot(y_insp,'r');hold off
%     
%     subplot(212);plot(ExpFlow);
     
   
end
function [tau_mech,tau_neural,y]=findExpTau(ExpFlow,Fs)

    FlowData=ExpFlow;
    start=[0.2 0.2]; % A reasonable estimate of a starting physiological value for each parameter
    lower=[0 0]; % Minimum conveivable physiological value for each parameter
    upper=[1 1]; % Maximum conveivable physiological value for each parameter

    OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',4000,'MaxFunEvals',150,'Algorithm','interior-point');
    Parameters=start;
    [FitParams] = fmincon(@(Parameters) SingleBreathModelExp(Parameters, FlowData ,Fs),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]

    [~,y]=SingleBreathModelExp(FitParams,FlowData,Fs);
    y=y';
    tau_neural=FitParams(1);
    tau_mech=FitParams(2); 
    if 0
    plot((0:length(ExpFlow)-1)/Fs,ExpFlow);
     
    title(['tau_mec=' num2str(tau_mech) ', ' 'tau_neu=' num2str(tau_neural)])
    end
end