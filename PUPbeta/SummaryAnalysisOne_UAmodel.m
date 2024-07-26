function [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,Vcomp,VcompCI]=SummaryAnalysisOne_UAmodel(NwindowsForUA,BreathDataTable3,ArThresPSG1,varlist,ploton,settings)
global SAfontsize

if settings.getCIs
   Nbootstrap = settings.Nbootstrap
else
   Nbootstrap = 0;
end

%NwindowsForUA=length(DataOut);
if NwindowsForUA>0
        x=BreathDataTable3.Vdr_est+1;    %later expand to BreathDataTable3.(DriveSignal), where DriveSignal = 'VdriveModel' and BreathDataTable2.DriveModel = BreathDataTable2.Vdr_est+1;
        
        t=BreathDataTable3.Time_start;

        y = BreathDataTable3.VI;
        a = BreathDataTable3.AR3;
        %veup = BreathDataTable3.Veup;
        hyp = BreathDataTable3.hypnog_B;
        %winX = floor(BreathDataTable3.UniqueID/1000);
        %t2=BreathDataTable3.Time_end - BreathDataTable3.Time_start;
        %e=BreathDataTable3.E1;
        nota = BreathDataTable3.notAR3; %make this earlier
        

%     if 0 %original 
%     [t x y a win nota veup hyp ttot] = VEVdriveArray(DataOut,1,criteriaUA); %DataOut{1}(criteriaUA==1)
%     else
%     [t x y a win nota veup hyp ttot] = VEVdriveArray2(DataOut{1}(criteriaUA==1)); %DataOut{1}(criteriaUA==1)
%     end
%   [t,x,y,a,~,nota,~,hyp] = VEVdriveArray2(DataOut); %DataOut{1}(criteriaUA==1)
%    BreathDataTable3.notaOrig=nota;
    %Veupnea = median(veup)*60; %L/min
    x_=1;

    arthres = ArThresPSG1/100; % taken from window summary data, median of means
   
    if 1 %set minimum for arthres PER SE
        if arthres<1.00
            arthres=1.00;
        end
    end
    arthresPSG2 = arthres;
    if arthresPSG2<1.05
        arthresPSG2=1.05;
    end
       if 0
           x = x(nota==1);
           y = y(nota==1);
       else
           x = x(nota==1&hyp~=4);
           y = y(nota==1&hyp~=4);
       end
       
       %xnew=x;
       %xorig=x;
       
       % criteria for sleep stages in breath analysis.
       if 0 && settings.breathlevelinclusion %hypok not passed into here; settings.breathlevelinclusion is default 1 for physiology
           %hyp = [0 1 2 3 4 0 1 2 4 5]';
           %hypok = [0 1 2 3];
           criteriabreath = sum(hyp==hypok,2)>0; % criteriaUA for sleep stages in breath analysis.
           %criteriabreath=criteriabreath(nota==1&hyp~=4);
           x = x(criteriabreath);
           y = y(criteriabreath);
       end
   
    if ploton % set up the space for this plot
        figure(1);
        subplot(1,3,1);
        plot(1,1,'marker','none');  hold('on');
        plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
        plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
    end
    
     if ~isfield(settings,'SummaryEndogramFilt121') 
         filt121 = 1;  %default on, since 20201027  in SummaryAnalysisOnePes, and since 20210712 in SummaryAnalysisOne
     else
         filt121 = settings.SummaryEndogramFilt121;
     end
        
    % note: more plotting is done in called Fn
    [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI]=...
        VEVdriveAnalysis(100*x,100*y,100*arthresPSG2,ploton,100,settings.Nciles,Nbootstrap,filt121);
    VpassiveCI=VpassiveCI(:);
    VactiveCI=VactiveCI(:);
    
    if ploton
        figure(1);
        subplot(1,3,1); % finish up this plot space
        set(gca,'fontsize',SAfontsize)
        xlim([0 max([prctile(100*x,95) 300])]);
        ylim([0 max([prctile(100*y,95) 100])]);
        ylabel('Ventilation (%eupnea)'); xlabel('Vdrive model (%eupnea)');
        hold('off');
        try
        title([PlotString, ' AHI =',num2str(round(AHIstate),'%u')]);
        catch me
        end
        %title(['AHI=',num2str(round(AHIstate),'%u')]);
        
        subplot(2,6,10);
        xedges=0:0.1:1.4;
        dx=0.1;
        h=histogram(100*y(x>(1-dx)&x<(1+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
        box('off');
        xlim(100*[min(xedges) max(xedges)]);
        %ylim([0 max([0.2;h.Values(:)])]);
        xlabel('Vpassive (%)')
        ylims=get(gca,'YLim');
        hold('on');
        i=6;
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        subplot(2,6,11);
        h=histogram(100*y(x>(arthresPSG2-dx)&x<(arthresPSG2+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
        box('off');
        xlim(100*[min(xedges) max(xedges)]);
        %ylim([0 max([0.2;h.Values(:)])]);
        xlabel('Vactive (%)')
        ylims=get(gca,'YLim');
        hold('on');
        i=7;
        plot(eval(varlist{i})*[1 1],ylims,'k');
        plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        %saveas(fig, [settings.directoryout, settings2.settings.filename], 'png');
    end
    
    Vcomp = Vactive-Vpassive;
    VcompCI = [VactiveCI(1)-VpassiveCI(2);VactiveCI(2)-VpassiveCI(1)];
    VcompCI = sort(VcompCI);
    
        
else
    Vpassive=NaN;
    Vactive=NaN;
    Vcomp=NaN;
    VpassiveCI=[-Inf;Inf];
    VactiveCI=[-Inf;Inf];
    VcompCI=[-Inf;Inf];
    %VEdeciles=[];
    %Vdrivedeciles=[];
    VEdeciles=NaN(10,1); % DLM changed from [] to NaN, because SummaryAnalysisN was getting cranky about empty not matching left hand side
    Vdrivedeciles=NaN(10,1); % rather than '10', this should actually be 'settings.Nciles'
end