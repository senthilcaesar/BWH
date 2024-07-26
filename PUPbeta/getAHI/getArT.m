function [ArT,ArIInfo]=getArT(SigT,ChannelsList,EventsStr,Hypnogram,Hypnogram_t)


dt = SigT.Time(2) - SigT.Time(1);

ArT = EventRespSigtoRespT(SigT(:,find(ChannelsList==string(EventsStr))),SigT(:,1));

%problem:
ArT.TimeSinceLastEvent = [NaN; ArT.EventStart(2:end) - ArT.EventEnd(1:end-1)];

%generate Epochs signal:
ArT.Epochs = interp1(Hypnogram_t,Hypnogram,ArT.EventStart,'previous');

ArT.AASMarousal = ArT.EventDuration>3 & ArT.Epochs<4 & ArT.Epochs>=0;
ArT.nAASMarousal = ArT.EventDuration>3 & ArT.Epochs<4 & ArT.Epochs>=0;
ArT.AASMarousal(ArT.TimeSinceLastEvent<10)=0; %problem:

TSTmin = sum(Hypnogram<4 & Hypnogram>=0)/2;

ArI = sum(ArT.AASMarousal) / (TSTmin/60);
ArI_nAASM = sum(ArT.nAASMarousal) / (TSTmin/60);

%Keep Things
ArIInfo.TSTmin = TSTmin;
ArIInfo.ArI = ArI;
ArIInfo.ArI_nAASM = ArI_nAASM;

% check on TST
% Epochs = DataEventHypnog_Mat(:,find(ChannelsList=="Epochs"));
% TSTmincheck = sum(Epochs<4&Epochs>=0)*(Time(2)-Time(1))/60;



%%
try
    
    Time = SigT.Time;
    ArSig = SigT(:,find(ChannelsList==string(EventsStr)));
    
    if EventsStr=="EventsAr"  % use original arousals
        ArIntOr = SigT.ArIntensityOrig;
        ArInt = SigT.ArIntensity;
        
        % compare with autoscored arousals
        % 1. WPr
        WPr=SigT.WPr;
        ArPr=SigT.ArPr;
        [ArT.WPrMax,ArT.ArPrMax,ArT.WSBalanceMax,ArT.ArBalanceMax,ArT.WPrMax3,ArT.ArPrMax3,ArT.WSBalanceMax3,ArT.ArBalanceMax3]...
            =getArMax(ArT,WPr,ArPr,dt);
        % 2. WPrB
        clear WPr ArPr
        WPr=SigT.WPrB;
        ArPr=SigT.ArPrB;
        [ArT.WPrBMax,ArT.ArPrBMax,ArT.WSBalanceBMax,ArT.ArBalanceBMax,ArT.WPrBMax3,ArT.ArPrBMax3,ArT.WSBalanceBMax3,ArT.ArBalanceBMax3]...
            =getArMax(ArT,WPr,ArPr,dt);
        
        
    elseif EventsStr=="EventsArWS"  % use WS A
        ArIntOr = SigT.ArIntensityOrigWS;
        ArInt = SigT.ArIntensityWS;
        
        WPr=SigT.WPr;
        ArPr=SigT.ArPr;
        
        [ArT.WPrMax,ArT.ArPrMax,ArT.WSBalanceMax,ArT.ArBalanceMax,ArT.WPrMax3,ArT.ArPrMax3,ArT.WSBalanceMax3,ArT.ArBalanceMax3]...
            =getArMax(ArT,WPr,ArPr,dt);
        
        
    elseif EventsStr=="EventsArWSB"  % use WS B
        ArIntOr = SigT.ArIntensityOrigWSB;
        ArInt = SigT.ArIntensityWSB;
        
        WPr=SigT.WPrB;
        ArPr=SigT.ArPrB; %incorrect until 2/5/2021
        
        [ArT.WPrBMax,ArT.ArPrBMax,ArT.WSBalanceBMax,ArT.ArBalanceBMax,ArT.WPrBMax3,ArT.ArPrBMax3,ArT.WSBalanceBMax3,ArT.ArBalanceBMax3]...
            =getArMax(ArT,WPr,ArPr,dt);

    end
    
    [ArIInfo,ArT]=getArT2(ArT,Time,ArIntOr,ArInt,ArIInfo);
catch me
    
end
%%
ArTinfo=struct2table(ArIInfo(:));

