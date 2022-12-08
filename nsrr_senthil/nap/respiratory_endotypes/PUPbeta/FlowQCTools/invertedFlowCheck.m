function invertedFlowCheck(patientNumber, BreathDataTable)
    windowCount = 0;
    invertCount = 0;
    for i = 1:length(BreathDataTable) % iterate through each window
        breathTable = BreathDataTable{1,i};
        try
            isnan(breathTable);
        catch ME
            VTe = breathTable{:, 'VTe'};
            VTi = breathTable{:, 'VTi'};
            VT = breathTable{:, 'VT'};
            apneaBreaths = breathTable{:, 'ApneaB'};
            VTeprev = VTe(1:end-1);
            VTicurr = VTi(2:end);
            VTecurr = VTe(2:end);
            temp = apneaBreaths(2:end);
            VTeprev(temp==1)=[];
            VTicurr(temp==1)=[];
            VTecurr(temp==1)=[];
            diffcurr = median(abs((VTicurr-VTecurr)))/mean(VT);
            diffprev = median(abs((VTicurr-VTeprev)))/mean(VT);
            if diffcurr>diffprev
                %disp(['Warning: Flow signal appears upside down, F=' num2str(-100*(diffprev-diffcurr),2)]);
                invertCount = invertCount + 1;
            end
            windowCount = windowCount + 1;
        end
    end
    invertPercentage = invertCount / windowCount * 100;
    disp(['Patient ' num2str(patientNumber) ' inverted flow warnings ' num2str(invertCount) '/' num2str(windowCount) ' windows ' num2str(invertPercentage) '%']);
end