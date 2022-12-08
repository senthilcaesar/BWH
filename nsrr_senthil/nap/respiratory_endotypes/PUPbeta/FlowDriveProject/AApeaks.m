function [NumberofPeaks,fp1,fp2,fp3,AA_IsTerminalPeak, Flow_peaks, Flow_troughs]= ...
    AApeaks(InspFlow,BB_I_end,PIF) 
    doAApeaks=0; %default, overwritten below if peaks found
    AA_IsTerminalPeak=0;
    
    PeakFindThr=0.01*max(InspFlow);
    [Flow_peaks,Flow_troughs] = peakdetOriginal(InspFlow,PeakFindThr); %%max e min of flow inspiration [index, value]
    Flow_troughs=[1 InspFlow(1);Flow_troughs]; %%add start and end of insp as minima
    Flow_troughs=[Flow_troughs; BB_I_end InspFlow(end)];

    if ~isempty(Flow_peaks)
        while Flow_peaks(end,1)>=Flow_troughs(end,1)%%last peak must be before a min
            Flow_peaks(end,:)=[];
            if isempty(Flow_peaks)
                break
            end
        end
    end   
    if ~isempty(Flow_peaks)
        while Flow_peaks(1,1)<=Flow_troughs(1,1) %%first peak must be after a min
            Flow_peaks(1,:)=[];
            if isempty(Flow_peaks)
                break
            end
        end     
        jj_tr=[];
        for iii=2:size(Flow_troughs,1) %%find each peak before and after the first min (excluded the starting point of the breath)
            fltemp=find(Flow_peaks(:,1)>Flow_troughs(iii-1,1) & Flow_peaks(:,1)<Flow_troughs(iii,1),1);
            if isempty(fltemp)
                jj_tr=[jj_tr;iii-1];
            end
        end
        Flow_troughs(jj_tr,:)=[];
        pkTro1=Flow_peaks(:,2)-Flow_troughs(1:end-1,2);
        pkTro2=Flow_peaks(:,2)-Flow_troughs(2:end,2);
        pkTros=sort([pkTro1;pkTro2],'descend');
        
        if length(pkTros)>=8
            PeakFindThr = max(pkTros(8),0.03);
        else
            PeakFindThr = 0.03;
        end
        itr=0; doAApeaks=1;
        while 1
            [minpkTro1,minpkTro1_i] = min(pkTro1);
            [minpkTro2,minpkTro2_i] = min(pkTro2);
            maxpkTro1 = max(pkTro1);
            maxpkTro2 = max(pkTro2);
            maxpkTro = max([maxpkTro1,maxpkTro2]);
            [minpkTro,Pkpattern] = min([minpkTro1,minpkTro2]);
            itr=itr+1;
            if isempty(minpkTro) || itr>100
                break
            else
                if minpkTro>PeakFindThr
                    break
                end
            end
            if maxpkTro<PeakFindThr
                doAApeaks=0;
                break
            end
            if Pkpattern==2
                Flow_peaks(minpkTro2_i,:)=[];
                Flow_troughs(minpkTro2_i+1,:)=[];
            elseif Pkpattern==1
                Flow_peaks(minpkTro1_i,:)=[];
                Flow_troughs(minpkTro1_i,:)=[];
            end
            %recalculate
            pkTro1=Flow_peaks(:,2)-Flow_troughs(1:end-1,2);
            pkTro2=Flow_peaks(:,2)-Flow_troughs(2:end,2);
        end
    end   
  
    fp1=NaN;fp1_idx=NaN; fp2=NaN;fp2_idx=NaN; fp3=NaN;fp3_idx=NaN;
    NumberofPeaks=0;
    
    if doAApeaks==1
        if ~isempty(Flow_peaks)
            i_temp=find(Flow_peaks(:,1)<=BB_I_end/3);
            if ~isempty(i_temp)
                NumberofPeaks=NumberofPeaks+1;
                [fp1,idx_temp]=max(Flow_peaks(i_temp,2));
%                 fp1_idx=Flow_peaks(i_temp(idx_temp),1)/Fs;
            end
            i_temp=find(Flow_peaks(:,1)>BB_I_end/3 & Flow_peaks(:,1)<=2*BB_I_end/3);
            if ~isempty(i_temp)
                NumberofPeaks=NumberofPeaks+1;
                [fp2,idx_temp]=max(Flow_peaks(i_temp,2));
%                 fp2_idx=Flow_peaks(i_temp(idx_temp),1)/Fs;
            end
            i_temp=find(Flow_peaks(:,1)>2*BB_I_end/3 & Flow_peaks(:,1)<3*BB_I_end/3);
            if ~isempty(i_temp)
                NumberofPeaks=NumberofPeaks+1;
                [fp3,idx_temp]=max(Flow_peaks(i_temp,2));
%                 fp3_idx=Flow_peaks(i_temp(idx_temp),1)/Fs;
            end
%             [peakFlow,idx_temp]=max(Flow_peaks(:,2));
%             TPIF_=Flow_peaks(idx_temp,1)/Fs;
            if ~isnan(fp3) && fp3==PIF
                AA_IsTerminalPeak=1;
            else 
                AA_IsTerminalPeak=0; % set to zero to differentiate 
            end
        end
    end