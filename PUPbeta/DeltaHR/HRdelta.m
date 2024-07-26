function HRresponse=HRdelta(HR,Start,End,Time)
%% Ali's HR summary code   

     percntiles = prctile(reshape(HR,[],1),[1 99.99]); %1th and 99.99th percentile
%     
%     %remove outlier values
     HR(HR<percntiles(1) | HR > percntiles(2)) = NaN;
% 
     baseline1=nanmean(HR,2); % average of all 200s
     baseline2=nanmin(HR(:,round(nanmean(Start))+1:round(nanmean(End))+1),[],2); % min value during events
%     
%     %%%Default search window is between [event end - eventDur/2 : event
%     end + 50 seconds]  SS thinks this would get the next event bring it
%     in a little e.g. 15 s
     defaultSrchWin=round(nanmean(End))-round(nanmean(End)-nanmean(Start))/2+1:151;
    
    %%%Modify the end of search window using average HR
    %%%Ensembled avg HR
    AvgHR=nanmean(HR,1); %Ensembles.HRfilt

    if sum(~isnan(baseline1))>=5
        [maxS,minS] = peakdetOriginal(AvgHR,0.25,Time);
        %%%Max response occurs somewhere between event end -10 and event end
        %%%+50 seconds
        MaxWin=50;
        if ~isempty(maxS) & ~isempty(minS)
            max_resp=maxS(maxS(:,1)>=-10 & maxS(:,1)<=MaxWin,:);
            max_resp=max_resp(max_resp(:,2)==max(max_resp(:,2)),:);

            if ~isempty(max_resp)
                min_pre=minS(minS(:,1)<max_resp(1,1),:);
                if ~isempty(min_pre)
                    min_pre=min_pre(end,:);
                end
                min_post=minS(minS(:,1)>max_resp(1,1),:);
                if ~isempty(min_post)
                    min_post=min_post(1,:);
                end
                if ~isempty(min_pre) && ~isempty(min_post)
                    defaultSrchWin=find(Time>=min_pre(1,1) & Time <= min_post(1,1));
                end
            end
        end
    end
    defaultSrchWin=round(defaultSrchWin);

    maxHR=nanmax(HR(:,defaultSrchWin),[],2); % max value during search window
    
    HRresponse=[];
    HRresponse.DeltaFromMean=maxHR-baseline1;
    HRresponse.DeltaFromMin=maxHR-baseline2;
    HRresponse.MeanBaseline=baseline1;
    HRresponse.MinBaseline=baseline2;
    HRresponse.SearchLeft = repmat(defaultSrchWin(1)-101,length(baseline1),1);
    HRresponse.SearchRight = repmat(defaultSrchWin(end)-101,length(baseline1),1);
    HRresponse = struct2table(HRresponse);

    
    
%     dHR_meanbsline=nanmean(maxHR-baseline1);
%     dHR_minbsline=nanmean(maxHR-baseline2) %this is the dHR in AJRCCM
%     IndHRRep_mean=maxHR-baseline1;
%     IndHRRep_min=maxHR-baseline2; %this is the dHR in AJRCCM

    