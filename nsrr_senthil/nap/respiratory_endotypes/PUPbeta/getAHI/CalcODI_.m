function [ODI3,ODI4,ODI3_TotDur,ODI4_TotDur,Duration_hr,Timebelow90,Timebelow90_per,Timeabove98,Timeabove98_per,Flag]=CalcODI(SaO2,dt,Include)

ploton = 1;
% if ploton==1
%     hold('on');
Time = (0:dt:dt*(length(SaO2)-1))';
% end
fig = get(groot,'CurrentFigure');
if ~isempty(fig)
    close(figure(1));
end
%% smooth data bypassing NaN's
% SaO2(isnan(SaO2))=nanmedian(SaO2);
SaO2=fillmissing(SaO2,'nearest'); % fill the NaNs  nearest non-missing value
[SaO2Mins,SaO2Maxs]=peakdet(-SaO2,0.5);
SaO2MaxIdx=SaO2Maxs(:,1);
SaO2MinIdx=SaO2Mins(:,1);
while SaO2MaxIdx(1)<0
    SaO2MaxIdx(1)=[];
end
while SaO2MinIdx(end)>=SaO2MaxIdx(end)
    SaO2MinIdx(end)=[];
end
while SaO2MinIdx(1)<SaO2MaxIdx(1)
    SaO2MinIdx(1)=[];
end
%checkallpos1 = SaO2MinIdx(1:end)-SaO2MaxIdx(1:end-1); %good
%checkallpos2 = SaO2MaxIdx(2:end)-SaO2MinIdx(1:end); %good
if ploton
    figure (1); plot(Time,SaO2); hold on
end
% plot (Time(SaO2MinIdx),SaO2(SaO2MinIdx),'*r'); plot(Time(SaO2MaxIdx),SaO2(SaO2MaxIdx),'*g');

MagDown = SaO2(SaO2MaxIdx(1:end-1))-SaO2(SaO2MinIdx(1:end)); % magnitude of drop
MagUp = SaO2(SaO2MaxIdx(2:end))-SaO2(SaO2MinIdx(1:end)); % Magnitude of recovery

MagAv_thres = 2;
while 1
    [minMagDown,minMagDown_i] = min(MagDown); % taking the overall minimum in MagDown array
    [minMagUp,minMagUp_i] = min(MagUp); % taking the overall minimum in MagUp array
    [minMag,minMagpattern] = min([minMagDown,minMagUp]); % find whether minMagdown or minMagUp has the lowest value.
    
    if minMag>MagAv_thres % if minimum exceeds average threshold,stop.
        break
    end
    
    if minMagpattern==1 %down
        i=minMagDown_i;
        SaO2MaxIdx(i)=[]; % if MagDown (drop) has the lowest value, make previous max idx and current min idx to zero.
        SaO2MinIdx(i)=[];
    elseif minMagpattern==2 %up
        i=minMagUp_i;
        SaO2MaxIdx(i+1)=[]; % if Magup (recovery) has the lowest value, make next max idx and current min idx to zero.
        SaO2MinIdx(i)=[];
    end
    
    %recalculate
    MagDown = SaO2(SaO2MaxIdx(1:end-1))-SaO2(SaO2MinIdx(1:end)); % recalculate after deleting lowest values
    MagUp = SaO2(SaO2MaxIdx(2:end))-SaO2(SaO2MinIdx(1:end));
    %pause
    if isempty(SaO2MaxIdx) || isempty(SaO2MinIdx)
        break
    end
end
MagDown = SaO2(SaO2MaxIdx(1:end-1))-SaO2(SaO2MinIdx(1:end));
MagUp = SaO2(SaO2MaxIdx(2:end))-SaO2(SaO2MinIdx(1:end));

% if ploton
% plot(Time(SaO2MaxIdx),[SaO2(SaO2MaxIdx)],'o',Time(SaO2MinIdx),[SaO2(SaO2MinIdx)],'o')
% end

%% search range of 2min for pre and post baseline
SaO2Min=SaO2(SaO2MinIdx);
SaO2Min_Dur=zeros(length(SaO2Min),1);
searchmaxleftrange = 120;
searchmaxrightrange = 120;

searchmaxleftrangei = round(searchmaxleftrange/(dt));
searchmaxrightrangei = round(searchmaxrightrange/(dt));

SpO2prei = SaO2MaxIdx(1:end-1);
SpO2posti = SaO2MaxIdx(2:end);
for i=1:length(SaO2MinIdx)
    ileft = SaO2MinIdx(i)-SaO2MaxIdx(i);
    iright = SaO2MaxIdx(i+1)-SaO2MinIdx(i);
    if ileft>searchmaxleftrangei
        ileft=searchmaxleftrangei;
    end
    if iright>searchmaxrightrangei
        iright=searchmaxrightrangei;
    end
    [valL,~]=max(SaO2(SaO2MinIdx(i)-ileft:SaO2MinIdx(i)));
    Il = find(SaO2(SaO2MinIdx(i)-ileft:SaO2MinIdx(i))==valL,1,'last'); % find the pre baseline value before sao2minidx
    %     [valMed]=median(SaO2(SaO2MinIdx(i)-ileft:SaO2MinIdx(i)));
    %     Imx = find(SaO2(SaO2MinIdx(i)-ileft:SaO2MinIdx(i))==valL,1,'last'); % find the pre baseline value before sao2minidx
    %     Imo= find(SaO2(SaO2MinIdx(i)-ileft:SaO2MinIdx(i))<=valMed,1,'last');
    %     if Imo>Imx, Il=Imo; else Il=Imx; end
    SpO2prei(i) = (SaO2MinIdx(i)-ileft)+Il-1;
    
    valMean=prctile(SaO2(SaO2MinIdx(i)-ileft:SaO2MinIdx(i)),95); % find threshold value for calculating drop
    [valR,~]=max(SaO2(SaO2MinIdx(i):SaO2MinIdx(i)+iright));
    Ir = find(SaO2(SaO2MinIdx(i):SaO2MinIdx(i)+iright)==valR,1,'first'); % find the post baseline after sao2minidx
    SpO2posti(i) = SaO2MinIdx(i)+Ir-1;
    SaO2Min_Dur(i)= (SpO2posti(i)-SpO2prei(i))/(1/dt);
    %     Drop_3=valL-3;
    %     Id1=find(SaO2(SpO2prei(i):SpO2posti(i))<=Drop_3,1,'first'); % finding the duration for which sao2 dropped below 3%
    %     Id2=find(SaO2(SpO2prei(i):SpO2posti(i))<=Drop_3,1,'last');
    %
    %     if ~isempty(Id1) && ~isempty(Id2) % if there is a real 3% drop
    %         if Id1==Id2
    %             SaO2Min_Dur(i)=1/(1/dt); % drop only for 1 sample
    %         else
    %             SaO2Min_Dur(i)=(Id2-Id1)/(1/dt);
    %         end
    %     elseif isempty(Id2) && ~isempty(Id1) % if 3% drop without a recovery above baseline
    %         SaO2Min_Dur(i)=((SpO2posti(i))-(Id1+SpO2prei(i)-1))/(1/dt);
    %     elseif isempty(Id1) && ~isempty(Id2) % if 3% drop with respect to recovery
    %            SaO2Min_Dur(i)=((SpO2prei(i)+Id2-1)-(SpO2prei(i)))/(1/dt);
    %
    %     else
    %         SaO2Min_Dur(i)=NaN;
    %     end
    
end

MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);

Is = [SpO2prei SaO2MinIdx SpO2posti];
% if ploton
%     plot(Time(Is),SaO2(Is),'*')
% end

%%
Flag=0;
J=(SaO2Min_Dur<5|isnan(SaO2Min_Dur)); % deleting all sao2 drops with duration less than 5 seconds
if ~isempty(J)
    Flag=1;
    SaO2MinIdx(J)=[];
    SpO2prei(J)=[];
    SpO2posti(J)=[];
    SaO2Min_Dur(J)=[];
    MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
    MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);
    Is = [SpO2prei SaO2MinIdx SpO2posti];
else
    Flag=0;
    disp('There are no SaO2 desats which are greater than 5 s in duration');
    ODI3=NaN;ODI4=NaN; ODI3_TotDur=NaN; ODI4_TotDur=NaN; Duration_hrsum(Include==1)*dt/3600;
end

if Flag==1
    I = MagDown<3; % deleting all the sao2 drops less than 3%
    SaO2MinIdx(I)=[];
    SpO2prei(I)=[];
    SpO2posti(I)=[];
    SaO2Min_Dur(I)=[];
    
    MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
    MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);
    Is = [SpO2prei SaO2MinIdx SpO2posti];
    if ploton
        hold on; plot(Time(Is),SaO2(Is),'*')
    end
    Duration_hr=sum(Include==1)*dt/3600;
    if Duration_hr<1 % taking minutes in cases where there is not even 1 hr of sao2 data
        Duration_hr=Duration_hr*60;
        Duration_hr(:,2)=1;
    else
        Duration_hr(:,2)=0;
    end
    Ndesats = sum(Include(SaO2MinIdx,:));
    ODI3 = Ndesats./Duration_hr(:,1);
    ODI3_TotDur=sum(SaO2Min_Dur(Include(SaO2MinIdx,:)==1))/60; % in minutes
    %%
    I = MagDown<4; % deleting all the sao2 drops less than 4%
    SaO2MinIdx(I)=[];
    SpO2prei(I)=[];
    SpO2posti(I)=[];
    SaO2Min_Dur(I)=[];
    
    MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
    MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);
    Is = [SpO2prei SaO2MinIdx SpO2posti];
    if ploton
        hold on; plot(Time,SaO2,Time(Is),SaO2(Is),'o');
    end
    
    Ndesats = sum(Include(SaO2MinIdx,:));
    ODI4 = Ndesats./Duration_hr(:,1);
    ODI4_TotDur=sum(SaO2Min_Dur(Include(SaO2MinIdx,:)==1))/60;
    
    SaO2temp=SaO2(Include==1);
    Timebelow90=(length(SaO2temp(SaO2temp<=90))/(1/dt)/60);
    Timeabove98=(length(SaO2temp(SaO2temp>=98))/(1/dt)/60);
    % if Duration_hr(:,2)==0
    Timebelow90_per=((Timebelow90/60)/Duration_hr(:,1))*100;
    Timeabove98_per=((Timeabove98/60)/Duration_hr(:,1))*100;
    % else
    %     Timebelow90_per=((Timebelow90)/Duration_hr(:,1))*100;
    %     Timeabove98_per=((Timeabove98)/Duration_hr(:,1))*100;
    % end
    pause(3)
    
end



