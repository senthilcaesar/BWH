function [ODI3,ODI4]=CalcODI(SaO2,dt,Include)

ploton = 0;
if ploton==1
    hold('on');
    Time = (0:dt:dt*(length(SaO2)-1))';
end

% span = 10;
% searchrange = span/2;
% N=length(SaO2);


%% smooth data in sections, bypassing NaN's

%SaO2x = SaO2;
%SaO2x(isnan(SaO2))=nanmedian(SaO2);
try
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

MagDown = SaO2(SaO2MaxIdx(1:end-1))-SaO2(SaO2MinIdx(1:end)); % magnitude of drop
MagUp = SaO2(SaO2MaxIdx(2:end))-SaO2(SaO2MinIdx(1:end)); % Magnitude of recovery
% MagAv = (MagDown+MagUp)/2;

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

searchmaxleftrange = 120; % search range of 2min for pre and post baseline 
searchmaxrightrange = 120;

searchmaxleftrangei = round(searchmaxleftrange/dt);
searchmaxrightrangei = round(searchmaxrightrange/dt);
%%
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
    SpO2prei(i) = (SaO2MinIdx(i)-ileft)+Il-1;
    [valR,~]=max(SaO2(SaO2MinIdx(i):SaO2MinIdx(i)+iright));
    Ir = find(SaO2(SaO2MinIdx(i):SaO2MinIdx(i)+iright)==valR,1,'first'); % find the post baseline after sao2minidx
    SpO2posti(i) = SaO2MinIdx(i)+Ir-1;
end
    
MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);

Is = [SpO2prei SaO2MinIdx SpO2posti];
if ploton
    plot(Time(Is),SaO2(Is),'.')
end

%%
I = MagDown<3; % deleting all the sao2 drops less than 3%

SaO2MinIdx(I)=[];
SpO2prei(I)=[];
SpO2posti(I)=[];
MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);
Is = [SpO2prei SaO2MinIdx SpO2posti];
if ploton
    plot(Time(Is),SaO2(Is),'.')
end
Duration_hr=sum(Include==1)*dt/3600;
Ndesats = sum(Include(SaO2MinIdx,:));
ODI3 = Ndesats./Duration_hr;
%%
I = MagDown<4; % deleting all the sao2 drops less than 4%

SaO2MinIdx(I)=[];
SpO2prei(I)=[];
SpO2posti(I)=[];
MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);
Is = [SpO2prei SaO2MinIdx SpO2posti];
% if ploton
% plot(Time,SaO2,Time(Is),SaO2(Is),'o');
% end

Duration_hr=sum(Include==1)*dt/3600;
Ndesats = sum(Include(SaO2MinIdx,:));
ODI4 = Ndesats./Duration_hr;

catch
    ODI3=NaN;
    ODI4=NaN;
end

