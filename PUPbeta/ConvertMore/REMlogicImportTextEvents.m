
default='K:\HeartBeat\Source\*.*'
[filename,directory]=uigetfile(default);


clear eventstemp
fid = fopen([directory filename]);
i=1;
while 1
    eventstemp{i,:} = fgetl(fid);  %read in the data
    if eventstemp{i,:}==-1
        eventstemp(i,:) = [];
        break
    end
    i=i+1;
end
fclose(fid);   %close the file

option=2;
switch option
    case 1 %fails because headers are not exactly consistent
        HeaderTextExact = {['Position' char(9) 'Time [hh:mm:ss.xxx]' char(9) 'Event' char(9) 'Duration[s]']}
        NPreHeaderlines=find(strcmp(eventstemp,HeaderTextExact))-1;
    case 2
        HeaderTextList = {'Event','Duration'};
        I=[];
        for ii=1:length(HeaderTextList)
            I=[I,contains(eventstemp,HeaderTextList{ii})];
        end
        I=sum(I,2)/length(HeaderTextList);
        NPreHeaderlines=find(I==1)-1
        HeaderTextExact = eventstemp{find(I==1)}
end

% Import Table (now that we know where the data starts)
EventsTable = readtable([directory filename],'Delimiter','\t','ReadVariableNames',1,'HeaderLines',NPreHeaderlines);

clear timetemp
for i=1:size(EventsTable,1)
            eventtemp=EventsTable.Time_hh_mm_ss_xxx_{i}; 
            if 0
            I=find(eventtemp=='	');
            if length(I)==4
                I(1)=[];
            end
            end
            %Time is in sec since the previous day's midnight.
            tempstr=eventtemp;%eventtemp(I(1)+1:I(2)-1);
            I2=find(tempstr==':');
            timetemp(i,1)=86400*(tempstr(I2(2)+4)=='A')+43200*(tempstr(I2(2)+4)=='P')+3600*mod(str2num(tempstr(1:I2(1)-1)),12)+60*str2num(tempstr(I2(1)+1:I2(1)+2)) + 1*str2num(tempstr(I2(2)+1:I2(2)+2)) + str2num(tempstr(end-3:end));
end
timetemp(timetemp<43200)=timetemp(timetemp<43200)+86400;
