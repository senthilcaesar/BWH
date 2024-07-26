function info = NDFread1(filename)

fid=fopen(filename);

fseek(fid, 0, -1); %fseek(fileID, offset, origin)
%c = fread(fid);
aa=fread(fid,'*int8'); 
N = length(aa); 
%c = char(aa(1:N))';

%new
c = char(aa(1:end))';   
    clear aa;
d = c(1:2:end);   
k2 = strfind(string(d),"T"); %search for end of channel header
k2(k2<9|k2>(length(d)-12))=[];
iSessionA = k2(find(d(k2-8)=='2' & d(k2+7)=='.')); %Expected string series: 2-------T------.------

%find which ones are actually numbers
Imatrix = iSessionA(:)+[-8:13];
kx = d(Imatrix); %view this to debug - each row should be a sensible date
%kx((find(d(k2-8)=='2' & d(k2+7)=='.')),:)
clear kxgood
for i=1:size(kx,1)
    tempnum = str2num(kx(i,1:8));
    tempnum2 = str2num(kx(i,10:15));
    tempnum3 = str2num(kx(i,17:22));
    kxgood(i)=0; %default
    if ~isempty(tempnum) && tempnum>19000000 && tempnum<50000000 && ~isempty(tempnum2) && ~isempty(tempnum3)
    kxgood(i)=1;
    end
end
iSession = iSessionA(logical(kxgood));

Nsessions = length(iSession);



% 
% 
% %c = fscanf(fid,'%c',N);
% k2_ = strfind(string(c),"</Channel>"); %search for end of channel header
% if isempty(k2_)
    I=1:2:N;
% else
%     I=1:1:N;
% end
%d = c(I);
k2 = strfind(string(d),"</Channel>"); %search for end of channel header
info=[];
varList={'Format','Scale','Offset','SamplingRate','Label','SourceType','Unit','ChannelNumber','IsDC','Type'};
doubleList=[0 1 1 1 0 0 0 1 1 0 ];

for i=1:length(varList)
    try
        k1 = strfind(string(d),strcat("<" ,varList{i} ,">"));
        k11 = strfind(string(d),strcat("</" ,varList{i} ,">"));
        info = setfield(info,varList{i},d(k1+length(varList{i})+2:k11-1));
        if doubleList(i)==1
            info = setfield(info,varList{i},str2num(getfield(info,varList{i})));
        end
    catch
    end
end

%k2 = strfind(string(d),"</Channel>"); %search for end of channel header
%k4 = strfind(string(d),"T"); %find the first "T" after the end of header
%k4(k4<k2)=[];
%k4=k4(1);
k4=iSession(1);
info.startTimeStr = [d(k4+1:k4+2) ':'  d(k4+3:k4+4) ':'  d(k4+5:k4+10)];



ScaleFactor=info.Scale;
ScaleOffset=info.Offset; 

if length(iSession)==1
    k5 = k4 + 23; %first data point is after the end of the header text here
    idx = I(k5)+1; %first data point is here
    fseek(fid, idx, -1); %go to first data point
    data=fread(fid,['*int' info.Format(end-1:end)]); %e.g. data=fread(fid,'*int16');
else
    data=[];
    for j=1:length(iSession)
        if string(info.Format(end-1:end))=="32"
            Nsamplesperinx = 2;
        else
            Nsamplesperinx = 1;
        end
        iSessionActual = I(iSession+23)+1;
        fseek(fid, iSessionActual(j), -1); %go to first data point
        Lsession = round(((diff([iSession N-1])-1)-34)/Nsamplesperinx); %magic number 34 is from empiricial testing (diff offset was 17 for 32-bit data)
        Lsession(end)=Inf;
        %remainders = mod((diff([iSession N-1])-1)/Nsamplesperinx,1) - checked these were zero
        data1=fread(fid,Lsession(j),['*int' info.Format(end-1:end)]); %e.g. data=fread(fid,'*int16');
        data=[data;data1];
        %info.data=double(data1)*ScaleFactor + info.Offset;
        info.Nsessions=length(iSession);
        info.SessionLength(j)=length(data1);
        k4=iSession(j);
        info.SessionStart{j}=[d(k4+1:k4+2) ':'  d(k4+3:k4+4) ':'  d(k4+5:k4+10)];
    end
   %    figure(90)
%    plot(double(data)*ScaleFactor + info.Offset) %29665
end


formatIn = 'HH:MM:SS.FFF';
StartTime = mod(datenum(info.startTimeStr,formatIn),1)*86400; % SO 20230328, get start and end of EDF, needed below
if StartTime<86400/2,StartTime=StartTime+86400; end
info.StartTime=StartTime;

if 1 %need to check how this works in a channel with offset still, should be one of the following ways:
info.data=double(data)*ScaleFactor + info.Offset;
else
info.data=(double(data)+info.Offset)*ScaleFactor;
end

if 0
   figure(89)
   plot(info.data) 
end


info.fsTrue=info.SamplingRate;

info.Label(any(info.Label==[ ' !-?*/()[]']',1))='_';
if ~isempty(str2num(info.Label(1)))
    info.Label = strcat('AA_',info.Label);
end

info.Nsessions=Nsessions;

fclose(fid);

%%
