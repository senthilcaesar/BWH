function signalfiltered = nanfilter(B,A,signal,zerophase)

% input order must be: B,A,signal not signal,B,A. If there is an error
% please update the code outside this function to the B,A,signal order.


if ~exist('zerophase','var')
    zerophase=0;
end

% cutoff = 1/(2*pi*0.05); %lowpass, smooth
% filter_order = 15;
% [B,A] = butter(filter_order,[cutoff]/(1/dt/2),'low');
% zerophase=1;
% signal = SignalAllrs;

% signal = [ 0 0 NaN NaN]'

signalfiltered = NaN*signal; %copies over NaN and initializes; actual data will be overwritten in loop

Ix=~isnan(signal);
if sum(Ix)==length(Ix)
    I1=1;
    I2=length(Ix);
else
    I1=find([NaN;diff(Ix)]==1);
    I2=find([diff(Ix)]==-1);
    [I1,I2] = TidyStartEndEventList(I1,I2,length(signal));
end
%[I1,I2] = TidyStartEndEventList(I2,I1,length(signal));

Ilength = I2-I1+1;
samplelimit = 1+3*(max([length(A) length(B)])-1); %possibly imprecise, only needed for filtfilt

for i=1:length(I1)
    if Ilength(i)<samplelimit
        continue
    end
    if zerophase
        signalfiltered(I1(i):I2(i)) = filtfilt(B,A,signal(I1(i):I2(i)));
    else
        signalfiltered(I1(i):I2(i)) = filter(B,A,signal(I1(i):I2(i)));
    end
end



if 0 %this version was in EdiFilt subdirectory, unclear which was being used in PUPbeta
    
    signalfiltered = NaN*signal; %copies over NaN and initializes; actual data will be overwritten in loop

    Ix=isnan(signal);
    I1=find([diff(Ix)]==1);
    I2=find([NaN;diff(Ix)]==-1);
    [I1,I2] = TidyStartEndEventList(I2,I1,length(signal));

    Ilength = I2-I1+1;
    samplelimit = 1+3*(max([length(A) length(B)])-1); %possibly imprecise, only needed for filtfilt

    for i=1:length(I1)
        if Ilength(i)<samplelimit
            continue
        end
        if zerophase
            signalfiltered(I1(i):I2(i)) = filtfilt(B,A,signal(I1(i):I2(i)));
        else
            signalfiltered(I1(i):I2(i)) = filter(B,A,signal(I1(i):I2(i)));
        end
    end

end

