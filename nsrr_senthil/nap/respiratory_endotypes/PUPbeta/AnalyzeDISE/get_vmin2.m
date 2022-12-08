% Save this function as "get_vmin2.m"
% This gets the volume minimum (end exp lung volume) from the tidal volume

function [vmin,ivmin] = get_vmin2(vol,ti,te,incr,fs)


% allocate memory for collection vectors
ivmin = []; 
vmin = [];
ivmax = [];
vmax = [];

% find the first minimum of the vol signal
[volmin,indmin] = min(vol(1:1+round(fs*(ti+te+ti))));  %search ti+te+ti seconds ahead
%[volmin,indmin] = min(vol(1:1+round(fs*incr*te)));  %search incr*te seconds ahead
vmin = [vmin volmin]; % find the minimum volume
ivmin = [ivmin indmin]; % find the index of the minimum

while (indmin+(fs*ti*2)+(fs*2*te) < length(vol)) % stop looping when indmin is within 2 breaths from the end
    [volmax,indmax] = max(vol(indmin:indmin+round(fs*ti*incr))); % search ahead by incr ti's
    indmax = indmin + indmax-1;
    vmax = [vmax volmax];
    ivmax = [ivmax indmax];

    [volmin,indmin] = min(vol(indmax:indmax+round(fs*incr*te))); % search ahead by incr te's    
    indmin = indmax + indmin-1;
    vmin = [vmin volmin];
    ivmin = [ivmin indmin];
end

vmin = vmin';
ivmin = ivmin';