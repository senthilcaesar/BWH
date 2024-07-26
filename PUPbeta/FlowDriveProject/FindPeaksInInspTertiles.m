function ThePeaks = FindPeaksInInspTertiles(InspFlow)
%% Find the largest peak within each tertile of the Inspiratory flow

% find all the peaks
PeakFindThr=0.01*max(InspFlow);
[Flow_peaks,~] = peakdetOriginal(InspFlow,PeakFindThr);

% alternate method, use peak prominence
% [pks,locs,w,p] = findpeaks(InspFlow(range(m,1):range(m,2))); 

% set up the tertile ranges
flowlength = length(InspFlow);
tertile_length = round(flowlength/3,0);
range=([1 tertile_length; tertile_length+1 tertile_length*2; (tertile_length*2)+1 length(InspFlow)]);

% find the biggest peak in each tertile
try
PkInd = NaN(3,1);
for m=1:3
    rowsinrange = find(Flow_peaks(:,1)>=range(m,1) & Flow_peaks(:,1)<=range(m,2));
    if ~isempty(rowsinrange) % only search for peaks if in range
        [~, indx] = max(Flow_peaks(rowsinrange,2));
        PkInd(m) = rowsinrange(indx);
    end 
end
catch me
    me.getReport;
end
% trim any NaN result values
PkInd(isnan(PkInd)) = [];

% return up to three peaks
ThePeaks = Flow_peaks(PkInd,:);

% [Pk1, Idx1] = max(InspFlow(1:tertile_length));
% [Pk2, Idx2] = max(InspFlow(tertile_length+1:tertile_length*2));
% Idx2 = Idx2+tertile_length;
% [Pk3, Idx3] = max(InspFlow((tertile_length*2)+1:end));
% Idx3 = Idx3+tertile_length*2;
% ThreePeaks = [Pk1 Idx1; Pk2 Idx2; Pk3 Idx3];
end