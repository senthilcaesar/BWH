function [deltaci,ci,medianx]=medianCI(data)
%data=[5.5 5.5 5.5 5.5 5.5 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10];
ci=bootci(length(data)*20,@nanmedian,data);
deltaci=abs(nanmedian(data)-ci);
medianx = nanmedian(data);

