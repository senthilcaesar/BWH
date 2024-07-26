function [deltaci]=medianCI(data)
%data=[5.5 5.5 5.5 5.5 5.5 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10];
ci=bootci(length(data)*20,@median,data);
deltaci=abs(median(data)-ci);

