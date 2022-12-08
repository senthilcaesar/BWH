function [ output ] = BasicFilt( input, MA )
%% BASICFILT is a very very basic zero phase moving average filter
% it's quick, it's dirty... it's also very cheap.
%
% input - is the input stream
% MA - is the number of samples to apply moving average filter to
% output - is the filtered output

B = 1/MA*ones(single(MA),1);
output = filtfilt(B,1,input); % all-zero filter (FIR)
%output = filtfilt(1,B,input); % all-pole filter


end
