function [signal]=normalise_minus1to1(signal)
%NORMALISE_minus1TO1 Return the input signal normalised to range min=-1 and max=1
signal = -1 + 2.*(signal-min(signal))/(max(signal)-min(signal));
% alternative method to normalize to -1...1
% signal = ((signal-min(signal))./(max(signal)-min(signal)) - 0.5 ) *2;
end

% to "de-normalize", apply the calculations in reverse
% original_signal = (signal./2+0.5) * (max(signal)-min(signal)) + min(signal)