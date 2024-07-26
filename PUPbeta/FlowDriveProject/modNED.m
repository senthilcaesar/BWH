function [NED1, NED2] = modNED(InspFlow, BB_I_end)

%Find peak between 0 and 30% of inspiratory time
peak0to30 = max(InspFlow(1:BB_I_end/3));

%Find min between 25 and 75% of inspiratory time
min25to75 = min(InspFlow(BB_I_end/4:3*BB_I_end/4));

%Find peak from 70 to 100% of inspiratory time
peak60to100 = max(InspFlow(2*BB_I_end/3:end));

%Compute NED from first peak to minimum
NED1 = (peak0to30 - min25to75)/peak0to30;

%Compute NED from last peak to minimum
NED2 = (peak60to100 - min25to75)/peak60to100;