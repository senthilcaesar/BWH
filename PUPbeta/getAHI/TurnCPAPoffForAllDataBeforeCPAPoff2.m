function CPAPoff = TurnCPAPoffForAllDataBeforeCPAPoff2(CPAPoff)
CPAPoff_original = CPAPoff;

% Find end of long period of CPAPon (i.e. not off)    
CPAPoff_diff = diff(CPAPoff);
CPAPswitchOff = find(CPAPoff_diff == 1);
CPAPswitchOn = find(CPAPoff_diff == -1);

% if CPAP starts OFF (first index of ON lower than OFF): add 1 
if CPAPswitchOn(1) < CPAPswitchOff(1)
    CPAPswitchOff = [1;CPAPswitchOff];
end

% if CPAP ends OFF add length(CPAPoff) to CPAPon
if CPAPswitchOn(end) < CPAPswitchOff(end)
    CPAPswitchOn = [CPAPswitchOn;length(CPAPoff)];
end

durationCPAPoff = CPAPswitchOn - CPAPswitchOff; % find biggest chunk of cpapoff and remove all data before
[~,maxdurationOFF_i] = max(durationCPAPoff);
CPAPoff(1:CPAPswitchOff(maxdurationOFF_i)) = false;

if 0
    figure(1), plot(CPAPoff_original), ylim([-5,5]), hold on 
    plot(CPAPoff), hold off
end