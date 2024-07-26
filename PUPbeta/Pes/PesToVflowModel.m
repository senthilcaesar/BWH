function [Vdot_intended]=PesToVflowModel(parameters,xdata)

Pmus = xdata.data;
time = xdata.time;
dt = xdata.dt;

Pmus = Pmus - prctile(Pmus,parameters(5));

%initialize with zeros
Vdot_intended=0*time;
Vol=0*time;

R = parameters(1); %respiratory system resistance
E = parameters(2); %respiratory system elastance
Vol(1) = parameters(3); %initial volume at time 0
Vleak = parameters(4); %flow leak or volume offset

%Loop to simulate intended flow / intended volume traces 
for i=1:length(time)
    Vdot_intended(i) = (-(-Pmus(i))-(E*Vol(i)))/R; %intended flow
    if i<length(time)
        Vol(i+1) = Vol(i) + Vdot_intended(i)*dt; %intended vol
    end
end

Vdot_intended = Vdot_intended + Vleak; %intended flow including flow leak (to best match recorded data if needed)

%plot
if 1
figure(101); 
ax101(1)=subplot(3,1,1); plot(time,-Pmus); ylabel('Pes');
ax101(2)=subplot(3,1,2); plot(time,Vdot_intended); ylabel('PesFlow');
ax101(3)=subplot(3,1,3); plot(time,Vol); ylabel('PesVol');
linkaxes(ax101,'x');
end


