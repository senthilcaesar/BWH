if 1
    l=10370;
    r=10387;
    range = find(Time1>l&Time1<r);
else
    range = find(Time1>-Inf&Time1<Inf);
end

[I,~,~,~,~]=Vflowanalysis2(FlowEdi1(range),Time1(range),0.008,0);

Ccw = 0.2;

%%
[Vdot_intended,parameters,rsquared]=PesToVflowFit(Flow1(range),Pes1(range),Time1(range),Ccw);

%%
range = find(Time1>-Inf&Time1<Inf);
[Vdot_intended]=PesToVflowRun(Flow1(range),Pes1(range),Time1(range),Ccw,parameters);