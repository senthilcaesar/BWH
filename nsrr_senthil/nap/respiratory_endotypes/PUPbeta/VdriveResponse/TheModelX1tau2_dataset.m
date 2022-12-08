function [SSres_LG1,Error,Rsq,Vchem,Varousal,Penalized_Error] = TheModelX1tau2_dataset(Parameters,DataM,i_1,delayed_VEM,polyfitorder,ParametersFixed)
%not operational yet
%aim is to eventually call this function instead of TheModelX1tau2 for
%running a simultaneous fit to a full set of flow/drive data, e.g. one best
%fit to all events.

for m=1:length(DataM)
    delayedVE = delayed_VEM{m};
    Data = DataM{m};
    [SSres(m),Error{m},Rsq(m),Vchem{m},Varousal{m},Penalized_Error{m}] = TheModelX1tau2(Parameters,Data,i_1,delayedVE,polyfitorder,ParametersFixed)
end

SSres_LG1 = mean(SSres);

