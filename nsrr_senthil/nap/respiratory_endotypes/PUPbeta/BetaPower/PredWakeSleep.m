function [WSpredlogit,Tselect,EEGRef] = PredWakeSleep(Tselect,mdlA,RefTable,Exclude)
logit = @(p) log(p./(1-p));
[EEGRef,Tselect] = EEGreference(Tselect,RefTable,Exclude);
WSpredlogit = logit(predict(mdlA,Tselect));                
Tselect.WSpredlogit = WSpredlogit;