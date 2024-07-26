%% Prepare and Save Model

DataTable = table(LG1,LGn,delay,VRA1,ArThres,Vpassive,Vactive,Vcomp);
DataTable(Exclude==1,:)=[];

DataTable.ArThres = DataTable.ArThres*100; 

arthres_=DataTable.ArThres; arthres_(arthres_<100)=100;
DataTable.ArThres0p5 = 100+(100*((arthres_-100)/100).^0.5);

Vpassive_=DataTable.Vpassive; Vpassive_(Vpassive_>100)=100;
DataTable.Vpassive0p5 = 100-(100*((100-Vpassive_)/100).^0.5);

Ilist = [10 2 8 9 4];
Model.Ilist = Ilist;
Model.varlist=varlist;
Model.MeanValuesTable = mean(DataTable{:,Ilist});
Model.includequadraticterms=1;
Model.ignoresquaredterms=0;
Model.Ilisttest=Ilisttest;
Model.Beta=Btemp;
Model.Cutoff=thresopt;
%Model.Coefficients etc = 

save Model Model





