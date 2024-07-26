function [Pred,Y] = EndotypePredictOAEfficacy(DataTable,BamagoosModel)

%DataTable expects data in the following order

arthres_=DataTable.ArThres; arthres_(arthres_<100)=100;
DataTable.ArThres0p5 = 100+(100*((arthres_-100)/100).^0.5);

Vpassive_=DataTable.Vpassive; Vpassive_(Vpassive_>100)=100;
DataTable.Vpassive0p5 = 100-(100*((100-Vpassive_)/100).^0.5);

%%
varlist = DataTable.Properties.VariableNames(Model.Ilist);

%%
labels = varlist;
Ain = DataTable{:,Model.Ilist} - Model.MeanValuesTable;
if Model.includequadraticterms
J=length(Model.Ilist);
Q0=tril(ones(J,J),-Model.ignoresquaredterms); %myf2 = @(A) K + [A]*L + sum(([A]*Q) .* [A], 2);
%A2 contains quadratic terms e.g. a^2 + b^2 + ab (for J=2)
A2=[];

for j=1:J
    for i=1:J
        if Q0(i,j)~=0
            A2=[A2,Ain(:,i).*Ain(:,j)];
            labels = [labels [varlist{j} '_x_' varlist{i}]];
        end
    end
end
Ain = [Ain A2]; 
end 

%%
DataTableLarge = array2table(Ain(:,Model.Ilisttest));
DataTableLarge.Properties.VariableNames = labels(Model.Ilisttest);

A = [ones(size(DataTableLarge,1),1) DataTableLarge{:,:}];
Y = A*Model.Beta;
Pred = Y>Model.Cutoff;


%%
