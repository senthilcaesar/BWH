function [EEGRef,Tselect_out,Tcentiles] = EEGreference(Tselect,RefTable,Exclude)

newcollist = {'Pbeta','Palpha','Ptheta','Pdelta'};
centilexlist = [5 25 50 75 95];

Tcentiles = table();
for i=1:length(newcollist)
    for j=1:length(centilexlist)
        variablename = [newcollist{i} num2str(centilexlist(j)) 'p'];
        Tcentiles{:,variablename} = prctile(Tselect{~Exclude,newcollist{i}},centilexlist(j));
    end
end

%RefTable = table(Constants,Coeffs)
%How to apply:
EEGRef=0;
for i=2:length(RefTable.Constants)
    EEGRef = EEGRef + RefTable.Coeffs(i)*eval(['Tcentiles.' RefTable.Constants{i}]);
end
EEGRef = EEGRef + RefTable.Coeffs(1);

Tselect_out = Tselect;
for i=1:length(newcollist)
    eval(['Tselect_out.' newcollist{i} '_ref = Tselect.' newcollist{i} ' - EEGRef;' ]);
end
