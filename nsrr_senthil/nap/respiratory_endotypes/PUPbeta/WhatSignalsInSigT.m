clear Tvars
for i=1:size(settings.patients,1)
    filename = [settings.ConvertedDirectory settings.patients{i,1} '.mat']
    clear SigT
    if exist(filename)==2
        disp('loading SigT')
        load(filename,'SigT');
        if ~exist('Tvars')
            Tvars = SigT(1,:);
            Tvars{i,:}=1;
            Tvars{1:i-1,:}=NaN;
            Tvars.Subj=[1:height(Tvars)]';
        else
            Tvars2 = SigT(1,:);
            Tvars2{1,:}=1;
            Tvars2.Subj=i;
            Tvars = outerjoin(Tvars,Tvars2,'MergeKeys',true);
        end
    end
end
clear Tvars2 SigT i
Tvars = sortrows(Tvars,'Subj','ascend');
Tvars = movevars(Tvars, 'Subj','Before',Tvars.Properties.VariableNames{1});
save([settings.SummaryDirectory 'Tvars.mat'],'Tvars');