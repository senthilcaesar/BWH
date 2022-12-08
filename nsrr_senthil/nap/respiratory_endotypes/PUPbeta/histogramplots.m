
varlist={'LG1','LGn','ArThres','Vpassive','Vactive','Vcomp'}
figure (111); clf (111);
subplot(2,3,1);
histogram(SummaryAnalysisTable.LG1,'Orientation', 'horizontal');
ylabel(varlist{1})
subplot(2,3,2);
histogram(SummaryAnalysisTable.LGn,'Orientation', 'horizontal');
ylabel(varlist{2})
subplot(2,3,3);
histogram(SummaryAnalysisTable.ArThres,'Orientation', 'horizontal');
ylabel(varlist{3})
subplot(2,3,4);
histogram(SummaryAnalysisTable.Vpassive,'Orientation', 'horizontal');
ylabel(varlist{4})
subplot(2,3,5);
histogram(SummaryAnalysisTable.Vactive,'Orientation', 'horizontal');
ylabel(varlist{5})
subplot(2,3,6);
histogram(SummaryAnalysisTable.Vcomp,'Orientation', 'horizontal');
ylabel(varlist{6})
