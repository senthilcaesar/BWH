clear all

loadfile=['C:\PSG_Data\FlowDrive\Analyzed\AA_Testing_2'];
load(loadfile,'settings')
[Data,DataCI,DataN,varlist,AHItotal,Fstates] = SummaryAnalysisOne_DLM(2,settings);

fig = gcf;
str = ['C:\PSG_Data\FlowDrive\Analyzed\AA_Testing_2_Plot'];
saveas(fig, str, 'png'); 