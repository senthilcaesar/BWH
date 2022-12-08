function success = MakeReport(Subjects)
global settings

success = zeros(length(Subjects),1);
for i=1:length(Subjects)
    Subject = Subjects(i);
    try
    A = getData(Subject);
    A = A(Subject,:);
    
    if ~(exist([settings.workdir,'Reports'], 'dir') == 7)
        mkdir([settings.workdir,'Reports']);
    end
    
    PSGReportfilename = [settings.patients{Subject}(1:end-4)];
    fdestination = [settings.workdir 'Reports' filesep  PSGReportfilename '_Report.xlsx'];
    
    statustemp = copyfile([settings.codedir '\PSGReport\getDataReportTemplate.xlsx'],fdestination);
    
    cell0='D80'; %for filename
    cell1='D81'; %for analysis
    xlswrite(fdestination,{PSGReportfilename},'PSGReport',cell0);
    xlswrite(fdestination,A{:,:}','PSGReport',cell1);
    success(i)=1;
    catch
    success(i)=0;    
    end
    
end