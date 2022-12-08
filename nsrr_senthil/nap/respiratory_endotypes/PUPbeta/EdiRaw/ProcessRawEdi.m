% Analyze Luo's data March 19, 2019

clear; clc; %close all;
filenamelist = {'placeholderfilename','placeholderfilename2'};

pathpiece = 'E:\Dropbox (Personal)\PhenotypeDrive2018\';
%pathpiece = 'C:\Users\dw46\Desktop\';
addpath(genpath([pathpiece 'Luos Data\ADIMAT'])); %read ADI toolset
%addpath(genpath([pathpiece 'PUPbeta20190404\'])); %PUPbeta
addpath(genpath(['E:\Dropbox (Personal)\PUPbeta_git\PUPbeta20190629']))
addpath(genpath([pathpiece 'Luos Data\EditDist\'])); %external string recognition function


for subject=3%2:length(filenamelist)
    
    filename = filenamelist{subject};
    f = adi.readFile([pathpiece 'Luos Data\Source\' filename '\' filename '.adicht']);
        
    for EdiOpt=1:5
        
        clear EMGdi Err
        [EMGdi,Err]=EdiAnalysis(f,EdiOpt);
        
        savevarlist = {'EMGdi','Err'};
        for i=1:length(savevarlist)
            eval([savevarlist{i} num2str(EdiOpt) '=' savevarlist{i} ';']);
            savevarlist2{i}=[savevarlist{i} num2str(EdiOpt)];
        end
        
        %if saving the figure (slow, ~600MB):
        %saveas(4,['subject' num2str(subject) '_' num2str(EdiOpt) '.fig']);
        
        if exist(['subject' num2str(subject) '.mat'],'file')==0
            save(['subject' num2str(subject) '.mat'],savevarlist2{:},'-v7.3');
        else
            save(['subject' num2str(subject) '.mat'],savevarlist2{:},'-append','-v7.3');
        end
        
        pause(5)
        close all
    end
end

