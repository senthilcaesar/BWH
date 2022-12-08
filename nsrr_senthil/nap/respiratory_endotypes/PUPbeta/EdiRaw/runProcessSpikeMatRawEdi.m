% RawEdi Analysis
clear; clc; %close all;

%Select file format (Luo / Spike data)
SignalFormat = 'Spike'; %'Spike' %'Luo'

% if strcmp(SignalFormat,'Luo')
%     filenamelist = {'placeholderfilename','placeholderfilename2'};
%     pathpiece = 'E:\Dropbox (Personal)\PhenotypeDrive2018\';
% elseif strcmp(SignalFormat,'Spike')
    pathpiece = pwd;
    %pathpiece = 'J:\PEOPLE\FACULTY\SANDS\O2PSG\AdPhenotype\_Scored\';
    %pathpiece = 'E:\Dropbox (Partners HealthCare)\PAtO\Studies\Scored\';
    %pathpiece = 'C:\Users\dw46\Desktop\';
    [spikefile,path]=uigetfile('*.*','Select original spike file',pathpiece);
    
% end

% addpath(genpath([pathpiece 'Luos Data\ADIMAT'])); %read ADI toolset
%addpath(genpath([pathpiece 'PUPbeta20190404\'])); %PUPbeta
% addpath(genpath('E:\Dropbox (Personal)\PUPbeta_git\PUPbeta20190629'));
% addpath(genpath([pathpiece 'Luos Data\EditDist\'])); %external string recognition function

f = load([path 'RawEdi' '.mat']);
    
Nsigs=5;
    for EdiOpt=1:Nsigs
        
        clear EMGdi Err
        [EMGdi,Err]=EdiAnalysisSpike(f,EdiOpt,SignalFormat);
               
        savevarlist = {'EMGdi'};
        for i=1:length(savevarlist)
            eval([savevarlist{i} num2str(EdiOpt) '=' savevarlist{i} ';']);
            savevarlist2{i}=[savevarlist{i} num2str(EdiOpt)];
        end
        
        %if saving the figure (slow, ~600MB):
        %saveas(4,['subject' num2str(subject) '_' num2str(EdiOpt) '.fig']);
        if exist([path 'FiltEMGdi' '.mat'],'file')==0
        %    save(['subjectA' num2str(subject) '.mat'],savevarlist2{:},'-v7.3');          
            save([path 'FiltEMGdi' '.mat'],savevarlist2{:},'-v7.3');
        else
            save([path 'FiltEMGdi' '.mat'],savevarlist2{:},'-append','-v7.3');
        end
        
        pause(5)
        close all
    end
    
%%

OneSignalFrom5
%SaveToSpike
