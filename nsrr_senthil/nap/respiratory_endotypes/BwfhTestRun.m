clear all; close all; clc;

sourcefolder='U:\working\puptest\apr14\'; 
%'U:\working\bwfh\';
files = dir(fullfile(sourcefolder, '*.edf'));

for jj=1:length(files)
    jj
    clear id
    resp_output='U:\working\puptest\apr14\';
    id=extractBefore(files(jj).name,'.edf');
    edfname1=[sourcefolder id '.edf'];
    xmlname1=[sourcefolder id '.annot'];
    codedir= 'U:/luna/testing/nsrr/nap/respiratory_endotypes/';
    StartHere(resp_output,id,edfname1,xmlname1)
    
end
