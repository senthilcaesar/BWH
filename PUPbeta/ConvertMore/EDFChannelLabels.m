function [Label,Transducer,Fs] = EDFChannelLabels(fname)
fid = fopen(fname,'r');
fseek(fid,252,-1);
M  = sscanf(char(fread(fid,4,'char')'),'%d');
for m=1:M
    fseek(fid,(m-1)*16+256,-1);
    Label{m} = char(fread(fid,16,'char')'); %%%%%%%%%
end

for m=1:M
    fseek(fid,(m-1)*80+256+M*16,-1); 
    Transducer{m} = fread(fid,80,'*char')';
end

fseek(fid,244,-1); 
Block = str2double(fread(fid,8,'*char')');
for m=1:M
    fseek(fid,(m-1)*8+256+M*16+M*80+M*120,-1); 
    Fs(m) = str2double(fread(fid,8,'*char')')/Block;
end
fclose(fid); % Close file
