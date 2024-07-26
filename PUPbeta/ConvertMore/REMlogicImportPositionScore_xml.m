

% Instructions from Convert for position table info
% PositionT (optional) is a table that can be created from position event
% lists. Must contain "Time" (same units as the Time signal, describes start times of each position 'event') and "Codes"
% (must be Default codes i.e. use SystemPos='Default' to translate back to
% text positions: 1=Supine,2=Left,3=Right,4=Prone,5=Unknown,6=Upright
clear all

study='MARIPOSA-1040-015-PSG2'
path = 'C:\Users\Apnimed Analysis\Apnimed Dropbox\Apnimed\F -- R&D\AD109 Trial APC005-MARIPOSA\08 Data and Reports\All EDF and scoring files\Source'
filename=[path '\' study '\' study '.edf.xml'] 

temp=readstruct(filename);

for i=1:length(temp.Events.Event)
idx(i)=contains(temp.Events.Event(i).Type.Text  ,'POSITION');
end

varList = {'POSITION-LEFT'...
'POSITION-PRONE'...
'POSITION-RIGHT'...
'POSITION-SUPINE'...
'POSITION-UNKNOWN'...
'POSITION-UPRIGHT'};
PosCodes = [2 4 3 1 5 6];

idx2=find(idx==1)
PositionT=table
for j = 1:length(idx2)
PositionT.Position(j) = temp.Events.Event(idx2(j)).Type.Text;
PositionT.PosTime(j) = temp.Events.Event(idx2(j)).StartTime.Text;
idxP = strcmp(PositionT.Position(j), varList);
PositionT.Codes(j)=PosCodes(idxP);
end

PositionT.Time = mod(datenum(PositionT.PosTime),1)*86400;
PositionT.Time(PositionT.Time<86400/2)=PositionT.Time(PositionT.Time<86400/2)+86400 %this is fix for time past midnight
                
                