function getSigTfromConvert(MrangeOverride)
% RUN StartHere.m first

% this function creates Signal table from converted files to replace DataEventHypnog_Mat .

global AMasterSpreadsheet

t_startGetData = clock;

%% read excel spreadsheet

[~,txt] = xlsread(AMasterSpreadsheet,1,'AD4:AE5000');
M = size(txt,1);

ID=1:M;

% override if mrange present
if exist('MrangeOverride')
    ID = MrangeOverride(:);
end

success=zeros(length(ID),1);

for i=1:length(ID)
    try
        clear W
        
        W = load([txt{ID(i),2} txt{ID(i),1}]);
        disp(['processing: ' txt{ID(i),1}])
        
        hasDataEventMat = isfield(W, 'DataEventHypnog_Mat');
        
        if hasDataEventMat
            W.SigT=array2table(W.DataEventHypnog_Mat);
            W.SigT.Properties.VariableNames = W.ChannelsList;
            % Field is there.  Remove it.
            W = rmfield(W, {'DataEventHypnog_Mat','ChannelsList'});
             
            save([txt{ID(i),2} txt{ID(i),1}],'-struct','W','-v7.3');
            success(i)=1;
            disp(['success: SigT created for: ' txt{ID(i),1}]);
            
        else
            % Field is not there
            disp('DataEventHypnog_Mat not present');
        end
       
    catch
        i
        disp(['failed creating SigT: ' txt{ID(i),1}]);
    end
end