function getEvtPRAddon(MrangeOverride)

% RUN StartHere.m first

global AMasterSpreadsheet

t_startGetData = clock;

%%

%[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26'); % only used incase of analyzed file

[~,txt] = xlsread(AMasterSpreadsheet,1,'AD4:AE5000');
M = size(txt,1);

ID=1:M;
if exist('MrangeOverride')
    ID = MrangeOverride(:);
end

success=zeros(length(ID),1);

for i=1:length(ID)
    try
        clear W Evts
        
        %Setup for getAddon mode
        W = load([txt{ID(i),2} txt{ID(i),1}]);
        disp(['processing: ' txt{ID(i),1}])
        
        if isfield(W, 'DataEventHypnog_Mat')
            SigT=array2table(W.DataEventHypnog_Mat);
            SigT.Properties.VariableNames = W.ChannelsList;
            % Field is there.  Remove it.
            W = rmfield(W, {'DataEventHypnog_Mat','ChannelsList'});
        else
            SigT=W.SigT;
        end
        
        
        Evts=W.Evts;
        
        Evts=getEvtPRMainFn(Evts,SigT);
        
        
        W.Evts=Evts; %overwrite updated.
        figure(23);
        save([txt{ID(i),2} txt{ID(i),1}],'-struct','W','-v7.3');
        success(i)=1;
        disp(['success:Pulse rate added to: ' txt{ID(i),1}]);
    catch
        i
        disp(['failed adding pulse rate to: ' txt{ID(i),1}]);
    end
end




