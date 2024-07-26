function HRresponse = HRresponseFromSignalT(SignalsT,LocalEvtsT,HRsignal,ArMode,plotaxish)
%HRresponseFromSignalT(SignalsT,Evts.RespT,'HRfilt',0,gca)

        
        if ArMode
            LocalEvtsT.EventEnd = LocalEvtsT.EventStart;
            LocalEvtsT.EventStart = LocalEvtsT.EventStart-10;     
        end
        [Boxes,Ensembles]=EventAnalysisRun(SignalsT,[],LocalEvtsT);
        
        HR = getfield(Boxes,HRsignal);
        End = repmat(101,size(HR,1),1);
        Start = End - round(LocalEvtsT.EventDuration);
        Time = Ensembles.Time';
        
        HRresponse = HRdelta(HR,Start,End,Time);
        
        if exist('plotaxish') && ~isempty(plotaxish)
            axes(plotaxish);
            plot(Ensembles{:,1},Ensembles{:,2});
            set(plotaxish,'box','off','tickdir','out');
            xlim([-45 45]);
        end
        