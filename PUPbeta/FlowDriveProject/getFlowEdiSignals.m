function [Flow, Edi] = getFlowEdiSignals(LocalSignals, n, win)
    Flow = LocalSignals{n}{win}(:,1);
    Edi = LocalSignals{n}{win}(:,2);
end
