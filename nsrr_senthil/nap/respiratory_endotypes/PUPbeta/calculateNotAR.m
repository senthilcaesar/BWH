function nota = calculateNotAR(AR,winID,NExcludepostAR)
    %nota1 = calculateNotAR(BreathDataTable3.AR3,BreathDataTable3.Time0,2);
    nota=1-AR;
    
    Istartnewwin = [1;(find(diff(winID)~=0)+1)];
    for i=1:(NExcludepostAR-1)
        Istartnewwin = [Istartnewwin; Istartnewwin+1];
    end
    Istartnewwin = sort(Istartnewwin);
    Istartnewwin(Istartnewwin>length(AR))=[];
    nota(Istartnewwin)=0;
    
    Iendarousal = [(find(diff(AR)==-1)+1)];
    for i=1:(NExcludepostAR-1)
        Iendarousal = [Iendarousal; Iendarousal+1];
    end
    Iendarousal = sort(Iendarousal);
    Iendarousal(Iendarousal>length(AR))=[];
    
    nota(Iendarousal)=0;
    
%     notaT = table;
%     notaT.AR=AR;
%     notaT.winID=winID;
%     notaT.nota=nota;
    
    
    