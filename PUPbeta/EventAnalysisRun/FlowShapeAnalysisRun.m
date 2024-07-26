function Ftrs = FlowShapeAnalysisRun(BrDataTable,BrFLDataTable)

% Calculate mean value of flow shape features for hypopneas

BrFLDataTable_ = BrFLDataTable(BrDataTable.FDuplicated==0 &...
    BrDataTable.Etype==4,:);

FtrsTable = varfun(@nanmean,BrFLDataTable_);
Ftrs = table2struct(FtrsTable);

