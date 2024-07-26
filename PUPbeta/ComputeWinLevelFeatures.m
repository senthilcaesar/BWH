function [WinDataTable,WinSnoreTable,PSDwin,EnvelpWin,EnvelpWin_n] = ComputeSnoreFeaturesNew(SnoreStruct,DBthresh)

SnoreCellArray = SnoreStruct.SnoreCellArray(1,:);
WinSnoreTable = [];
EnvelpWin = [];
EnvelpWin_n = [];
for ii = 1:length(SnoreCellArray)
    SnoreArray = SnoreStruct.SnoreCellArray(ii);
    EnvelpWin = SnoreStruct.lpcCellArray(1,ii);
    EnvelpArray_n = SnoreStruct.lpcCellArray(2,ii);
    if ii == 1
        WinSnoreTable = SnoreArray(:,2:end);
        EnvelpWin = EnvelpArray;
        EnvelpWin_n = EnvelpArray_n;
    else
        WinSnoreTable = vertcat(WinSnoreTable,SnoreArray(:,2:end));
        EnvelpWin = vertcat(EnvelpWin,EnvelpArray);
        EnvelpWin_n = vertcat(EnvelpWin_n,EnvelpArray_n);
    end
end