function TableOut = SignaltoBreathStats(S,BB_i_start,BB_i_mid,BB_i_end,siglabel,keepcols)

%before 12/6 was error in 95th centile (was actually 5th centile)

for i=1:length(BB_i_start)
    I = BB_i_start(i):BB_i_end(i);
    Smax(i,1) = max(S(I));
    Smin(i,1) = min(S(I));
    Smean(i,1) = nanmean(S(I));
    Smedian(i,1) = nanmedian(S(I));
    S5thcentile(i,1) = prctile(S(I),5);
    S10thcentile(i,1) = prctile(S(I),10);
    S90thcentile(i,1) = prctile(S(I),90);
    S95thcentile(i,1) = prctile(S(I),95);
    I = BB_i_start(i):BB_i_mid(i);
    SmaxI(i,1) = max(S(I));
    SminI(i,1) = min(S(I));
    SmeanI(i,1) = nanmean(S(I));
    SmedianI(i,1) = nanmedian(S(I));
    S5thcentileI(i,1) = prctile(S(I),5);
    S10thcentileI(i,1) = prctile(S(I),10);
    S90thcentileI(i,1) = prctile(S(I),90);
    S95thcentileI(i,1) = prctile(S(I),95);
    I = BB_i_mid(i):BB_i_end(i);
    SmaxE(i,1) = max(S(I));
    SminE(i,1) = min(S(I));
    SmeanE(i,1) = nanmean(S(I));
    SmedianE(i,1) = nanmedian(S(I));
    S5thcentileE(i,1) = prctile(S(I),5);
    S10thcentileE(i,1) = prctile(S(I),10);
    S90thcentileE(i,1) = prctile(S(I),90);
    S95thcentileE(i,1) = prctile(S(I),95);
end

Table1 = table(Smax,Smin,Smean,Smedian,S5thcentile,S10thcentile,S90thcentile,S95thcentile);
TableI = table(SmaxI,SminI,SmeanI,SmedianI,S5thcentileI,S10thcentileI,S90thcentileI,S95thcentileI);
TableE = table(SmaxE,SminE,SmeanE,SmedianE,S5thcentileE,S10thcentileE,S90thcentileE,S95thcentileE);

TableOut = [Table1(:,keepcols) TableI(:,keepcols) TableE(:,keepcols)];

temp = TableOut.Properties.VariableNames;
for i=1:length(temp)
    temp{i}=[siglabel temp{i}(2:end)];
end
TableOut.Properties.VariableNames = temp;





