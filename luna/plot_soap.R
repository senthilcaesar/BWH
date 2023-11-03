library(data.table)

sp <- read.table("/Users/sq566/Desktop/saps/summ.soap",header=T,stringsAsFactors=F)
sp <- setDF( dcast( setDT( sp ) , ID ~ SCH , value.var = names(sp)[-(1:2)] ) )
sp$site <- substr(sp$ID, start = 1, stop = 4)

png(file="/Users/sq566/Desktop/saps/STAGES-SOAP.png", height=800, width=2600, res=100)

par(mfcol=c(1,2))
plot( sp$K3_C3 , col = as.factor(sp$site), pch=16 , ylab="SOAP kappa (C3)" , xlab="Individual", cex = 0.8)
legend("bottomleft", legend=levels(as.factor(sp$site)), col=1:14, pch=16)
title(main = "STAGES / C3 / SOAP", cex.main = 1.5)
plot( sp$K3_C4 , col = as.factor(sp$site), pch=16 , ylab="SOAP kappa (C4)" , xlab="Individual", cex = 0.8)
legend("bottomleft", legend=levels(as.factor(sp$site)), col=1:14, pch=16)
title(main = "STAGES / C4 / SOAP", cex.main = 1.5)

dev.off()
