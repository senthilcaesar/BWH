library(data.table)
sp <- read.table("/Users/sq566/Desktop/dreem/summ.soap",header=T,stringsAsFactors=F)
sp <- setDF( dcast( setDT( sp ) , ID ~ SCH , value.var = names(sp)[-(1:2)] ) )

my_vector <- 1:length(sp$ID)

png(file="/Users/sq566/Desktop/dreem/03-summ-soap.png", height=1600, width=3500, res=300 )
par(mfcol=c(1,2))
plot( sp$K3_C3 , pch=20 , ylim=c(0, 1), ylab="SOAP kappa C3" , xlab="Individual", col = "black", cex = 0.8, bg="navy"  )
#label_indices <- which(sp$K3_C3 < 0.5)
#text(my_vector[label_indices], sp$K3_C3[label_indices], labels=sp$ID[label_indices], pos=3, col=rgb(1, 0.5, 0.5), cex=0.8)

plot( sp$K3_C4 , pch=20 , ylim=c(0, 1), ylab="SOAP kappa C4" , xlab="Individual", col = "black", cex = 0.8, bg="navy"  )
#label_indices <- which(sp$K3_C4 < 0.4)
#text(my_vector[label_indices], sp$K3_C4[label_indices], labels=sp$ID[label_indices], pos=3, col=rgb(1, 0.5, 0.5), cex=0.8)

mtext("Dreem 03", side=3, outer=TRUE, line=-2, cex=1.5)
dev.off()

# Get the IDs
sp2 <- sp[ sp$K3_C4 < 0.4, ]
sp2[,c("ID", "K3_C4")]
sp3 <- sp[ sp$K3_C3 < 0.4, ]
sp3[,c("ID", "K3_C3")]

