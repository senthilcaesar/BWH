#library(luna)

#k <- ldb("/Users/sq566/Desktop/mesa-5/out.db")
#ss <- lx(k , "STAGE" , "E" )$STAGE 

#d <- lx(k , "CORREL" , "CH1" , "CH2" , "E") 
d <- read.csv("/Users/sq566/Desktop/mesa-5/corr.csv", sep = "\t")

average_d <- aggregate(R ~ ID, data = d, FUN = mean, na.rm = TRUE)

average_d$ID_numeric <- as.numeric(as.factor(average_d$ID))

#png(file="/Users/sq566/Desktop/mesa-5/plots/mesa5-corr.png", width = 2400, height = 1200, res = 300)

plot(average_d$ID_numeric, average_d$R, main = "N2 Sleep (HR - DHR)",
     xlab = "Individuals", ylab = "Average Correlation", pch = 19, col = "blue", ylim=c(-1,1))

abline( h=0 , lty=2 )
