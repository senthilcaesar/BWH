h1 <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.hypno1",header=T,stringsAsFactors=F)
h2 <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.hypno2",header=T,stringsAsFactors=F)
sp <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.soap",header=T,stringsAsFactors=F)
st <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.stats",header=T,stringsAsFactors=F)
hj <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.hjorth",header=T,stringsAsFactors=F)
pw <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.psd",header=T,stringsAsFactors=F)

library(data.table)
h2 <- setDF( dcast( setDT( h2 ) , ID ~ SS , value.var = names(h2)[-(1:2)] ) )
sp <- setDF( dcast( setDT( sp ) , ID ~ SCH , value.var = names(sp)[-(1:2)] ) )
st <- setDF( dcast( setDT( st ) , ID ~ CH , value.var = names(st)[-(1:2)] ) )
hj <- setDF( dcast( setDT( hj ) , ID ~ CH , value.var = names(hj)[-(1:2)] ) )
pw <- setDF( dcast( setDT( pw ) , ID ~ CH + F  , value.var = names(pw)[-(1:3)] ) )

d <- read.table("/Users/sq566/Desktop/mesa-5/demo.txt",header=T,stringsAsFactors=F, sep='\t')


# merge
h1 <- merge(h1, d , by="ID" )
h2 <- merge(h2, d , by="ID" )
sp <- merge(sp, d , by="ID" )
hj <- merge(hj, d , by="ID" )
st <- merge(st, d , by="ID" )
pw <- merge(pw, d , by="ID" )

png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-macro1.png", height=400, width=800, res=100 )
par(mfcol=c(1,2))
plot( h1$TST , pch=21 , ylab="TST (in Minutes)" , col = "black", cex = 0.6, bg="navy" )
plot( h1$TRT , pch=21 , ylab="TRT (in Minutes)" , col = "black", cex = 0.6, bg="navy" )
dev.off()

png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-macro2.png", height=400, width=800, res=100 )
par(mfcol=c(1,4))
plot( h2$PCT_N1 , pch=21 , ylab="% N1", col = "black", cex = 0.6, bg="navy"  )
plot( h2$PCT_N2 , pch=21 , ylab="% N2", col = "black", cex = 0.6, bg="navy"  )
plot( h2$PCT_N3 , pch=21 , ylab="% N3", col = "black", cex = 0.6, bg="navy"  )
plot( h2$PCT_R  , pch=21 , ylab="% R", col = "black", cex = 0.6, bg="navy"  )
dev.off()

# Get a list of channel labels (n=17)
chs <- gsub( "MEAN_","", names(st)[grep( "MEAN_" , names(st) ) ] )
png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-means.png", height=800, width=1000, res=100 )
par(mfcol=c(4,4) , mar=c(1,4,1,1) )
for (ch  in chs ) {
  # To remove outliers from the data, we can use several methods depending on how we define an outlier. 
  # A common approach is to use the interquartile range (IQR). Here's how we can do it:

  # Calculate the IQR: This is the difference between the 5th and 95th percentiles. 
  # Values outside 1.5 times the IQR below the 5th percentile or above the 95th percentile are 
  # typically considered outliers.

  # Filter the Data: Exclude data points that are outside these bounds.

  column_name <- paste("MEAN", ch, sep="_")

  # Calculate Q1 and Q3 and then IQR
  Q1 <- quantile(st[, column_name], 0.5, na.rm = TRUE)
  Q3 <- quantile(st[, column_name], 0.95, na.rm = TRUE)
  IQR <- Q3 - Q1
  # Define bounds
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  filtered_data <- st[!is.na(st[[column_name]]) & st[[column_name]] >= lower_bound & st[[column_name]] <= upper_bound, ]
  
  plot( filtered_data[, column_name], pch=20, cex=0.6 , ylab=ch)
}
dev.off()

library(luna)
png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-means2.png", height=800, width=1000, res=100 )
par(mfcol=c(4,4) , mar=c(1,4,1,1))
for (ch  in chs ) {
  mn <- outliers(st[,paste("MEAN",ch,sep="_")])
  p5 <- st[,paste("P05",ch,sep="_")]
  p95 <- st[,paste("P95",ch,sep="_")]
  p5 <- p5[ ! is.na( mn ) ]
  p95 <- p95[ ! is.na( mn ) ]
  mn <- mn[ ! is.na( mn ) ]
  ylim <- range( c(p5,p95) , na.rm=T )
  plot( mn , ylim = ylim ,col="black" , pch=20, cex=0.6 , ylab=ch)
  points( p5, pch=20, cex=0.6 , col="blue" )
  points( p95, pch=20, cex=0.6 , col="red" )
}
dev.off()

png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-soap.png", height=500, width=800, res=100 )
par(mfcol=c(1,1))
plot( sp$K3_C4 , pch=20 , ylab="SOAP kappa (C4-M1)" , xlab="Individual", col = "black", ylim=c(0, 1), cex = 0.6, bg="navy"  )
dev.off()


library(data.table)
sp <- read.table("/Users/sq566/Desktop/mesa-5/res/summ.soap",header=T,stringsAsFactors=F)
sp <- setDF( dcast( setDT( sp ) , ID ~ SCH , value.var = names(sp)[-(1:2)] ) )
d <- read.table("/Users/sq566/Desktop/mesa-5/demo.txt",header=T,stringsAsFactors=F, sep='\t')
sp <- merge(sp, d , by="ID" )
sp2 <- sp[ sp$K3_C4 < 0.4, ]
sp2[,c("ID", "K3_C4")]

freq <- rep( seq( 0.5 , 64 , 0.25 ) ,2 )
fidx <- unique( freq )
vars.c4 <- names(pw)[2:256]
ylim <- c(-30, 60)
ids <- unique( pw$ID )
sites <- unique( pw$maternalrace )

png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-psd1.png", width = 800, height = 600, res = 100 )
par(mfcol=c(1,1))
plot( fidx , fidx , type="n" , ylim = ylim , xlab="Frequency (Hz)", ylab="Power" , main="C4_M1" )
for (id in ids ) lines( fidx, as.numeric( pw[pw$ID == id,vars.c4] ) , col = rgb(0,100,100,25,max=255), lwd=1 )
dev.off()

freq <- rep( seq( 0.5 , 64 , 0.25 ) ,2 )
fidx <- unique( freq )
vars.c4 <- names(pw)[2:100]
ylim <- c(-30, 60)
fidx <- fidx[ fidx <= 25 ]
ids <- unique( pw$ID )

png(file="/Users/sq566/Desktop/mesa-5/plots/pats-summ-psd2.png", width = 800, height = 700, res = 100 )
par(mfcol=c(1,1))
plot( fidx , fidx , type="n" , ylim = ylim , xlab="Frequency (Hz)", ylab="Power" , main="C4_M1" )
for (id in ids ) lines( fidx, as.numeric( pw[pw$ID == id,vars.c4] ) , col = rgb(0,100,100,25,max=255), lwd=1 )

dev.off()


# Plot PSD average by site
png(file="/Users/sq566/Desktop/mesa-5/plots/mesa5-psd-mean-site.png", width = 800, height = 700, res = 100 )
frqs <- seq(0.5, 60, 0.5)
col <- c("brown", "green", "blue", "pink", "red", "purple")
ch = "C4_M1"
plot(range(frqs), c(-50, 50), type="n", xlab = "Hz", ylab = "PSD", main = "Mesa5 PSD mean (C4_M1)")
for (i in 3:8) {
  vars <- paste0(ch, "_", frqs)
  M <- colMeans(pw[pw$site == i, vars])
  lines(frqs, M, col=col[i-2], lwd=1)
}
legend("topright", legend = paste("site",3:8), lwd=2, col=col)
dev.off()


# Plot PSD by site and individuals
png(file="/Users/sq566/Desktop/mesa-5/plots/mesa5-psd-site.png", width = 1200, height = 800, res = 150 )
frqs <- seq(0.5, 60, 0.5)
col <- c("brown", "green", "blue", "pink", "red", "purple")
ch = "C4_M1"
par(mfcol=c(2,3), mar=c(3,3,3,1))  # Adjust margins as needed
sites <- 3:8

# Looping through each site to plot data
for (site in sites) {
  
  plot(range(frqs), c(-50, 50), type="n", xlab = "Hz", ylab = "PSD")
  vars <- paste0(ch, "_", frqs) 
  # Plot data for each site
  
  for (id in unique(pw$ID[pw$site == site])) {
    data_to_plot <- as.numeric(pw[pw$ID == id & pw$site == site, vars])
    lines(frqs, data_to_plot, col=col[site-2], lwd=1)
  }
  average_series <- colMeans(pw[pw$site == site, vars])
  lines(frqs, average_series, col = "black", lwd = 2)
  num_individuals <- length(unique(pw$ID[pw$site == site]))
  legend("topright", legend = paste("Site", site, "-", num_individuals, "IDs"), col = col[site-2], lwd = 1, bty = "n")
}
mtext("Mesa-5 PSD (C4_M1)", outer = TRUE, cex = 1, line = -2)
# Closing the device to save the file
dev.off()






