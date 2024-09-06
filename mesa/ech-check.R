d <- read.table("/Users/sq566/Desktop/mesa-5/res/res.sec.gz", header=T, stringsAsFactors=F)
demo <- read.table("/Users/sq566/Desktop/mesa-5/demo.txt",header=T,stringsAsFactors=F, sep='\t')
d <- merge(d, demo , by="ID" )
ch="ECG"

png(file="/Users/sq566/Desktop/mesa-5/plots/mesa5-ECG.png", width = 2400, height = 1200, res = 300 )
plot(tapply( d$MD[ d$CH == ch ], d$SEC[ d$CH == ch ], mean ), type="l", main =ch, ylab="", ylim = c(-2,2), col = "red", lwd = 2)
ids <- unique(d$ID[ d$CH == ch ])
for (id   in ids) lines(d$MD[ d$CH == ch & d$ID == id ], col=rgb(100,100,100,100,max=255), lwd=0.2)
lines(tapply( d$MD[ d$CH == ch ], d$SEC[ d$CH == ch ], mean ), type="l", main =ch, ylab="", ylim = c(-2,2), col = "black", lwd = 1)
dev.off()


png(file="/Users/sq566/Desktop/mesa-5/plots/mesa5-ECG-site.png", width = 2400, height = 1200, res = 300 )

col <- c("brown", "green", "blue", "pink", "red", "purple")
ch = "ECG"
par(mfcol=c(2,3), mar=c(3,3,3,1))  # Adjust margins as needed
sites <- 3:8

# Looping through each site to plot data
for (site in sites) {
  
  plot(tapply(d$MD[d$CH == ch & d$site == site], d$SEC[d$CH == ch & d$site == site], mean), type="l", 
       xlab="", ylab="", ylim = c(-2,2), col = "black", lwd = 1)
  
  for (id in unique(d$ID[d$site == site])) {
    data_to_plot <- d$MD[d$ID == id & d$site == site]
    lines(data_to_plot, col=col[site-2], lwd=0.2)
  }
  
  lines(tapply(d$MD[d$CH == ch & d$site == site], d$SEC[d$CH == ch & d$site == site], mean), type="l", 
       xlab="", ylab="", ylim = c(-2,2), col = "black", lwd = 1)
  
  num_individuals <- length(unique(d$ID[d$site == site]))
  legend("topright", legend = paste("Site", site, "-", num_individuals, "IDs"), col = col[site-2], lwd = 1, bty = "n")
  
  trans_black = rgb(0, 0, 0, alpha = 0.5)
  abline(h = 0, col = trans_black, lty = "dashed", lwd = 1)
}

mtext("ECG", outer = TRUE, cex = 1, line = -2)
# Closing the device to save the file
dev.off()
