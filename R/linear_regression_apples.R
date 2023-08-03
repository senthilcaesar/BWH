library(luna)
library(data.table)

# we consider a positive control analysis, relating
# demographic variables (age and sex) to spectral power during N2 sleep.
# What is the relationship between N2 spectral power and demographic variables (Age, Sex and Race) ?
# The values of N2 power is a function of (Age, Sex and Race)
pw <- read.table("/Users/sq566/Desktop/apples/summ.psd",header=T,stringsAsFactors=F)
pw <- setDF( dcast( setDT( pw ) , ID ~ CH + F  , value.var = names(pw)[-(1:3)] ) )
d <- read.table("/Users/sq566/Desktop/apples/demo.txt",header=T,stringsAsFactors=F, sep='\t')
d$ID <- gsub('id_','apples-',d$ID)
pw <- merge(pw, d , by="ID" )

# get all PSD variables ( power 0.25 to 64 Hz for C3, then C4)
vars <- names(pw)[2:511] 

f1 <- function(x) {
  as.numeric(coef(summary(lm(outliers(x) ~ I(age/10) + sex + as.factor(site) , data = pw)))[2:3,c(1,4)])
}

f2 <- function(x) {
  regression <-  coef(summary(lm(outliers(x) ~ I(age/10) + sex + as.factor(site) , data = pw)))
  print(regression)
}

k <- data.frame(t(sapply(pw[,vars] , f1)))

# Age Estimate, Sex Estimate, Age P-value, Sex P-value
names(k) <- c("e.age","e.male","p.age","p.male") 
k$z.age <- sign( k$e.age ) * -log10( k$p.age )
k$z.male <- sign( k$e.male ) * -log10( k$p.male )
k$F <- rep( seq( 0.5 , 64 , 0.25 ) ,2 ) 
k$CH <- ifelse( grepl( "C3" , rownames(k) ), "C3","C4" )

# Plot
fidx <- unique( k$F ) 
png(file="/Users/sq566/Desktop/apples/apples-summ-n2-power-demo.png", height=500, width=1000, res=100 )
par( mfrow=c(1,2) )
plot( fidx, k$z.age[k$CH == "C3"] , lwd=2, col = rgb(100,0,100,100,max=255), ylim=c(-25,5),
      main = "N2 PSD ~ age" , ylab="Signed -log10(p)" , xlab="Frequency (Hz)" , type="l" )
lines( fidx, k$z.age[ k$CH == "C4" ] , lwd=2, col = rgb(0,100,100,100,max=255) )
abline(h=0); abline(h=c(3,-3),lty=2,col="gray")

# sex diffs
plot( fidx, k$z.male[k$CH == "C3"] , lwd=2, col = rgb(100,0,100,100,max=255), ylim=c(-25,5),
      main = "N2 PSD ~ male" , ylab="Signed -log10(p)" , xlab="Frequency (Hz)" , type="l" )
lines( fidx, k$z.male[ k$CH == "C4" ] , lwd=2, col = rgb(0,100,100,100,max=255) )
abline(h=0); abline(h=c(3,-3),lty=2,col="gray")
dev.off()


# Extract row 2 and 3, Column 1 and 4
# regression[2:3,c(1,4)]
"                   Estimate   Std Error     t value     Pr(>|t|)
(Intercept)       -74.80190952 0.10543575 -709.454909 0.000000e+00
I(age/10)          -0.04673366 0.01730100   -2.701212 7.019757e-03
sexM               -0.11889224 0.04416064   -2.692267 7.209507e-03
as.factor(site)SL   0.99433015 0.07096236   14.012077 5.009579e-41
as.factor(site)SM   1.02156677 0.08263421   12.362516 7.069294e-33
as.factor(site)SU   0.67933345 0.06245144   10.877787 3.433029e-26
as.factor(site)UA   0.74675191 0.06137771   12.166499 5.872305e-32"






