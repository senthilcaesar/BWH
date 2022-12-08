#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# -------------------------------------------------------------------------------
# expects 3 args:
#  1. script folder (containing this file)
#  2. NAP resource folder (currently unused)
#  3. NAP output for this EDF
# -------------------------------------------------------------------------------

if ( length(args) != 3 )
 stop( "usage: coda2.R <nap-script-folder> <nap-resources-folder> <nap-output-folder>" )

# ensure trailing folder delimiter

nap.dir             <- paste( args[1] , "/" , sep="" )
nap.resources.dir   <- paste( args[2] , "/" , sep="" )
nap.output.dir      <- paste( args[3] , "/" , sep="" )

# -------------------------------------------------------------------------------
#
# Attach library dependencies
#
# -------------------------------------------------------------------------------

nap.paths <- .libPaths()

# add any nonstandard library paths here: 
nap.paths <- c( "/data/nsrr/lib/" , nap.paths )

.libPaths( nap.paths )

suppressMessages( library( data.table ) )


# -------------------------------------------------------------------------------
#
# Helper functions
#
# -------------------------------------------------------------------------------

source( paste( nap.dir , "codaf.R" , sep="" ) )


# ------------------------------------------------------------------------------------------------------
#
# By convention, we have special filename encoding to recognize particular outputs, e.g. HEADER-CH.txt
#
# nb. if a vector of transpose values is given, currently this transposes based on alphabetic table order
# e.g. for SOAP: there is the base ("") table only (based on alphabetic sort order) , so "_E", "_SS" then ""
# 
# ------------------------------------------------------------------------------------------------------


try( fextract( nap.output.dir , "luna_core_CANONICAL" , "Canonical signals" , tables = "_CS" , transpose = F ) )

try( fextract( nap.output.dir , "luna_core_FLIP" , "EEG polarity flips" , tables = c( "_CH" , "_CH_METHOD") , transpose = F ) )

try( fextract( nap.output.dir , "luna_spec_PSD" , "Power spectral density" , tables = c( "_B_CH_SS-N2" , "_B_CH_SS-N3" ) ) )

try( fextract( nap.output.dir , "luna_spso_SPINDLES" , "Spindles/SO" , tables = c( "_F_CH_SS-N2" , "_F_CH_SS-N23" , "_CH_SS-N2" , "_CH_SS-N23"  ) , transpose = T ) )

try( fextract( nap.output.dir , "luna_macro_HYPNO" , "NREM cycles" , tables = c( "_C" ) ,  transpose = T ) )

# SOAP
try( fextract( nap.output.dir , "luna_suds_SOAP" , "SOAP" ,
               tables = c( "_E" , "_NSS_PRED_OBS" , "_SS" , "" ) ,
               transpose = c(F,F,F,T) ) )

# POPS
try( fextract( nap.output.dir , "luna_suds_POPS" , "POPS" ,
               tables = c( "_E" , "_PRED_OBS" , "_SS" , "" ) ,
	       transpose = c(F,F,F,T) ) )

# Respiratory Endotype Summary
try( fextract( nap.output.dir , "resp_endotype_summary" , "Endotype Summary" ,
               tables = c( "_SS_ALL" , "_SS_NREM" ) ,
               transpose = c(F,F) ) )

# Respiratory Endotype Flow QC
try( fextract( nap.output.dir , "resp_qc" , "FlowQC" ,
               tables = c( "_flow" ) ,
               transpose = F ) )

# ------------------------------------------------------------------------------------------------------
#
# SIGSTATS and STATS processing for Signals View (on the HARMONIZED EDFs) :
# this is repeated in coda1.R for the Original EDFs 
#
# ------------------------------------------------------------------------------------------------------

# Hjorth
sigstats.filename <- paste( nap.output.dir , "luna_stats_SIGSTATS_E_CH_HARM-1.txt" , sep="" )

# Means
stats.filename <- paste( nap.output.dir , "luna_stats_STATS_E_CH_HARM-1.txt" , sep="" )   

df1 <- df2 <- data.frame()

if ( file.exists( sigstats.filename ) ) {
df1 <- nap.read.table( sigstats.filename )
df1 <- df1[ , c( "ID" , "CH" , "E" , "H1" , "H2" ) ]
df1$H1 <- log( df1$H1 )
names(df1) <- c("ID","CH","E","S1","S2")
df1$S1[ ! is.finite( df1$S1 ) ] <- NA
df1$S2[ ! is.finite( df1$S2 ) ] <- 0
}

if ( file.exists( stats.filename ) ) {
df2 <- nap.read.table( stats.filename )
df2 <- df2[ , c( "ID" , "CH" , "E" , "MEAN" ) ]
df2$S2 <- NA
names(df2) <- c("ID","CH","E","S1","S2")
df2$S1[ ! is.finite( df2$S1 ) ] <- NA
}

# merge
if ( dim(df1)[1] > 0 ) {
 df <- df1
 if ( dim(df2)[1] > 0 ) df <- rbind( df , df2 )
} else {
 df <- df2
}
# save
if ( dim(df)[1]>0) 
 saveit( df , "harm.sigstats" , file= paste( nap.output.dir , "nap.harm.sigstats.RData" , sep="" ) )



# ------------------------------------------------------------------------------------------------------
#
# PSD spectra
#
# ------------------------------------------------------------------------------------------------------

try( {
 psd <- list()
 sss <- "N2"
 psd[[ "N2" ]] <- nap.read.table( paste( nap.output.dir , "luna_spec_PSD_F_CH_SS-N2.txt" , sep="/" ) )
 if ( file.exists( paste( nap.output.dir , "luna_spec_PSD_F_CH_SS-N3.txt" , sep="/" ) ) )
 {
  sss <- c( sss , "N3" )  
  psd[[ "N3" ]] <- nap.read.table( paste( nap.output.dir , "luna_spec_PSD_F_CH_SS-N3.txt" , sep="/" ) ) 
 }
 # already log-scaled, etc 
 psd.spec <- list( desc = "Welch EEG PSD (QC+ epochs)" )
 for (ss in sss )
 {
  dt <- psd[[ ss ]]
  chs <- unique( dt$CH )
  for (ch in chs )
  { 
   png( file = paste( nap.output.dir , paste( "psd-" , ch , "-" , ss , ".png" , sep="" ) , sep="/" ) , res=100 , width=500  , height = 300 )
   x <- dt$F[ dt$CH == ch ]
   y <- dt$PSD[ dt$CH == ch ]
   nx <- length(unique(x))
   ny <- length(unique(y))
   plot( x , y , lwd=2 , type="l" , xlab="Frequency (Hz)" , ylab="Absolute log(power)" )
   dev.off()
   psd.spec[[ paste( ch , ss ) ]] <- list( desc = paste( ch, ss ) , figure = paste( "psd-" , ch , "-" , ss , ".png" , sep="" )  )
 }} # next channel / SS
 if ( length( psd.spec ) > 0 ) 
  saveit( psd.spec , str = "psd" , file = paste( nap.output.dir , "psd-fig.RData" , sep="" ) )
}) # end of PSD block



# ------------------------------------------------------------------------------------------------------
#
# MTM segment track
#
# ------------------------------------------------------------------------------------------------------

# based on nap1.sh command;
# luna    -s 'MTM tw=5 segment-sec=5 segment-inc=1 epoch max=45 sig=csC4|csC3|csF3|csF4 dB'
#    time half-bandwidth (nw) = 5
#    number of tapers         = 9
#    spectral resolution      = 2Hz
#    segment duration         = 5s
#    segment step             = 1s
#    FFT size                 = 512
#    number of segments       = 40406
#    adjustment               = none

if ( file.exists( paste( nap.output.dir , "luna_spec_MTM_F_CH_SEG_SPSD-1.txt.gz" , sep="" ) ) )
{

# should be for a single canonical channel only
d <- read.table( paste( nap.output.dir , "luna_spec_MTM_F_CH_SEG_SPSD-1.txt.gz" , sep="" ) , header=T, stringsAsFactors=F ) 

# helper functions
fma <- function(x,o=30) {
 pn <- floor(o/2)
 xx <- c( rep(x[1],pn) , x , rep(x[length(x)],pn) )
 xx <- as.numeric( filter(xx,rep(1,o)/o) )
 xx[ (pn+1):(pn+length(x)) ] 
}

fwin <- function(x,p=0.10) {
 q <- quantile( x , c(p,1-p) , na.rm=T)
 x[ x< q[1] ] <- q[1] ;  x[ x> q[2] ] <- q[2]
 x
}

fstd <- function(x) { ( x-min(x,na.rm=T) ) / ( max(x,na.rm=T) - min(x,na.rm=T) )  }

fpad <- function(x,p) { c( rep(NA,p) , x , rep(NA,p) ) }

fdb <- function(x) { 10*log10(x) }

fbin <- function(x,nb) { as.numeric( tapply( x , floor( ( (1:length(x)) - 1) /  nb ) , mean , na.rm=T ) ) } 

# step 1: 
# sum as bands; log-scale; pad to full # of seconds
# nb. the 5-second MTM window will drop psec=2 seconds each side, so scale w/ fpad()
# (if MTM has segment-sec=5 , that is ) 

psec=2
slow   <- fpad( fdb( tapply( d$MTM[ d$F < 1 ] , d$SEG[ d$F < 1 ] , sum ) ) , psec  )
delta  <- fpad( fdb( tapply( d$MTM[ d$F >= 1  & d$F < 4 ]  , d$SEG[ d$F >= 1  & d$F < 4 ] , sum ) ) , psec )
theta  <- fpad( fdb( tapply( d$MTM[ d$F >= 4  & d$F < 8 ]  , d$SEG[ d$F >= 4  & d$F < 8 ] , sum ) ) , psec )
alpha  <- fpad( fdb( tapply( d$MTM[ d$F >= 8  & d$F < 11 ] , d$SEG[ d$F >= 8  & d$F < 11 ] , sum ) ) , psec )
sigma  <- fpad( fdb( tapply( d$MTM[ d$F >= 11 & d$F < 15 ] , d$SEG[ d$F >= 11 & d$F < 15 ] , sum ) ) , psec )
beta   <- fpad( fdb( tapply( d$MTM[ d$F >= 15 & d$F < 30 ] , d$SEG[ d$F >= 15 & d$F < 30 ] , sum ) ) , psec )
gamma  <- fpad( fdb( tapply( d$MTM[ d$F >= 30 & d$F < 45 ] , d$SEG[ d$F >= 30 & d$F < 45 ] , sum ) ) , psec )

# step 2:
#  winsorize at 5%,
#  take mwsec=5 second moving average
#  bin into: 1s, 5s, 10s, 30s
#  then winsorize each at 1%
#  then 0..1 scale

# typical study: say 10 hours ( 36000 seconds)
#   d1    up to 10 minutes  --> 1 sec resol ( up to 600 datapoints;  one epoch = 30 data points) 
#   d5    up to 1 hour      --> 5 sec resol ( up to 720 datapoints)
#   d10   up to 2 hours     --> 10 sec resol ( up to 720 datapoints)
#   d30   over 2 hours      --> 30 sec resol ( e.g. 1200 data points for 10 hour study)
#                                                   but typically this should work well
w1    <- 0.05
w2    <- 0.01
mw    <- 5

# 5 bins MA :   5 seconds for 1 sec bins
#               25 seconds for 5 sec bins
#               50 seconds for 10 sec bins
#               150 seconds for 30 second bins

slow1  <- fstd( fwin( fma( fwin(slow,w1) , mw ) , w2 ) )  
delta1 <- fstd( fwin( fma( fwin(delta,w1) , mw ) , w2 ) ) 
theta1 <- fstd( fwin( fma( fwin(theta,w1) , mw ) , w2 ) )
alpha1 <- fstd( fwin( fma( fwin(alpha,w1) , mw ) , w2 ) )
sigma1 <- fstd( fwin( fma( fwin(sigma,w1) , mw ) , w2 ) )
beta1  <- fstd( fwin( fma( fwin(beta,w1) , mw ) , w2 ) )
gamma1 <- fstd( fwin( fma( fwin(gamma,w1) , mw ) , w2 ) )

slow5  <- fstd( fwin( fma( fbin( fwin(slow,w1) , 5), mw ) , w2 ) )  
delta5 <- fstd( fwin( fma( fbin( fwin(delta,w1) , 5 ), mw ) , w2 ) ) 
theta5 <- fstd( fwin( fma( fbin( fwin(theta,w1) , 5 ), mw ) , w2 ) )
alpha5 <- fstd( fwin( fma( fbin( fwin(alpha,w1) , 5), mw ) , w2 ) )
sigma5 <- fstd( fwin( fma( fbin( fwin(sigma,w1) , 5), mw ) , w2 ) )
beta5  <- fstd( fwin( fma( fbin( fwin(beta,w1) , 5), mw ) , w2 ) )
gamma5 <- fstd( fwin( fma( fbin( fwin(gamma,w1) , 5), mw ) , w2 ) )

slow10  <- fstd( fwin( fma( fbin( fwin(slow,w1) , 10), mw ) , w2 ) )  
delta10 <- fstd( fwin( fma( fbin( fwin(delta,w1) , 10 ), mw ) , w2 ) ) 
theta10 <- fstd( fwin( fma( fbin( fwin(theta,w1) , 10 ), mw ) , w2 ) )
alpha10 <- fstd( fwin( fma( fbin( fwin(alpha,w1) , 10), mw ) , w2 ) )
sigma10 <- fstd( fwin( fma( fbin( fwin(sigma,w1) , 10), mw ) , w2 ) )
beta10  <- fstd( fwin( fma( fbin( fwin(beta,w1) , 10), mw ) , w2 ) )
gamma10 <- fstd( fwin( fma( fbin( fwin(gamma,w1) , 10), mw ) , w2 ) )

slow30  <- fstd( fwin( fma( fbin( fwin(slow,w1) , 30), mw ) , w2 ) )  
delta30 <- fstd( fwin( fma( fbin( fwin(delta,w1) , 30 ), mw ) , w2 ) ) 
theta30 <- fstd( fwin( fma( fbin( fwin(theta,w1) , 30 ), mw ) , w2 ) )
alpha30 <- fstd( fwin( fma( fbin( fwin(alpha,w1) , 30), mw ) , w2 ) )
sigma30 <- fstd( fwin( fma( fbin( fwin(sigma,w1) , 30), mw ) , w2 ) )
beta30  <- fstd( fwin( fma( fbin( fwin(beta,w1) , 30), mw ) , w2 ) )
gamma30 <- fstd( fwin( fma( fbin( fwin(gamma,w1) , 30), mw ) , w2 ) )

spsd <- list(
 d1=data.frame(SLOW=slow1,DELTA=delta1,THETA=theta1,ALPHA=alpha1,SIGMA=sigma1,BETA=beta1,GAMMA=gamma1),
 d5=data.frame(SLOW=slow5,DELTA=delta5,THETA=theta5,ALPHA=alpha5,SIGMA=sigma5,BETA=beta5,GAMMA=gamma5),
 d10=data.frame(SLOW=slow10,DELTA=delta10,THETA=theta10,ALPHA=alpha10,SIGMA=sigma10,BETA=beta10,GAMMA=gamma10),
 d30=data.frame(SLOW=slow30,DELTA=delta30,THETA=theta30,ALPHA=alpha30,SIGMA=sigma30,BETA=beta30,GAMMA=gamma30) )

spsd$d1[ is.na(spsd$d1) ] <- 0
spsd$d5[ is.na(spsd$d5) ] <- 0
spsd$d10[ is.na(spsd$d10) ] <- 0
spsd$d30[ is.na(spsd$d30) ] <- 0

saveit( spsd , "spsd" , file= paste( nap.output.dir , "nap.spsd.RData" , sep="" ) )

}



