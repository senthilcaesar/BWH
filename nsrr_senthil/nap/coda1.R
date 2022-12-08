#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# -------------------------------------------------------------------------------
# expects 5 args:
#  1. script folder (containing this file)
#  2. NAP resource folder (currently unused)
#  3. NAP output for this EDF
#  4. Optional path for R libraries
#  5. 0/1 flag for whether to load MTM 
# -------------------------------------------------------------------------------

if ( length(args) != 5 )
 stop( "usage: coda1.R <nap-script-folder> <nap-resources-folder> <nap-output-folder> <nap-rlibs> <mtm=0/1>" )

# ensure trailing folder delimiter

nap.dir             <- paste( args[1] , "/" , sep="" )
nap.resources.dir   <- paste( args[2] , "/" , sep="" )
nap.output.dir      <- paste( args[3] , "/" , sep="" )
nap.rlibs           <- paste( args[4] , "/" , sep="" )
nap.domtm           <- args[5]

# -------------------------------------------------------------------------------
#
# Attach library dependencies
#
# -------------------------------------------------------------------------------

nap.paths <- .libPaths()

# add any nonstandard library paths here: 
if ( nap.rlibs != "" ) nap.paths <- c( "/data/nsrr/lib/" , nap.rlibs , nap.paths )

.libPaths( nap.paths )

suppressMessages( library( data.table ) )
suppressMessages( library( luna ) )


# -------------------------------------------------------------------------------
#
# Helper functions
#
# -------------------------------------------------------------------------------

source( paste( nap.dir , "codaf.R" , sep="" ) )


# ------------------------------------------------------------------------------------------------------
#
# Annotation mapping
#
# ------------------------------------------------------------------------------------------------------

annot.alias.filename <- paste( nap.output.dir , "luna_core_ALIASES_ANNOT.txt" , sep="" )
annot.map.filename <- paste( nap.output.dir , "luna_core_ANNOTS_ANNOT_INST_A.txt" , sep="" )

if ( file.exists( annot.alias.filename ) ) {
df <- nap.read.table( annot.alias.filename )
df <- df[, c(3,2) ] # ORIG , ANNOT
names(df) <- c( "ORIGINAL" , "ALIAS" ) 
if ( dim(df)[1]>0)
 saveit( df , "annot.alias" , file= paste( nap.output.dir , "nap.annot.alias.RData" , sep="" ) )
}

if ( file.exists( annot.map.filename ) ) {
df <- nap.read.table( annot.map.filename )
df <- df[, c(3,4,2) ] # CLASS, INST , MAPPED
names(df) <- c( "ANNOT" , "INST" , "MAPPED" )
if ( dim(df)[1]>0)
 saveit( df , "annot.map" , file= paste( nap.output.dir , "nap.annot.map.RData" , sep="" ) )
}


# ------------------------------------------------------------------------------------------------------
#
# Primary channel mapping (original --> NSRR harmonized set)
#
# ------------------------------------------------------------------------------------------------------

sigs.map1.filename <- paste( nap.output.dir , "luna_core_CANONICAL_CS_EDF-HARM.txt" , sep="" )
sigs.map2.filename <- paste( nap.output.dir , "luna_core_CANONICAL_CH_EDF-HARM.txt" , sep="" )

if ( file.exists( sigs.map1.filename ) ) {
df <- nap.read.table( sigs.map1.filename )
if ( dim(df)[1]>0) saveit( df , "sigs.harm.map1" , file= paste( nap.output.dir , "nap.sig.map1.RData" , sep="" ) )
}

if ( file.exists( sigs.map2.filename ) ) {
df <- nap.read.table( sigs.map2.filename )
if ( dim(df)[1]>0) saveit( df , "sigs.harm.map2" , file= paste( nap.output.dir , "nap.sig.map2.RData" , sep="" ) )
}


# ------------------------------------------------------------------------------------------------------
#
# Secondary channel mapping (harmonized --> canonical)
#
# ------------------------------------------------------------------------------------------------------

sigs.map1.filename <- paste( nap.output.dir , "luna_core_CANONICAL_CS_EDF-BASE.txt" , sep="" )
sigs.map2.filename <- paste( nap.output.dir , "luna_core_CANONICAL_CH_EDF-BASE.txt" , sep="" )

if ( file.exists( sigs.map1.filename ) ) {
df <- nap.read.table( sigs.map1.filename )
if ( dim(df)[1]>0) saveit( df , "sigs.base.map1" , file= paste( nap.output.dir , "nap.sig.base.map1.RData" , sep="" ) )
}

if ( file.exists( sigs.map2.filename ) ) {
df <- nap.read.table( sigs.map2.filename )
if ( dim(df)[1]>0) saveit( df , "sigs.base.map2" , file= paste( nap.output.dir , "nap.sig.base.map2.RData" , sep="" ) )
}


# ------------------------------------------------------------------------------------------------------
#
# SIGSTATS and STATS processing for Signals View (on the ORIGINAL EDFs) ;
# this is repeated in coda2.R for the Harmonized EDFs  --> 
#
# ------------------------------------------------------------------------------------------------------

# Hjorth
sigstats.filename <- paste( nap.output.dir , "luna_stats_SIGSTATS_E_CH.txt" , sep="" )

# Means
stats.filename <- paste( nap.output.dir , "luna_stats_STATS_E_CH.txt" , sep="" )   

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
 saveit( df , "sigstats" , file= paste( nap.output.dir , "nap.sigstats.RData" , sep="" ) )



# ------------------------------------------------------------------------------------------------------
#
# MTM spectrograms
#
# ------------------------------------------------------------------------------------------------------

if ( nap.domtm == 1 ) {
try( { 
 mtm.file <- paste( nap.output.dir , "luna_spec_MTM_F_CH_SEG.txt.gz" , sep="/" ) 
 if ( file.exists( mtm.file ) ) {
   mtm <- nap.read.table( mtm.file )
   mtm$MTM <- 10*log10( mtm$MTM )
   mtm <- mtm[ mtm$F >= 0.5 , ] 
   spectrograms <- list( desc = "MTM EEG spectrograms" )
   chs <- unique( mtm$CH )
   for (ch in chs)
   { 
   png( file = paste( nap.output.dir , paste( "mtm-" , ch , ".png" , sep="" ) , sep="/" ) , res=100 , width=1000  , height = 300 )
   par(mar=c(1.5,2.5,1.5,0.5))
   x <- mtm$SEG[ mtm$CH == ch ]
   y <- mtm$F[ mtm$CH == ch ]
   z <- mtm$MTM[ mtm$CH == ch ]
   z <- lwin( z , p = 0.1 )
   nx <- length(unique(x))
   ny <- length(unique(y))
   nz <- length(z)
   if (nz != nx * ny) stop("requires square data")
   d <- data.frame(x, y, z)
   d <- d[order(d$y, d$x), ]
   m <- matrix(d$z, byrow = T, nrow = ny, ncol = nx)
   hmcols <- lturbo(100)
   image(t(m[1:ny, ]), col = hmcols, xaxt = "n", yaxt = "n")
   axis(1,0.95,labels="Epochs",tick=F,line=-1,cex.axis=0.8)
   axis(2,0.98,labels=paste(round(max(y)),"Hz ") ,line=-1,tick=F,las=2,cex.axis=0.8)
   axis(2,0.02,labels="0.5 Hz " ,tick=F,las=2,line=-1,cex.axis=0.8)
   axis(3,0.05,labels=paste( "Channel:" , ch ) , tick=F , line=-1, cex.axis=1 )
   dev.off()
   spectrograms[[ ch ]] <- list( desc = ch , figure = paste( "mtm-" , ch , ".png" , sep="" )  )
 }  # next EEG channel

 if ( length( spectrograms ) > 0 ) 
   saveit( spectrograms , str = "spec" , file = paste( nap.output.dir , "mtm-spectrograms-fig.RData" , sep="" ) )
}
}) # end of MTM block
}