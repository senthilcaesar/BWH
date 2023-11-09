library(stringr)


# Set source path working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# APPLES
d <- read.csv("apples-fileid-demographics-smp-cut-20230215-site.csv", stringsAsFactors=F)
d <- d[ , c( "fileid" , "nsrr_age" , "nsrr_sex", "site") ]
d$site <- substr( d$site , 4, 6 )
names(d) <- c("ID","age","male", "site")
d$male <- ifelse( d$male == "male" , 1 , 0 )
write.table( d, file="files/apples.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# CCSHS
d <- read.csv("ccshs-trec-dataset-0.7.0.csv", stringsAsFactors=F)
d$ID <- paste("ccshs-trec", formatC(d$nsrrid), sep = "-")
d <- d[ , c( "ID" , "age" , "male") ]
write.table( d, file="files/ccshs.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# CFS
d <- read.csv("cfs-visit5-dataset-0.7.0.csv", stringsAsFactors=F)
d$ID <- paste( "cfs-visit5" , d$nsrrid , sep="-" )
d <- d[ , c( "ID" , "age" , "male") ]
write.table( d, file="files/cfs.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# CHAT baseline
d <- read.csv("chat-baseline-dataset-0.12.0.csv", stringsAsFactors=F)
d$ID <- paste( "chat-baseline" , d$nsrrid , sep="-" )
d <- d[ , c( "ID" , "ageyear_at_meas" , "male") ]
names(d) <- c("ID","age","male")
write.table( d, file="files/chat-baseline.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# CHAT followup
d <- read.csv("chat-followup-dataset-0.12.0.csv", stringsAsFactors=F)
d$ID <- paste( "chat-followup" , d$nsrrid , sep="-" )
d <- d[ , c( "ID" , "ageyear_at_meas" , "male") ]
names(d) <- c("ID","age","male")
write.table( d, file="files/chat-followup.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# CHAT baseline nonrandomized
d <- read.csv("chat-nonrandomized-dataset-0.12.0.csv", stringsAsFactors=F)
d$ID <- paste( "chat-baseline" , d$nsrrid , sep="-" )
d <- d[ , c( "ID" , "age_nr") ]
names(d) <- c("ID","age")
write.table( d, file="files/chat-nonrandomized.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# MESA
d <- read.csv("mesa-sleep-harmonized-dataset-0.6.0.csv", stringsAsFactors=F)
d$mesaid <- str_pad(d$mesaid, width=4, side="left", pad="0")
d$ID <- paste( "mesa-sleep" , d$mesaid , sep="-" )
d <- d[ , c( "ID" , "nsrr_age" , "nsrr_sex") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == "male" , 1 , 0 )
write.table( d, file="files/mesa.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# MROS Visit1
d <- read.csv("mros-visit1-dataset-0.6.0.csv", stringsAsFactors=F)
d$ID <- paste( "mros-visit1" , tolower(d$nsrrid) , sep="-" )
d <- d[ , c( "ID" , "vsage1" , "gender") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == 2 , 1 , 0 )
write.table( d, file="files/mros-visit1.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# MROS Visit2
d <- read.csv("mros-visit2-dataset-0.6.0.csv", stringsAsFactors=F)
d$ID <- paste( "mros-visit2" , tolower(d$nsrrid) , sep="-" )
d <- d[ , c( "ID" , "vs2age1" , "gender") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == 2 , 1 , 0 )
write.table( d, file="files/mros-visit2.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# MSP
d <- read.csv("msp-harmonized-dataset-0.1.1.csv", stringsAsFactors=F)
d <- d[ , c( "fileid" , "nsrr_age" , "nsrr_sex") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == "male" , 1 , 0 )
write.table( d, file="files/msp.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# NCHSDB
d <- read.csv("nchsdb-dataset-harmonized-0.3.0.csv", stringsAsFactors=F)
d <- d[ , c( "filename_id" , "nsrr_age" , "nsrr_sex") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == "male" , 1 , 0 )
write.table( d, file="files/nchsdb.dat" , sep="\t" , row.names=F, col.names=T , quote=F )


# SHHS 1
d <- read.csv("shhs1-dataset-0.20.0.csv", stringsAsFactors=F)
d$ID <- paste( "shhs1" , d$nsrrid , sep="-" )
d <- d[ , c( "ID" , "age_s1" , "gender") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == 1 , 1 , 0 )
write.table( d, file="files/shhs1.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

# SHHS 2
d <- read.csv("shhs2-dataset-0.20.0.csv", stringsAsFactors=F)
d$ID <- paste( "shhs2" , d$nsrrid , sep="-" )
d <- d[ , c( "ID" , "age_s2" , "gender") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == 1 , 1 , 0 )
write.table( d, file="files/shhs2.dat" , sep="\t" , row.names=F, col.names=T , quote=F )


# SOF
d <- read.csv("sof-visit-8-dataset-0.8.0.csv", stringsAsFactors=F)
d$sofid <- str_pad(d$sofid, width=5, side="left", pad="0")
d$ID <- paste( "sof-visit-8" , d$sofid , sep="-" )
d <- d[ , c( "ID" , "v8age" , "gender") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == 1 , 0 , 1 )
write.table( d, file="files/sof.dat" , sep="\t" , row.names=F, col.names=T , quote=F )


# STAGES
d <- read.csv("/Users/sq566/Desktop/saps/csv/stages-harmonized-dataset-0.3.0.csv", stringsAsFactors=F)
d <- d[ , c( "subject_code" , "nsrr_age" , "nsrr_sex") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == "male" , 1 , 0 )
write.table( d, file="/Users/sq566/Desktop/saps/csv/files/stages.dat" , sep="\t" , row.names=F, col.names=T , quote=F )


# WSC
d <- read.csv("/Users/sq566/Desktop/saps/csv/wsc-dataset-0.6.0.csv", stringsAsFactors=F)
d$ID <- paste( "wsc-visit", d$wsc_vst, "-", d$wsc_id , "-nsrr", sep="")
d <- d[ , c( "ID" , "age" , "sex") ]
names(d) <- c("ID","age","male")
d$male <- ifelse( d$male == "M" , 1 , 0 )
write.table( d, file="/Users/sq566/Desktop/saps/csv/files/wsc.dat" , sep="\t" , row.names=F, col.names=T , quote=F )

