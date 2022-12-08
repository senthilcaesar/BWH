
# template code to make the _pheno-*.RData and _derived-*.RData files that Luna/shiny will load
# and display in either the Phenotypes or Sample metrics tabs

# this should be run in the parent folder of the main `nap/` output folder

# nb. codes p1, p2, d1, d2, etc do not show in the viewer

# --------------------------------------------------------------------------------
#
# Convenience function 
#
# --------------------------------------------------------------------------------

# for lload() -- to load and flatten tables (one row = one INDIV)
source( "http://zzz.bwh.harvard.edu/dist/luna/lload2.R" )

freader <- function( d , f , factors = NULL , vars = NULL ) {
 df <- lload( paste( d , f , sep="/" ) , factors = factors , variables = vars )$df
 if ( names(df)[1] != "ID" ) stop( "expecting column 1 to be ID" )
 df
}



# --------------------------------------------------------------------------------
#
# Derived metrics
#
# --------------------------------------------------------------------------------

# g <- list( desc = "group description" ,
#            d1 = list( desc = "description1" , data = data.frame() ) ,
#            d2 = list( desc = "description2" , data = data.frame() ) )

# name of folder in which derived metrics were compiled

f <- "derived1"

# 1) EDF header info
# X) TODO: annot count/duration for key annot types (e.g. arousal, apnea, etc) - putting 0 in as needed
# 2) Macro arch
# 3) SOAP/SUDS
# 4) EEG artifact
# 5) Micro-arch

#
# 1) EDF headers
#

g1 <- list(
      desc = "Original/harmonized datasets" , 

      d1 = list( desc = "EDF durations (secs)" , 
                 data = freader( f , "nap.nums" , vars = c("ORIG_SEC" , "HARM_SEC" , "DROP_SEC" ) ) ) , 

      d2 = list( desc = "Number of signals" , 
                 data = freader( f , "nap.nums" , vars = c("ORIG_NS" , "HARM_NS" , "DROP_NS" ) ) ) , 

      d3 = list( desc = "Number of annotations" , 
                 data = freader( f , "nap.nums" , vars = c("ORIG_NA" , "HARM_NA" , "DROP_NA" ) ) )  

)

save( g1 , file = "nap/_derived-HEADERS.RData" )


#
# Annotations
#


g2 <- list(
      desc = "Existing annotations" ,

      d1 = list( desc = "Apnea" ,
                 data = freader( f , "luna_core_ANNOTS_ANNOT_HARM_apnea.txt" ) ) , 

      d2 = list( desc = "Hypopnea" ,
                 data = freader( f , "luna_core_ANNOTS_ANNOT_HARM_hypopnea.txt" ) ) , 

      d3 = list( desc = "Desats" ,
                 data = freader( f , "luna_core_ANNOTS_ANNOT_HARM_desat.txt" ) ) , 

      d4 = list( desc = "Arousals" ,
                 data = freader( f , "luna_core_ANNOTS_ANNOT_HARM_arousal.txt" ) ) , 

      d5 = list( desc = "Movements" ,
                 data = freader( f , "luna_core_ANNOTS_ANNOT_HARM_movement.txt" ) ) , 

      d6 = list( desc = "Artifacts" ,
                 data = freader( f , "luna_core_ANNOTS_ANNOT_HARM_artifact.txt" ) )  

)

save( g2 , file = "nap/_derived-ANNOTS.RData" )


#
# 3) Macro architecture
#

g3 <- list(
      desc = "Macro-architecture" ,

      d0 = list( desc = "Coverage" , 
                 data = freader( f , "luna_core_SPANNING_HARM.txt" , vars = c( "VALID_N", "INVALID_N", "ANNOT_OVERLAP", "SPANNED_PCT",  "SPANNED_SEC" ) ) ) , 

      d1 = list( desc = "Absolute stage duration" ,
                 data = freader( f , "luna_macro_HYPNO.txt" , vars = c("MINS_N1","MINS_N2","MINS_N3","MINS_REM","TST") ) ) ,

      d2 = list( desc = "Relative stage duration" ,
                   data = freader( f , "luna_macro_HYPNO.txt" , vars = c("PCT_N1","PCT_N2","PCT_N3","PCT_REM") ) ) ,

      d3 = list( desc = "Efficiency/latency" ,
                 data = freader( f , "luna_macro_HYPNO.txt" , vars = c("WASO", "SLP_EFF","SLP_EFF2","SLP_MAIN_EFF","SLP_LAT","REM_LAT") ) ) ,

      d4 = list( desc = "Other metrics" ,
                   data = freader( f , "luna_macro_HYPNO.txt" , vars = c("TIB","TPST","TWT" ) ) ) ,

      d5 = list( desc = "NREM cycles" , 
                   data = freader( f , "luna_macro_HYPNO.txt" , vars = c( "NREMC", "NREMC_MINS") ) )

   )

save( g3 , file = "nap/_derived-MACRO.RData" )


#
# 4) SOAP/SUDS 
#

g4 <- list(
      desc = "SOAP/SUDS" , 

      d1 = list( desc = "SOAP metrics" ,
	         data = freader( f , "luna_suds_SOAP.txt" , vars = c( "K", "K3", "ACC", "ACC3" ) ) ) ,

      d2 = list( desc = "SOAP durations" ,
	           data = freader( f , "luna_suds_SOAP_SS.txt" , factors = "SS" , vars = c( "DUR_OBS" , "DUR_PRD" ) ) ) ,	

      d3 = list( desc = "SUDS metrics" ,
	         data = freader( f , "luna_suds_SUDS.txt" , vars = c( "K", "K3", "ACC", "ACC3" ) ) ) ,

      d4 = list( desc = "SUDS durations" ,
	           data = freader( f , "luna_suds_SUDS_SS.txt" , factors = "SS" , vars = c( "DUR_OBS" , "DUR_PRD" )  ) ) 

   )

save( g4 , file = "nap/_derived-SUDS.RData" )


#
# 5) Gross artifact (EEG)
#  - spectral peaks
#  - flat & clipped signals
#  - skew and DC offsets
#  - number of outlier epochs (w/ at least one EEG channel)

g5 <- list(
      desc = "Artifacts" , 

      d1 = list( desc = "Masked base EDF epochs" , 
                 data = freader( f , "bad.base.epochs.summary.txt" ) ) 

)

save( g5 , file = "nap/_derived-ARTIFACTS.RData" )


#
# 6) Micro-arch
#

g6 <- list(
      desc = "Micro-architecture" ,

      d1 = list( desc = "NREM spectral power" , 
                 data = freader( f , "luna_spec_PSD_B_CH_SS-N2.txt" , vars = c("PSD") , factors=c("B") ) ) ,

      d2 = list( desc = "N2 spindles" , 
                 data = freader( f , "luna_spso_SPINDLES_CH_F_SS-N2.txt" , vars = c("DENS") , factors="F" ) ) 

   )

save( g6 , file = "nap/_derived-MICRO.RData" )


# --------------------------------------------------------------------------------
#
# Phenotypes (templates; second usage extracts only those cols)
#
# --------------------------------------------------------------------------------

p1 <- freader( "." , "phe.txt" )
save( p1 , file = "nap/_pheno-p1.RData" )

#p2 <- freader( "folder" , "replace-me.txt" , vars = c("AGE","SEX","BMI") )
#save( p2 , file = "nap/_pheno-p2.RData" )

# or, if wanting a single file (doesn't matter from Luna/Shiny perspective)

#p <- merge( p1 , p2 , by="ID" , all = T )
#save( p , file = "nap/_pheno.RData" )

