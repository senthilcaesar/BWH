
# -------------------------------------------------------------------------------
# 
# Table/figure formats
#
# -------------------------------------------------------------------------------

# The Shiny app looks for files *-tab.RData  (tables)
#                           and *-fig.RData  (images)

# format expected by shiny app (w/ unique 'g's ):
#
# g <- list( desc = "group description" ,
#            d1 = list( desc = "description1" , data = data.frame() ) ,
#            d2 = list( desc = "description2" , data = data.frame() ) )


# g <- list( desc = "group description" ,
#            d1 = list( desc = "image description1" , figure = "fig1.png" ) , 
#            d2 = list( desc = "image description2" , figure = "fig2.png" ) )



# -------------------------------------------------------------------------------
# 
# Helper functions
#
# -------------------------------------------------------------------------------

# generic table reader

nap.read.table <- function( filename ) {
 read.table( filename , header = T , stringsAsFactors = F , comment.char = "" , sep="\t" )
}


# put this in h:m:s format

# save a dataframe with a given variable name
saveit <- function(dat, str, file) {
  x <- list(dat)
  names(x) <- str
  save(list=names(x), file=file, envir=list2env(x))
}



fextract <- function( outdir , cmd , desc , tables , transpose = rep( F , length( tables ) ) ) {

# if transposes not specified directly (i.e. differently) for each table:
# nb. otherwise, transpose assumes to match by **alphabetical** order of present tables

if ( length( transpose ) == 1 ) 
 transpose <- rep( transpose , length( tables ) )

# get all file names ending .txt or .txt.gz in the output folder
nap.files <- c( list.files( outdir , full.names = T , pattern = glob2rx( paste( cmd , "*.txt" , sep="" ) ) ) , 
                list.files( outdir , full.names = T , pattern = glob2rx( paste( cmd , "*.txt.gz" , sep="" ) ) ) ) 

# same, but ignore full paths (for display)
nap.files.short <- c( list.files( outdir , full.names = F , pattern = glob2rx( paste( cmd , "*.txt" , sep="" ) ) ) , 
                      list.files( outdir , full.names = F , pattern = glob2rx( paste( cmd , "*.txt.gz" , sep="" ) ) ) )
nap.files.short <- gsub( "*.txt" , "" , gsub( ".txt.gz" , "" , nap.files.short ) ) 

# only extract requested
inc <- nap.files.short %in% paste( cmd  , tables , sep="") 
nap.files <- nap.files[inc]
nap.files.short <- nap.files.short[inc]

# ensure 1+ files are given prior to output
if ( length( nap.files ) > 0 ) { 
 data  <- list( desc = desc ) 
 for (x in nap.files.short ) data[[ x ]] <- list( desc = x , data = NA )
 #lapply( nap.files.short , function( x ) { data[[ x ]] <<- list( desc = x , data = NA ) } )
 
 for (g in 1:length(nap.files)) {
  dt <- nap.read.table( nap.files[g] )
  if ( transpose[g] ) 
   { 
      dt <- as.data.frame( cbind( names(dt) , t( dt ) ) )  
      names(dt) <- c("ID" , paste("R",1:(dim(dt)[2]-1),sep="" ) ) 
   } 

# can remove ID (either col or row depending on transpose)
  if ( ! transpose[g] ) dt$ID <- NULL
  else dt <- dt[ dt$ID != "ID" , ]  

  data[[ nap.files.short[g] ]]$data <- dt

} 
 saveit( data , str = cmd , file = paste( outdir , cmd , "-tab.RData" , sep="" ) ) 
}
# all done
}

