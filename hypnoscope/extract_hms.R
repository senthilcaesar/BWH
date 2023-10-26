
fileName <- '/Users/senthilp/Desktop/test.annot'
n <- 3
hms_regex <- "^([01]?[0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9]$"
pat <- paste0('^([^:]+(?::[^:]+){',n-1,'}).*')

conn <- file(fileName,open="r")
linn <-readLines(conn)
for (i in 1:length(linn)){
  line <- gsub('\t', " ", linn[i], fixed = TRUE)
  line2 <- strsplit(line, " ")[[1]]
  startTime <- line2[4]
  startTime <- sub(pat, '\\1', startTime)
  
  if (grepl(hms_regex, startTime)) {
    print(startTime)
    break
  }
}
close(conn)

