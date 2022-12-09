library(luna)
library(data.table)

sub_file <-  '/Users/sq566/cases.txt'
sub_vector <- readLines(sub_file)

for(i in 1:length(sub_vector)){
  subname <- sub_vector[i]
  filename <-  paste0('/Users/sq566/nap/', subname, '/luna_suds_SOAP-tab.RData')
  output <-  get(load(filename))
  soap_hyp <- output$luna_suds_SOAP_E$data$PRED
  fwrite(list(soap_hyp), file = paste0('/Users/sq566/analysis/soap_eannots/', subname ,'.eannot'))
}
