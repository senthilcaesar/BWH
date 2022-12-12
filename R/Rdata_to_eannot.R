library(luna)
library(data.table)

sub_file <-  '/Users/sq566/cases.txt'
sub_vector <- readLines(sub_file)

# POPS prediction to eannot
for(i in 1:length(sub_vector)){
  subname <- sub_vector[i]
  RdataFile <-  paste0('/Users/sq566/nap/', subname, '/luna_suds_POPS-tab.RData')
  output <-  get(load(RdataFile))
  pops_hyp <- output$luna_suds_POPS_E$data$PRED
  fwrite(list(pops_hyp), file=paste0('/Users/sq566/analysis/pops_eannots/', subname ,'.eannot'))
}

# SOAP prediction to eannot
for(i in 1:length(sub_vector)){
  subname <- sub_vector[i]
  RdataFile <-  paste0('/Users/sq566/nap/', subname, '/luna_suds_SOAP-tab.RData')
  output <-  get(load(RdataFile))
  soap_hyp <- output$luna_suds_SOAP_E$data$PRED
  fwrite(list(soap_hyp), file=paste0('/Users/sq566/analysis/soap_eannots/', subname ,'.eannot'))
}

