library(GEOquery)
library(affy)

treatment <- function(gse) { 
  Meta(gse)[['source_name_ch1']]
}

write_csv <- function(gse_file_name) {
  # Assumes the .CEL and .soft files are present in a folder called new_data/
  gse <- getGEO(filename=paste0('new_data/', gse_file_name, '_family.soft'))
  pd <- data.frame(condition=as.factor(sapply(GSMList(gse), treatment)))
  celfiles <-  paste0('new-data/', rownames(pd), '.CEL')
  affydata <- read.affybatch(only_cel, phenoData = new("AnnotatedDataFrame", pd))
  eset <- rma(affydata)
  write.csv(eset@assayData[['exprs']], paste0(gse_file_name, '.csv'))
}

