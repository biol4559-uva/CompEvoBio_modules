### necessary libraries
  .libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(SRAdb)
  library(data.table)
  library(readxl)

### load in excel table
  sras <- read_xlsx("~/CompEvoBio_modules/data/SRA_accessions.xlsx")
  sras <- as.data.table(sras)

###  
  getSRA(search_terms=sras$accession[1])
  