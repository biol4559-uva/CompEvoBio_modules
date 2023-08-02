### install libraries
  BiocManager::install("SRAdb")
  BiocManager::install("gdsfmt")
  BiocManager::install("SeqArray")

### necessary libraries
  .libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(foreach)
  library(ggplot2)
  library(SeqArray)
  library(rlang)
  library(SRAdb)
  library(data.table)
  
