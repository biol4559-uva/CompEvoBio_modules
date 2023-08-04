### install libraries
  BiocManager::install("SRAdb", lib="~/biol4559-R-packages/")
  BiocManager::install("gdsfmt", lib="~/biol4559-R-packages/")
  BiocManager::install("SeqArray", lib="~/biol4559-R-packages/")
  install.packages("readxl", lib="~/biol4559-R-packages/")

### necessary libraries
  .libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(foreach)
  library(ggplot2)
  library(SeqArray)
  library(SRAdb)
  library(data.table)
  library(readxl)
