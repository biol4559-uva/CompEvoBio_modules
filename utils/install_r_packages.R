### install libraries
  BiocManager::install("SRAdb",     lib="/project/biol4559-aob2x/biol4559-R-packages/")
  BiocManager::install("gdsfmt",    lib="/project/biol4559-aob2x/biol4559-R-packages/")
  BiocManager::install("SeqArray",  lib="/project/biol4559-aob2x/biol4559-R-packages/")
  install.packages("readxl",        lib="/project/biol4559-aob2x/biol4559-R-packages/")
  install.packages("ggplot2",        lib="/project/biol4559-aob2x/biol4559-R-packages/")

### scrath
  install.packages("ggplot2",        lib="/scratch/aob2x/biol4559/biol4559-R-packages/")
  install.packages("R.utils", lib="/scratch/aob2x/biol4559/biol4559-R-packages")

  cp -R /scratch/aob2x/biol4559/biol4559-R-packages /standard/vol186/bergland-lab/biol4559-aob2x/


### necessary libraries
  .libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(foreach)
  library(ggplot2)
  library(SeqArray)
  library(SRAdb)
  library(data.table)
  library(readxl)


  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("Rsamtools",lib="/project/biol4559-aob2x/biol4559-R-packages/")
