### necessary libraries
  .libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(foreach)
  library(ggplot2)
  library(readxl)
  library(data.table)

### load excel sheet
  sras <- na.omit(as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions.xlsx")))

### how much space do we need?
  info <- foreach(acc=sras$accession)%do%{
    # acc <- "PRJNA657615"
    as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions.xlsx", sheet=acc))

  }
  info <- rbindlist(info[-1], fill=T)

### plots
  ggplot(data=info, aes(x=Model, y=bases)) + geom_point()
