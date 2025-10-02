### libraries
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(foreach)

### load metadata
  xl.fn <- "/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/ExpEvo_meta.xlsx"
  eem <- foreach(i=excel_sheets(xl.fn))%do%{
    read_excel(xl.fn, i)
  }
  eem <- rbindlist(eem, fill=T)
  eem
  eem[,locality:=tstrsplit(sampleId, "_")[[2]]]
  eem

### load DEST metadata
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

### combine
  sampse <- rbind(samps, eem, fill=T)
  dim(samps)
  dim(sampse)

### save(
  write.csv(sampse, file="/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/full_sample_metadata.90Sept2025_ExpEvo.csv", quote=F, row.names=F)
  
