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

### load DEST metadata
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv")

### combine
  sampse <- rbind(samps, eem, fill=T)
  dim(samps)
  dim(sampse)
  
