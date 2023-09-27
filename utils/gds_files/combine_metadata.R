### libraries
  library(data.table)
  library(readxl)

### load meta-data file
  samps <- fread("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/utils/gds_files/dest_v2.samps_8Jun2023.csv")

### expevo samples
  expevo <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/utils/gds_files/expevo_samps.xlsx"))
  expevo[long<0 & grepl("North", continent), long:=long*-1]
  table(expevo$long<0)

### combine
  samps <- rbind(samps, expevo, fill=T)
  samps[sampling_strategy=="ExpEvo"]
  table(samps$set)

### export
  write.csv(samps, quote=F, row.names=F, file="/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/utils/gds_files/biol4559_sampleMetadata.csv")
