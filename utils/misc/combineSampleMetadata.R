### libraries
  library(data.table)
  library(readxl)

### load data
  dest <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_3May2024.csv")
  g1 <- as.data.table(read_excel("/Users/alanbergland/BIOl 4020 - Sample Metadata.xlsx"))
  g2 <- as.data.table(read_excel("/Users/alanbergland/Metatdata_proj_paper_4020 (2).xlsx", range="A1:Y60")); g2 <- g2[,-"LibraryName", with=F]
  g3 <- as.data.table(read_excel("/Users/alanbergland/Martins_BIOL_4020_Group_Metadata.xlsx"))

### combine
  ee <- rbindlist(list(g1, g2, g3), fill=T)
  ee[,lat:=as.numeric(tstrsplit(lat, " ")[[1]])]
  ee[,long:=as.numeric(tstrsplit(long, " ")[[1]])]
  ee[continent=="North_America", long:=long*-1]
  ee <- ee[,-"subsample", with=F]

  
