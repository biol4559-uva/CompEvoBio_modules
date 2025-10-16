# ijob -A berglandlab_standard -c1 -p standard --mem=6G
### module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; R


### libraries
  library(data.table)
  library(readxl)
  library(foreach)
  library(SeqArray)

### load metadata
  xl.fn <- "/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/ExpEvo_meta.xlsx"
  #xl.fn <- "~/CompEvoBio_modules/data/ExpEvo_meta.xlsx"

  eem <- foreach(i=excel_sheets(xl.fn))%do%{
    read_excel(xl.fn, i)
  }
  eem <- rbindlist(eem, fill=T)
  eem
  eem[,locality:=tstrsplit(sampleId, "_")[[2]]]
  eem

### load DEST metadata
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

### semi-curated
  badSamp <- fread("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/full_sample_metadata.90Sept2025_ExpEvo.csv")
  eem <- merge(eem, badSamp, by="sampleId", all.x=T)
  eem[,Recommendation:="Pass"]
  eem

### rbidn
  samps2 <- rbind(samps, eem, fill=T)
  save(samps, file="/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/full_sample_metadata.90Sept2025_ExpEvo.csv", quote=F, row.names=F)









### combine
  sampse <- rbind(samps, eem, fill=T)
  dim(samps)
  dim(sampse)
  table(sapply(sampse$sampleId, function(x) x%in%s.dt$sampleId.gds))

### save(
  #write.csv(sampse, file="/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/full_sample_metadata.90Sept2025_ExpEvo.csv", quote=F, row.names=F)

### Get GDS file sample names
  genofile <- seqOpen("/standard/BerglandTeach/data/dest.all.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.gds")
  s.dt <- data.table(sampleId.gds=seqGetData(genofile, "sample.id"), gds=T)

### merge and simplify
  all <- merge(eem[,c("sampleId","nFlies")], s.dt[grepl("ExpEvo", sampleId.gds)], by.x="sampleId", by.y="sampleId.gds", all=T)
  all[grepl("0455", sampleId)]
  all[grepl("5429", sampleId)][gds==T][is.na(nFlies)]
  tmp <- all[grepl("5429", sampleId)][is.na(gds)][!is.na(nFlies)]
  tmp[,newName:=gsub("MM", "", sampleId)]

  sampse <- merge(sampse, tmp[,c("sampleId", "newName")], all.x=T)
  sampse[!is.na(newName), sampleId:=newName]

### check
  ### are all of the samples in the metadata in the ``
  table(sapply(sampse$sampleId, function(x) x%in%s.dt$sampleId.gds))
  all2 <- merge(sampse[,c("sampleId","nFlies")], s.dt, by.x="sampleId", by.y="sampleId.gds", all=T)
  all2[is.na(gds)]
  s.dt[grepl("_PRJNA657615", sampleId.gds)]
  all2[grepl("_PRJNA657615", sampleId)][is.na(gds)]


  table(sapply(sampse[!grepl("_PRJNA657615", sampleId) & !grepl("gen", sampleId)]$sampleId, function(x) x%in%s.dt$sampleId.gds))

  all2 <- merge(sampse[!grepl("_PRJNA657615", sampleId) & !grepl("gen", sampleId)][,c("sampleId","nFlies")], s.dt, by.x="sampleId", by.y="sampleId.gds", all=T)
  all2[is.na(gds)]

### all 3
  samps.ag <- all2[,list(.N), list(sampleId)]
  table(samps.ag$N)
  all3 <- all2[-as.numeric(sapply(samps.ag[N==2]$sampleId, function(x) which(all2$sampleId==x)[1]))]
  all3


  write.csv(all3, file="~/full_sample_metadata.90Sept2025_ExpEvo.csv", quote=F, row.names=F)


### fix it again
  sampse <- fread("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/full_sample_metadata.90Sept2025_ExpEvo.csv")
  samps <- merge(samps, sampse[,c("sampleId", "gds")], all=T)
  tmp <- samps[,.N,list(sampleId)]
  table(samps$Recommendation)
  samps[is.na(Recommendation), Recommendation:="Pass"]
  samps[is.na(locality),locality:=tstrsplit(sampleId, "_")[[2]]]
  samps
  table(tmp$N)
