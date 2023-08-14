### necessary libraries
  .libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(foreach)
  library(ggplot2)
  library(readxl)
  library(data.table)

### load excel sheet
  sras <- na.omit(as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions.xlsx")))

### approx how much space do we need?
  info <- foreach(acc=sras$accession)%do%{
    # acc <- "PRJNA657615"
    as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions.xlsx", sheet=acc))

  }
  info <- rbindlist(info[-1], fill=T)

### just return vector of individual accession IDs
  runs <- foreach(acc=sras$accession, .combine="rbind")%do%{
    # acc <- sras$accession[1]
    tmp <- as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions.xlsx", sheet=acc))
    data.table(project=acc, run=tmp$Run)
  }
  write.table(runs, file="/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/runs.csv", quote=F, row.names=F, col.names=F)

### plots
  ggplot(data=info, aes(x=Model, y=bases)) + geom_point()
  sum(info$size_MB)


### double check that all the files downloaded
  system("ls -lh > /scratch/aob2x//downloaded_sra")
  dl <- fread("~/downloaded_sra")
  dl <- dl[grepl("gz", V9)]
  dl[,run:=gsub(".fastq.gz", "", V9)]
  dl[,run:=gsub(".fastq", "", run)]
  dl[,run:=tstrsplit(run, "_")[[1]]]

  dl.ag <- dl[,.N,run]


  dl <- list

### merge
  rdl <- merge(runs, dl.ag, by="run", all.x=T)
  rdl.ag <- rdl[,list(missing=mean(is.na(N))), list(project)]
  table(is.na(rdl$N))


#### sras
  runs <- fread("~/CompEvoBio_modules/data/runs.csv", header=F)
  sras <- gsub(".sra", "", list.files("/scratch/aob2x/compBio/sra/"))
  runs[!runs$V2%in%sras]
