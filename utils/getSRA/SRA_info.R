# ijob -A berglandlab_standard -c1 -p standard --mem=6G
### module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; R


### necessary libraries
  #.libPaths(c("~/biol4559-R-packages/", .libPaths()))
  library(foreach)
  library(ggplot2)
  library(readxl)
  library(data.table)

### load excel sheet
  sras <- na.omit(as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions_v2.xlsx", sheet="Papers")))

### approx how much space do we need?
  info <- foreach(acc=sras$accession)%do%{
    # acc <- "PRJNA657615"
    as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions_v2.xlsx", sheet=acc))

  }
  info <- rbindlist(info[-1], fill=T)

### just return vector of individual accession IDs
  runs <- foreach(acc=sras$accession, .combine="rbind")%do%{
    # acc <- sras$accession[1]
    tmp <- as.data.table(read_xlsx("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/SRA_accessions_v2.xlsx", sheet=acc))
    data.table(project=acc, run=tmp$Run)
  }
  write.table(runs, file="/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/data/runs_v2.csv", quote=F, row.names=F, col.names=F)

### plots
  ggplot(data=info, aes(x=Model, y=bases)) + geom_point()
  sum(info$size_MB)


### double check that all the files downloaded
  dl <- foreach(prj=runs$V1, .combine="rbind", .errorhandling="remove")%do%{
    #prj <- runs$V1[2]
    fl <- list.files(paste("/scratch/aob2x/compBio/fastq/", prj, sep=""))
    dl <- data.table(proj=prj, run=fl)
    dl[,run:=gsub(".fastq.gz", "", run)]
    dl[,run:=gsub(".fastq", "", run)]

    dl[,run:=tstrsplit(run, "_")[[1]]]
    dl
  }

### merge
  dl.ag <- dl[,list(N=length(unique(run)), run=unique(run)),proj]

  rdl <- merge(runs, dl.ag, by.x="V1", by.y="proj", all.x=T)
  rdl.ag <- rdl[,list(missing=mean(is.na(N))), list(V1)]
  table(is.na(rdl$N))

  missing <- merge(runs[V1=="PRJNA657615"], dl.ag[proj=="PRJNA657615"], by.x="V2", by.y="run", all.x=T)[is.na(proj)]
  write.table(missing, file="~/CompEvoBio_modules/data/runs_missing.csv", quote=F, row.names=F, col.names=F, sep=",")

#### sras
  runs <- fread("~/CompEvoBio_modules/data/runs.csv", header=F)
  sras <- gsub(".sra", "", list.files("/scratch/aob2x/compBio/sra/"))
  runs[!runs$V2%in%sras]

### 
