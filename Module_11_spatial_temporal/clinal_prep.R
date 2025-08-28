# ijob -A berglandlab -c5 -p standard --mem=20G
#   module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; echo "R_LIBS_USER=~/R/goolf/4.3" > ~/.Renviron; R

### libraries
  library(data.table)


### load core20 data
  core20.orig <- fread("/project/berglandlab/alan/drosRTEC/drosRTEC/mel_clinal_uniquepops.glm")
  #core20.swap <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")

  core20.orig[,set:="cline"]
  #core20.swap[,set:="swap"]
  #core20 <- rbind(core20.orig, core20.swap)
  core20 <- core20.orig

  setnames(core20, "chrom", "chr")

### R5 -> R6 DGRP conversion table
  liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
  liftover <- fread(liftover.fn)
  liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### do liftover
  setnames(core20, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
  setkey(core20, dm3_chr, dm3_pos)
  setkey(liftover, dm3_chr, dm3_pos)

  core20 <- merge(core20, liftover)

  setnames(core20, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

### save
  write.csv(core20, file="/standard/BerglandTeach/clinal_table.csv", quote=F, row.names=F)



### seasonal prep

### load core20 data
  core20.orig <- fread("/project/berglandlab/alan/drosRTEC/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.glm")
  #core20.swap <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")

  core20.orig[,set:="cline"]
  #core20.swap[,set:="swap"]
  #core20 <- rbind(core20.orig, core20.swap)
  core20 <- core20.orig

  setnames(core20, "chrom", "chr")

### R5 -> R6 DGRP conversion table
  liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
  liftover <- fread(liftover.fn)
  liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### do liftover
  setnames(core20, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
  setkey(core20, dm3_chr, dm3_pos)
  setkey(liftover, dm3_chr, dm3_pos)

  core20 <- merge(core20, liftover)

  setnames(core20, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

### save
  write.csv(core20, file="/standard/BerglandTeach/seasonal_table.csv", quote=F, row.names=F)
