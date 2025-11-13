### libraries
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(poolfstat)

### load this function
  source("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/utils/misc/getData_function.R")

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/compBio_SNP_29Sept2025/dest.all.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.gds")
  genofile

### load meta-data file
  samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/full_sample_metadata.90Sept2025_ExpEvo.csv", fill=T)

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"),
                       refAllele=seqGetData(genofile, "$ref"),
                       altAllele=seqGetData(genofile, "$alt"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency

### get SNPs from best window
  ## 3L 7219979 7319979
  setkey(candidates, chr, pos)
  candidates <- merge(snp.dt, candidates)
  freqs <- foreach(i=1:dim(candidates)[1], .combine="rbind")%do%{
    getData(snps=merge(snp.dt, candidates[,c("chr", "pos")])[i], samples=samps)
  }

  ggplot(data=freqs[grepl("PRJEB5713", sampleId)], aes(x=exp_rep, y=af_nEff, group=variant.id)) + geom_line()

  freqs.ag <- freqs[grepl("PRJEB5713", sampleId)][,list(p=mean(af_nEff, na.rm=T)), list(exp_rep, variant.id)]
  ggplot(data=freqs.ag, aes(x=exp_rep, y=p, group=variant.id)) + geom_line()

  setkey(freqs, sampleId)
  freqs <- merge(freqs, samps[,c("Recommendation", "sampleId")])
  ggplot(data=freqs[Recommendation=="Pass"][continent=="Europe"], aes(x=as.numeric(as.character(lat)), y=af_nEff, group=variant.id)) + geom_line()
