### capture command line arguments
  args = commandArgs(trailingOnly=TRUE)
  slurmJob=as.numeric(args[1])

### libraries
  library(curl)
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(4)


### source a few functions
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_input_freqs.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_correct.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_nocorrect.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Likelihood%20functions%20for%20OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/FST%20functions.R")

### load this function
  source("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/utils/misc/getData_function.R")

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/dest.expevo.PoolSeq.PoolSNP.001.50.28Sept2024_ExpEvo.norep.gds")
  genofile

### load meta-data file
  samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/full_sample_metadata.28Sept2024_ExpEvo.csv")

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency

### build a dictionary of windows
  win.bp <-  5e5
  step.bp <- 5e5

  setkey(snp.dt, "chr")

  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),.combine="rbind", .errorhandling="remove")%dopar%{
    # chr.i <- "2L"
    tmp <- snp.dt[J(chr.i)]
    data.table(chr=chr.i,
               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }

  wins[,window:=1:dim(wins)[1]]
  dim(wins)

### iterate thorugh windows to calculate a summary of your data
  setkey(snp.dt, chr, pos)

  # slurmJob <- 50
  out <- foreach(i=slurmJob, .errorhandling="remove")%do%{
    #i <- 120
    # message(i)
    ### get data for your sample for this window
      focalSNPs <-snp.dt[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
      tmp.data <- getData(snps=focalSNPs, samples=samps[grepl("SRP002024", sampleId)])
      tmp.data[,window:=i]
      tmp.data.clean <- tmp.data[!is.na(af_nEff)]

      fst.out <- foreach(i=unique(tmp.data.clean$variant.id), .errorhandling="remove")%dopar%{
        #i <- unique(tmp.data$variant.id)[1]; i<-544048
        snp.data <- tmp.data.clean[variant.id==i & !is.na(af_nEff)]

        if(all(tmp.data[variant.id==i & !is.na(af_nEff)]$af_nEff==0) | all(tmp.data[variant.id==i & !is.na(af_nEff)]$af_nEff==1)) {
          return(NULL)
        }

        if(dim(snp.data)[1]==1) {
          return(NULL)
        }

        fst.cor <- outflank_mod_fst_correct(snp.data$af_nEff, snp.data$nEff, Ho=NULL)
        fst.nocor <- outflank_mod_fst_nocorrect(snp.data$af_nEff, snp.data$nEff, Ho=NULL)

        if(length(fst.cor)==1) return(NULL)

        # An FST data frame
        fstDF <- data.frame(variant.id=i, nPops=length(tmp.data[variant.id==i & !is.na(af_nEff)]$af_nEff),
                            cbind(as.data.frame(fst.nocor[c('FSTNoCorr', 'T1NoCorr', 'T2NoCorr')]), as.data.frame(fst.cor)))
        fstDF

      }
      fst.out <- rbindlist(fst.out)
      if(dim(fst.out)[1]==0) {
        fst.out <- data.table(variant.id=unique( tmp.data$variant.id))
      }

    ### format output
      merge(fst.out,
            tmp.data[,list(annotation=unique(annotation), gene=unique(gene), meanRD=mean(af_nEff, na.rm=T), meanAF=mean(af_nEff, na.rm=T)), list(variant.id, chr, pos)],
            by="variant.id", all.y=T)
  }
  out <- rbindlist(out, fill=T)

### save output
    outputdir="/scratch/aob2x/fstoutput/"
    if(!dir.exists(outputdir)) {
      dir.create(outputdir)
    }
    save(out, file=paste(outputdir, slurmJob, "window.Rdata", sep=""))
