
### libraries
  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### source a few functions
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_input_freqs.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_correct.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_nocorrect.R")
  #source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Likelihood%20functions%20for%20OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/FST%20functions.R")
  source("/Users/alanbergland/Documents/GitHub/CompEvoBio_modules/Module_9/OutFLANK.R")

### collect
  fl <- list.files("/scratch/aob2x/outflank_output", "window.Rdata", full.names=T)
  out <- foreach(fl.i=fl)%dopar%{
    load(fl.i)
    return(out)
  }
  out <- rbindlist(out, fill=T)
  setnames(out, "meanAF", "meanAlleleFreq")
  setnames(out, "variant.id", "LocusName")

### save big object
  save(out, file="/scratch/aob2x/outFLANK.Rdata")
  load("~/outFLANK.Rdata")
  hist(out$FST)

### outlier detection
  tmp <- out[!is.na(FST)][FST>0 & FST<1][!is.na(FSTNoCorr)]

  tmp <- tmp[,c("LocusName", "He", "FST", "T1", "T2", "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")]
  of <- OutFLANK(as.data.frame(tmp), LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.01, NumberOfSamples=2, qthreshold=0.05)
  of <- as.data.table(of)

### summarize
  table(of$results.OutlierFlag)
  of[which.min(results.pvaluesRightTail)]

  ggplot(data=of, aes(x=results.FST, y=log10(results.pvaluesRightTail))) + geom_point()

### some basic plots 
