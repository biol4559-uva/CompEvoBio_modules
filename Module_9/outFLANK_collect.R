# ijob -A berglandlab -c5 -p largemem --mem=250G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.0.3; R

### libraries
  .libPaths(c("/scratch/aob2x/biol4559-R-packages-newer")); .libPaths()
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
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Likelihood%20functions%20for%20OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/FST%20functions.R")


### collect
  fl <- list.files("/scratch/aob2x/fstoutput/", "window.Rdata", full.names=T)
  out <- foreach(fl.i=fl)%dopar%{
    load(fl)
    return(out)
  }
  out <- rbindlist(out, fill=T)

### save big object
  save(out, file="/scratch/aob2x/fstoutput/fstOutput.Radta")

### for RITA
  load("/scratch/aob2x/fstOutput.Radta")
  out[,nPops:=6]

### outlier detection
  of <- OutFLANK(out[!is.na(FST)][FST>0 & FST<1], LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=6, qthreshold=0.05)
  save(of, file="/standard/biol4559-aob2x/RitaW_fstoutput/outFLANK_output.Rdata")
### how many outliers?
  str(of)

### make output plot
  ggplot(data=out, aes(FST)) + geom_histogram()

### your turn (you can choose):
  1) make a Manhattan plot of fst or PBS
  2) make a sliding window average Fst plot
  3) does Fst vary between chromosomes? What about inside inverted regions?



  outliers <- out[FSTNoCorr>...][annotation==...]
  setkey(outliers, id)
  outliers_nodup <- outliers[!duplicated(id)]

  dim(outliers)
  dim(outliers_nodup)

  setkey(out, id)
  randomIds <- sample(out$id, 45)
  out[J(randomIds)]
  
