
### libraries
  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### collect
  fl <- list.files("/scratch/aob2x/Fst_out/", ".Rdata", full.names=T)
  out <- foreach(fl.i=fl)%dopar%{
    load(fl)
    return(out)
  }
  out <- rbindlist(out, fill=T)

### save big object
  save(out, file="/scratch/aob2x/big_fstOutput.Rdata")

### your turn
### generate plot of the histogram of Fst values across combinations of comparisions (e.g., control vs treatment) for each chromosome

### let's generate a sliding window analysis calculating the average Fst values between comparisons
  ### build a dictionary of windows
    win.bp <- 4e4
    step.bp <- 4e3 ### note how the step size is the NOT same as the window size. These are overlapping windows....

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
    setkey(out, chr)

    windowedFst <- foreach(i=1:dim(wins)[1])%dopar%{
      # i <- 1000
      tmp <- out[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
      tmp
    }

### how can we test if you have a bunch of outliers in your window?
### we can calculate for each comparision of interst the normalize Fst. Basically, is the Fst value in the top 1%, top 5%, etc.
### If each SNP is independent, then if we take 1000 SNPs we expect 5% of them to be in the top 5%.
### because of linkage, if a region is subjected to selection and has a bunch of high Fst sites then we expect more than 5% to be in the top 5%, etc.

  out.norm <- out[,list(normFst=1 - rank(aveFst)/(1+length(aveFst)), chr, pos, aveFst), list(trt1, trt2)]
  out.norm[aveFst>.25]
  setkey(out.norm, chr)
  windowedFst <- foreach(i=1:dim(wins)[1])%dopar%{
    # i <- 10
    tmp <- out.norm[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
    tmp.ag <- tmp[,list(nSNPs=.N, pr05=sum(normFst<=.05), meanPos=mean(pos), window=i), list(trt1, trt2, chr)]
    tmp.ag[,p_value:=pbinom(pr05, nSNPs, .05, lower.tail=F)]
    tmp.ag
  }
