### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(doMC)
  registerDoMC(5)
  #install.packages("dplyr")
  library(dplyr)

### load big object
  load(file="/scratch/aob2x/big_fstOutput.Rdata")

### are the same regions outliers for different treatments?
  fst.out.w <- dcast(fst.out, chr+pos~trt1+trt2, value.var="aveFst", fun.aggregate=mean, na.rm=T)
  fst.out.w

  tab <- table(fst.out.w$ContSys_Control>.3, fst.out.w$Control_VirSys>.3)
  fisher.test(tab)
  tab

  enrichment.test <- foreach(thr=seq(from=0, to=.4, by=.01), .errorhandling="remove", .combine="rbind")%dopar%{
    message(thr)
    tab <- table(fst.out.w$ContSys_VirSys>thr, fst.out.w$Control_VirSys>thr)
    fet <- fisher.test(tab)
    data.table(thr=thr, or=fet$estimate, p=-log10(fet$p.value))
  }
  ggplot(data=enrichment.test, aes(x=thr, y=or)) + geom_line()

  thr <- .3
  fst.out.w[ContSys_VirSys>thr & Control_VirSys>thr]

### clinal enrichment
  clinal <- fread("/standard/BerglandTeach/misc/clinal_table.csv")
  setkey(clinal, chr, pos)
  setkey(fst.out.w, chr, pos)
  cf <- merge(clinal, fst.out.w)

  thr <- .3
  tab <- table(cf$clinal.p<.005, cf$ContSys_Control>thr & cf$Control_VirSys>thr)
  fisher.test(tab)


### seasonal enrichment
  seasonal <- fread("/standard/BerglandTeach/misc/seasonal_table.csv")
  seasonal[,set:="seasonal"]
  setkey(seasonal, chr, pos)
  setkey(fst.out.w, chr, pos)
  sf <- merge(seasonal, fst.out.w)

  thr <- .2
  tab <- table(sf$seas.p<.005, cf$ContSys_Control>thr & cf$Control_VirSys>thr)
  fisher.test(tab)

### windowed
  out.norm <- fst.out[,list(normFst=1 - rank(aveFst)/(1+length(aveFst)), chr, pos, aveFst), list(trt1, trt2)]
  out.norm[aveFst>.25]
  setkey(out.norm, chr)
  windowedFst.pval <- foreach(i=1:dim(wins)[1])%dopar%{
    # i <- 10
    message(i)
    tmp <- out.norm[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
    tmp.ag <- tmp[,list(nSNPs=.N, pr05=sum(normFst<=.05), aveFst=mean(aveFst), meanPos=mean(pos), window=i), list(trt1, trt2, chr)]
    tmp.ag
  }
  windowedFst.pval <- rbindlist(windowedFst.pval)
  windowedFst.pval[,pvalue:=pbinom(pr05, nSNPs, .05, lower.tail=F)]
  ggplot(data=windowedFst.pval, aes(x=window, y=aveFst)) + geom_line() + facet_grid(trt1~trt2) ### your turn - work with this plot to make it look nicer
  ggplot(data=windowedFst.pval, aes(x=meanPos, y=-log10(pvalue))) + geom_line() + facet_grid(trt1+trt2~chr) ### your turn - work with this plot to make it look nicer


  
