### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(doMC)
  registerDoMC(5)
  #install.packages("dplyr")
  library(dplyr)

### collect data
fl <- list.files("/scratch/aob2x/Fst_out", "window", full.names=T)
fst.out <- foreach(i=fl)%do%{
  load(i)
  return(out)
}
fst.out <- rbindlist(fst.out)

### save big object
  save(fst.out, file="/scratch/aob2x/big_fstOutput.Rdata")
  load(file="/scratch/aob2x/big_fstOutput.Rdata")

### your turn
### generate plot of the histogram of Fst values across combinations of comparisions (e.g., control vs treatment) for each chromosome

### your turn
### generate a "Manhattan" plot of the raw pvalues, but only for the extreme values
  ggplot(data=fst.out[aveFst>.15], aes(x=pos, y=aveFst)) + geom_point() + facet_grid(trt1+trt2~chr)

### let's generate a sliding window analysis calculating the average Fst values between comparisons
### build a dictionary of windows
  win.bp <- 1e5
  step.bp <- 5e4 ### note how the step size is the NOT same as the window size. These are overlapping windows....

  setkey(fst.out, "chr")

  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),.combine="rbind", .errorhandling="remove")%dopar%{
    # chr.i <- "2L"
    tmp <- fst.out[J(chr.i)]
    data.table(chr=chr.i,
               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }

  wins[,window:=1:dim(wins)[1]]
  dim(wins)

### Your turn
### iterate thorugh windows to calculate a summary of your data. What attributes of Fst do you want to record?
  setkey(fst.out, chr)

  windowedFst <- foreach(i=1:dim(wins)[1])%dopar%{
    # i <- 1000
    tmp <- fst.out[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
    tmp
  }

### how can we test if you have a bunch of outliers in your window?
### we can calculate for each comparision of interst the normalize Fst. Basically, is the Fst value in the top 1%, top 5%, etc.
### If each SNP is independent, then if we take 1000 SNPs we expect 5% of them to be in the top 5%.
### because of linkage, if a region is subjected to selection and has a bunch of high Fst sites then we expect more than 5% to be in the top 5%, etc.
### play around with the threshold (can you figure out how to set multiple thresholds? Does the signal change?)

  out.norm <- fst.out[,list(q_normFst=1 - rank(aveFst)/(1+length(aveFst)), chr, pos, aveFst), list(trt1, trt2)]
  out.norm[aveFst>.25]
  setkey(out.norm, chr)
  windowedFst.pval <- foreach(i=1:dim(wins)[1])%dopar%{
    # i <- 10
    message(i)
    tmp <- out.norm[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
    foreach(j=c(.001, .005, .01, .05), .combine="rbind")%do%{
      tmp.ag <- tmp[,list(nSNPs=.N, TR=sum(normFst<=j), thr=j, aveFst=mean(aveFst), meanPos=mean(pos), window=i, medianFst=median(aveFst)), list(trt1, trt2, chr)]
      tmp.ag
    }
  }
  windowedFst.pval <- rbindlist(windowedFst.pval)
  windowedFst.pval[,pvalue:=pbinom(TR, nSNPs, thr, lower.tail=F)]

  ggplot(data=windowedFst.pval, aes(x=meanPos, y=-log10(pvalue), group=thr, color=as.factor(thr))) + geom_line() + facet_grid(trt1+trt2~chr) ### your turn - work with this plot to make it look nicer
  ggplot(data=windowedFst.pval, aes(x=meanPos, y=aveFst)) + geom_line() + facet_grid(trt1+trt2~chr) ### your turn - work with this plot to make it look nicer

  ggplot(data=windowedFst.pval[chr=="3L"], aes(x=meanPos, y=-log10(pvalue), group=thr, color=as.factor(thr))) + geom_line() + facet_grid(trt1+trt2~chr) ### your turn - work with this plot to make it look nicer


  ggplot(data=windowedFst.pval[chr=="3L"][thr==.001], aes(x=-log10(pvalue), y=aveFst)) + geom_point()
  windowedFst.pval[chr=="3L"][thr==.001][which.max(aveFst)]
  wins[1114]

  setkey(fst.out, chr, pos)
  candidates <- out.norm[chr=="3L"][pos>=7219979][pos<=7319979][trt1=="Control"][trt2=="VirSys"][normFst<.0001]


  setkey(candidates, chr, pos)
  fst.out[J(candidates)][trt1=="Control"][trt2=="VirSys"]



### are certain classes of SNPs enriched for high-fst sites?
  fst.out[,annotation_simple:=case_when(
    annotation=="5_prime_UTR_premature_start_codon_gain_variant" ~ "start_stop",
    annotation=="downstream_gene_variant"~"intergenic_region",
    annotation=="initiator_codon_variant" ~"splice",
    annotation=="initiator_codon_variant&splice_region_variant" ~ "splice",
    annotation=="missense_variant&splice_region_variant"~"missense_variant",
    annotation=="splice_acceptor_variant&intron_variant"~"splice",
    annotation=="splice_donor_variant&intron_variant"~"splice",
    annotation=="splice_region_variant"~"splice",
    annotation=="splice_region_variant&intron_variant"~"splice",
    annotation=="splice_region_variant&non_coding_transcript_exon_variant"~"splice",
    annotation=="splice_region_variant&stop_retained_variant"~"splice",
    annotation=="splice_region_variant&synonymous_variant"~"splice",
    annotation=="start_lost&splice_region_variant"~"start_lost",
    annotation=="stop_gained&splice_region_variant"~"stop_gained",
    annotation=="stop_lost&splice_region_variant"~"stop_gained",
    annotation=="upstream_gene_variant"~"intergenic_region",
    annotation%in%c("stop_gained", "stop_lost", "stop_retained_variant", "start", "start_lost")~"start_stop",
    annotation=="start_lost"~"start_stop",
    annotation=="stop_gained"~"start_stop",
    annotation=="intergenic"~"intergenic_region",
    .default=NA)]
  fst.out[is.na(annotation_simple), annotation_simple:=annotation]
  fst.out[annotation=="start_lost", annotation_simple:="start_stop"]
  fst.out[annotation=="stop_gained", annotation_simple:="start_stop"]
  fst.out[annotation=="initiator_codon_variant&non_canonical_start_codon", annotation_simple:="start_stop"]
  fst.out

  ggplot(data=fst.out[annotation_simple%like%"missense" | annotation_simple%like%"synon"] , aes(x=annotation_simple, y=aveFst)) + geom_boxplot() + facet_grid(trt1~trt2)
  ggplot(data=fst.out[annotation_simple%like%"missense" | annotation_simple%like%"synon"] , aes(x=annotation_simple, y=Htot)) + geom_boxplot() + facet_grid(trt1~trt2)

  pnps.fst <-  fst.out[,list(nOutlier=sum(aveFst>.1), .N), list(annotation_simple, trt1, trt2)]
  pnps.fst[,prop:=nOutlier/N]
  pnps.fst[annotation_simple%like%"missense" | annotation_simple%like%"synon"]
  pnps.fst.wide <- dcast(data=pnps.fst[annotation_simple%like%"missense" | annotation_simple%like%"synon"], trt1+trt2~annotation_simple, value.var="nOutlier")
  pnps.fst.wide[,odds:=missense_variant/synonymous_variant]
  pnps.fst.wide

### Are there any genes that stand out? What does Flybase say about these genes? Do you use the window approach or the SNP-based outlier approach?


### can you identify candidate SNPs and look at their specific allele frequency change across replicates
  table(out.norm[normFst<.000001]$chr)



### get max Fst
  fst.out[,maf:= (1+sqrt(1-2*Htot))/2]
  fst.out[maf>.5, maf:=1-maf]
  fst.out[,fstMax:=maf/(1-maf)]
  fst.out[,fstBar:=aveFst/fstMax]

  fst.out[which.max(fstBar)]

  library(SeqArray)

  focalSNPs <- snp.dt[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
  tmp.data <- getData(snps=focalSNPs, samples=samps[grepl("PRJEB5713", sampleId)])
