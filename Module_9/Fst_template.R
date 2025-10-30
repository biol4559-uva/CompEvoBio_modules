args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

## install packages
  # install.packages("poolfstat")

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

### build a dictionary of windows
  win.bp <- 1e5
  step.bp <- 1e5 ### note how the step size is the same as the window size. These are non-overlapping windows....

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

  system.time(out <- foreach(i=1:100, .errorhandling="remove")%dopar%{
    #i <- 100
    message(i)
    ### get data for your sample for this window
      focalSNPs <- snp.dt[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
      tmp.data <- getData(snps=focalSNPs, samples=samps[grepl("PRJEB5713", sampleId)])
      tmp.data[,window:=i]
      tmp.data[is.na(dp), dp:=0]
      tmp.data[is.na(ad), ad:=0]
      tmp.data.summary <- tmp.data[,list(nFlies=mean(nFlies)), list(sampleId)]
      tmp.data[,refDep:=dp-ad]

    ### format conversion for poolfstat
      tmp.data.ref.mat <- dcast(data=tmp.data[,c("sampleId", "variant.id", "refDep")], variant.id~sampleId, value.var="refDep")
      tmp.data.dp.mat <- dcast(data=tmp.data[,c("sampleId", "variant.id", "dp")], variant.id~sampleId, value.var="dp")


      setnames(focalSNPs, c("chr", "pos", "refAllele", "altAllele"), c("Chromosome", "Position", "RefAllele", "AltAllele"))

      tmp.pool <- new("pooldata",
                  npools=dim(tmp.data.summary)[1],
                  nsnp=dim(tmp.data.dp.mat)[1],
                  refallele.readcount=as.matrix(tmp.data.ref.mat[,-"variant.id"]),
                  readcoverage=as.matrix(tmp.data.dp.mat[,-"variant.id"]),
                  poolsizes= tmp.data.summary$nFlies * 2,
                  poolnames = tmp.data.summary$sampleId,
                  snp.info = focalSNPs[,c("Chromosome", "Position", "RefAllele", "AltAllele")])

    ### calculte pairwise Fst
      fst.out <- compute.pairwiseFST(x=tmp.pool, output.snp.values=T)

    ### format Fst output
      fst.snp.mat <- fst.out@PairwiseSnpFST
      rownames(fst.snp.mat) <- focalSNPs$id

      fst.snp.long <- as.data.table(reshape2::melt(fst.snp.mat))
      fst.snp.long[,samp1:=tstrsplit(Var2, ";")[[1]]]
      fst.snp.long[,samp2:=tstrsplit(Var2, ";")[[2]]]
      fst.snp.long <- merge(fst.snp.long, focalSNPs, by.x="Var1", by.y="id")
      setnames(fst.snp.long, "value", "Fst")

    ### format Htot output
      Q1.snp.mat <- fst.out@PairwiseSnpQ1
      rownames(Q1.snp.mat) <- focalSNPs$id
      Q1.snp.long <- as.data.table(reshape2::melt(Q1.snp.mat))

    ### format Htot output
      Q2.snp.mat <- fst.out@PairwiseSnpQ2
      rownames(Q2.snp.mat) <- focalSNPs$id
      Q2.snp.long <- as.data.table(reshape2::melt(Q2.snp.mat))

    ### combine
      qs <- merge(Q1.snp.long, Q2.snp.long, by=c("Var1", "Var2"))
      fst.snp.long <- merge(fst.snp.long, qs , by=c("Var1", "Var2"))

      fst.snp.long <- fst.snp.long[,-c("Var1", "Var2", "af")]
      fst.snp.long[,Htot:=1-value.y]
      fst.snp.long[,Hwith:=1-value.x]

    ### your turn
    ### summarize the Fst values.
    ### summarize by taking the average Fst within and between groups per SNP
    ### also remember to take the average Htot and Hwith values too.
      fst.snp.long[,trt1:=tstrsplit(samp1, "_")[[3]]]
      fst.snp.long[,trt2:=tstrsplit(samp2, "_")[[3]]]


    ### your turn
    ### make sure to return the data you want
      fst.snp.long[,list(aveFst=mean(Fst, na.rm=T)), list(trt1, trt2, chr=Chromosome, pos=Position)][!is.nan(aveFst)]

  })
  out <- rbindlist(out)
  system("mkdir /scratch/aob2x/Fst_out")
  save(out, file=paste("/scratch/aob2x/Fst_out/window", i, ".Rdata", sep=""))
  out
