### libraries
  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(10)

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
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency

### build a dictionary of windows
  win.bp <- 1e5
  step.bp <- 5e4

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
### note, that if you run this as a loop it will take about 10 minutes to run.
### However, you can easily run this script through a batch job in Rivanna to parallelize it
  setkey(snp.dt, chr, pos)

  ## windows.to.use <- 1:dim(wins)[1]     ### full set
   windows.to.use <- 1:10                 ### test set


  system.time(out <- foreach(i=wins.to.use, .errorhandling="remove")%dopar%{
    #i <- 100
    message(i)
    ### get data for your sample for this window
      focalSNPs <-snp.dt[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
      tmp.data <- getData(snps=focalSNPs, samples=samps[grepl("PRJEB5713", sampleId)])
      tmp.data[,window:=i]

    ### Now we summarize.
    ### Some summary statistics include: missing data per sample; average coverage; average frequency of mutations; averge difference in allele frequency between treatments? What else?
    ### the example below shows coverage and missing data rate.
      tmp.data.ag <- tmp.data[,list(coverage=mean(dp, na.rm=T), missing=mean(is.na(dp))), list(sampleId, window)]

    ### your turn You'll have to figure out how to return the object you want from each iteration of your foreach loop.


    ### return
      return(tmp.data.ag)
  })
  out <- rbindlist(out)
  setkey(out, window)
  out <- merge(out, wins)
