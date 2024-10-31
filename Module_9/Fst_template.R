### libraries

  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(4)

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
  setkey(snp.dt, chr, pos)

  system.time(out <- foreach(i=1:dim(wins)[1], .errorhandling="remove")%dopar%{
    #i <- 100
    # message(i)
    ### get data for your sample for this window
    focalSNPs <-snp.dt[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
    tmp.data <- getData(snps=focalSNPs, samples=samps[grepl("PRJEB5713", sampleId)])
    tmp.data[,window:=i]

    ### Now we summarize. Examples of summary statistics include: missing data per sample; average coverage; average frequency of mutations; averge difference in allele frequency between treatments? What else?
    ### your turn
    tmp.data[,list(coverge=mean(dp, na.rm=T), missing=mean(is.na(dp))), list(sampleId, window)]
    tmp1 <- tmp.data[,list(delta_contsysl_virsys=mean(af_nEff[exp_rep=="contsys"], na.rm=T) - mean(af_nEff[exp_rep=="virsys"], na.rm=T)), list(variant.id, window)]
    tmp1[,list(mean_delta_contsysl_virsys=mean(delta_contsysl_virsys, na.rm=T)), list(window)]

    ### merge with polymorphic sites
    tmp.data <- merge(tmp.data, bps, all.y=T, by=c("chr","pos"))
    tmp.data[is.na(af_nEff), af_nEff:=0]

    ### Now we summarize. Examples of summary statistics include: missing data per sample; average coverage; average frequency of mutations; average difference in allele frequency between treatments? What else?
    ### Try something simple first. Then, try to calculate differences in allele frequency between two groups, averaged within a window.
    # tmp.data[,list(coverge=mean(dp, na.rm=T), missing=mean(is.na(dp)), het=mean(2*af*(1-af), na.rm=T)), list(sampleId, window)]

    tmp.data[,list(coverge=mean(dp, na.rm=T), missing=mean(is.na(dp)),
                   het=mean(2*af*(1-af), na.rm=T),
                   het_correct=mean(2*af_nEff*(1-af_nEff))), list(sampleId, window)]



    tmp1 <- tmp.data[,list(delta_contsysl_virsys=mean(af_nEff[exp_rep=="contsys"], na.rm=T) - mean(af_nEff[exp_rep=="virsys"], na.rm=T)),
                     list(variant.id, window)]


  })
  out <- rbindlist(out)
  setkey(out, window)
  out <- merge(out, wins)
  ggplot(data=out, aes(x=window, y=mean_delta_contsysl_virsys, color=chr)) + geom_line()


### another version
out[,mid:=start/2 + end/2]


inversion.bp <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/main/Module_8/InversionsMap_hglft_v6_inv_startStop.txt")
ggplot(data=out, aes(x=mid, y=mean_delta_contsysl_virsys, color=chr)) + geom_line() + facet_grid(~chr, scales="free_x") +
  geom_vline(data=inversion.bp, aes(xintercept=start )) +
  geom_vline(data=inversion.bp, aes(xintercept=stop ))
