### libraries
  install.packages("curl")
  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### load this function
  getData <- function(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps) {
    # snps=snp.dt[pos==14617051 & chr=="2L"]; samples=samps[grepl("SRP002024", sampleId)]$sampleId

    ### filter to target
    seqSetFilter(genofile, variant.id=snps$id, sample.id=samples$sampleId)

    ### get annotations
    message("Annotations")
    tmp <- seqGetData(genofile, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snps$id, times=len1))

    # Extract data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                          list(variant.id=id)]

    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

    ### get frequencies
    message("Allele Freqs")

    ad <- seqGetData(genofile, "annotation/format/AD")
    dp <- seqGetData(genofile, "annotation/format/DP")

    if(class(dp)[1]!="SeqVarDataList") {
      dp <- list(data=dp)
    }


    af <- data.table(ad=expand.grid(ad$data)[,1],
                     dp=expand.grid(dp$data)[,1],
                     sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                     variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    ### tack them together
    message("merge")
    afi <- merge(af, snp.dt1.an, by="variant.id")
    afi <- merge(afi, snps, by.x="variant.id", by.y="id")

    afi[,af:=ad/dp]

    ### calculate effective read-depth
    afis <- merge(afi, samples[,c("sampleId", "set", "nFlies", "locality",
                                  "lat", "long", "continent", "country", "province", "city",
                                  "min_day", "max_day", "min_month", "max_month", "year", "jday",
                                  "bio_rep", "tech_rep", "exp_rep", "loc_rep", "subsample", "sampling_strategy",
                                  "SRA_Accession"), with=F], by="sampleId")

    afis[chr=="X|Y", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
    afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
    afis[,af_nEff:=round(af*nEff)/nEff]
    setnames(afis, "col", "annotation")
    ### return
    afis[,-c("n"), with=F]
  }

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/GDS/dest.expevo.PoolSNP.001.50.25Sept2023.norep.ann.gds")
  genofile

### load meta-data file
  samps <- fread("/scratch/aob2x/GDS/biol4559_sampleMetadata.csv")
  samps[is.na(set), set:="ExpEvo"]

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
