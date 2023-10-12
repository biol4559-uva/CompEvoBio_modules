### libraries
  install.packages("curl")

  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(5)

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
  genofile <- seqOpen("/scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.gds")
  genofile

### load meta-data file
  samps <- fread("/scratch/aob2x/GDS/biol4559_sampleMetadata.csv")
  samps[is.na(set), set:="ExpEvo"]

### load in inversion markers
  inv.markers <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/main/Module_8/inversion_markers_v6.txt")
  setnames(inv.markers, c("chrom", "position"), c("chr", "pos"))

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency

### your turn: our goal is to retrieve the variant.ids for each inversion marker SNP.
### That way we can query our GDS object for allele frequencies.
### you will need to use the "merge" command from the data.table package
### how many of the inversion markers do not intersect with our set? Do you feel comfortable with that level of intersection?
  markers <- merge(snp.dt, inv.markers, by=c("chr", "pos"))

### get the allele frequencies
  dat <- getData(snps=markers, samples=samps[grepl("PRJNA194129", sampleId)])

### your turn: write a figure legend for this figure.
  ggplot(data=dat, aes(x=inversion, y=af_nEff)) + geom_point() + facet_wrap(~sampleId)

### your turn:
### calculate the frequency of the inversion for each of your samples by averaging the frequency of all SNPs associated with the inversion
  dat.ag <- dat[,list(freq=mean(af_nEff, na.rm=T)),
                 list(inversion, exp_rep, sampleId, bio_rep)]

### your turn:
### plot the estimated allele frequency of each inversion across treatments
  ggplot(data=dat.ag[!is.na(freq)], aes(x=exp_rep, y=freq, color=inversion)) +
    geom_jitter(width=.15) + facet_grid(~inversion)

### some simple models to test for differences in frequency between the treatment groups
  mod <- lm(freq~exp_rep, dat.ag[inversion=="In(2L)t"])
  summary(mod)
  anova(mod)

### extracting the p-value
  str(anova(mod))

  pvals <- data.table(inversion="In(2L)t", p=anova(mod)$Pr[1])
  ggplot(data=dat.ag[!is.na(freq)], aes(x=exp_rep, y=freq, color=inversion)) +
    geom_jitter(width=.15) + facet_grid(~inversion) +
    geom_text(data=pvals, aes(x=1.5, y=.5, label=paste("p=", round(p, 10), collapse="")))

  
