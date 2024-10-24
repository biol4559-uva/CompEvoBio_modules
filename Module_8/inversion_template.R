### libraries

  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### load this function
  source("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/utils/misc/getData_function.R")

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/dest.expevo.PoolSeq.PoolSNP.001.50.28Sept2024_ExpEvo.norep.gds")
  genofile

### load meta-data file
  samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/full_sample_metadata.28Sept2024_ExpEvo.csv")

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

### get the allele frequencies
  dat <- getData(snps=markers, samples=samps.new)

### your turn: write a figure legend for this figure.

### your turn:
### calculate the frequency of the inversion for each of your samples by averaging the frequency of all SNPs associated with the inversion
  dat.ag <- dat[,list(freq=mean(af_nEff, na.rm=T)),
                 list(inversion, exp_rep, sampleId, bio_rep)]

### your turn:
### plot the estimated allele frequency of each inversion across treatments

### some simple models to test for differences in frequency between the treatment groups
  mod <- lm(freq~exp_rep, dat.ag[inversion=="In(2L)t"])
  summary(mod)
  anova(mod)

### your turn: write a foreach loop that implements the linear model above for all inversions and returns the p-value associated with the treatment effect

  inversion.names <- c("In(2L)t", "In(3R)C", "In(3R)K", "In(3R)Mo", "In(3R)Payne")

  ## hint
    anova(mod)$Pr[1] ### this gets the pvalue. But, double check that this p-value is the one associated with your model

### your turn: add those p-values to your plot of the frequency.
