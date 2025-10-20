
####################################################################################################
############## Move onto loading their own Drosophila GDS data and preform a PCA ###################
####################################################################################################

# Libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(SeqArray)
  library(ggforce)
  library(zoo)
  library(ggrepel)
  library(doMC)
  registerDoMC(4)

######### Load SNP data for larger metapopulation study of Drosophila melanogaster (DEST) #########

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/compBio_SNP_29Sept2025/dest.all.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.gds")
  genofile

### load meta-data file
  samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/full_sample_metadata.90Sept2025_ExpEvo.csv")

# Extract samples from GDS
  samps.gds <- seqGetData(genofile, "sample.id")

# Prepare SNP files - Select the samples you worked on!
# Enter your BioProject here:
  bioproj = "PRJNA285429"
  samps.new <- samps[grepl(bioproj, sampleId)]

# Get basic SNP table
  seqSetFilter(genofile, sample.id=samps.new$sampleId)

  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile))

# Choose biallelic sites - multiallelic positions are difficult to work with and violate assumptions
  snp.dt <- snp.dt[nAlleles==2]
  seqSetFilter(genofile, sample.id=samps.new$sampleId, variant.id=snp.dt$variant.id)

# Select sites 1) on autosomes
  subsamp.variantIds <- sample(snp.dt[chr%in%c("2L", "2R", "3L","3R")]$variant.id, 1000)
  seqSetFilter(genofile,
               sample.id=samps.new$sampleId,
               variant.id=subsamp.variantIds)

# Extract allele frequency data - this will take ~3 mins
  ad <- seqGetData(genofile, "annotation/format/AD") # Alternative allele depth
  dp <- seqGetData(genofile, "annotation/format/DP") # Full depth (coverage) at sites
  dat <- ad$data/dp

# List number of: ( samples & sites )
  dim(dat)

# Add samples to the rownames
  rownames(dat) <- seqGetData(genofile, "sample.id")

# Transpose
  dat_t <- data.table(t(dat))

# We need to clean the data in two ways. First, we need to remove sites that have 100% missing data, then we need to impute the missing data
  # first, covert this file to long format
    dat_t[,variant.id:=seqGetData(genofile, "variant.id")]
    dat.long <- melt(dat_t, id.vars="variant.id", variable.name="sampleId", value.name="af")

  # next, calcualte missing rate per SNP and also calulate average allele frequency per SNP
    dat.long.missing <- dat.long[,list(missing=mean(is.na(af)), af_mean=mean(af, na.rm=T)), list(variant.id)]

  # only keep sites with low missing data
    setkey(dat.long, variant.id)
    setkey(dat.long.missing, variant.id)

    dat.long.filtered <- merge(dat.long, dat.long.missing[missing<1])
    dat.long.filtered

  # impute
    dat.long.filtered[is.na(af), af:=af_mean]

# now, to conduct PCA we need to turn our long data back into a matrix
    dat.wide <- dcast(dat.long.filtered, variant.id~sampleId, value.var="af")
    dim(dat.wide)

# Normalize the data using scale() function with centering and scaling
  dat.matrix <- dat.wide[,-"variant.id", with=F]
  normalized_data <- scale(dat.matrix, center = TRUE, scale = TRUE)

# Perform PCA using PCA function
  pca.snps <- prcomp(t(normalized_data))

# variance proportion (%)
  summ.pca <- summary(pca.snps)
  prop.var <- data.table(var=c(summ.pca$importance[2,]*100),
                         PC=colnames(summ.pca$importance))

# plot the first and second PC
  pca.dt <- data.table(PC1=pca.snps$x[,1], PC2=pca.snps$x[,2], sample=dimnames(pca.snps$x)[[1]])
  pca.dt[,trt:=tstrsplit(sample, "_")[[3]]]

  pca.dt
  ggplot(data=pca.dt, aes(x=PC1, y=PC2), size=10) + geom_point()


# Now it is your turn to merge your samples with the metadata object "samps" and:
# 1) plot the PCA using data.table and ggplot2 to show your experimental treatment groups.
# Make sure to utilize the geom_mark_eclipse() function for your treament groups, do the groupings make sense?

# Use: geom_label_repel(data=dt, aes(x=PC1, y=PC2, label=treatment)) from ggrepel to add labels!

##############################################################################
################## Move onto Larger species-wide dataset #####################
##############################################################################

# Restrict to your samples
  bioproj = "PRJNA285429"

# Restrict to well-known populations across D. melanogaster
# Remove D. simulans sample in dgn dataset
  samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/full_sample_metadata.90Sept2025_ExpEvo.csv", fill=T)

  samps.new <- samps[grepl(bioproj, sampleId)]
  samps <- rbind(samps[set=="DrosRTEC"],
                 samps[set=="DrosEU"],
                 samps[set=="dgn"][!sampleId=="SIM_SIM_w501_1_NA-MM-DD"],
                 samps.new)
  samps.n <- samps[,list(.N), list(sampleId)]
  table(samps.n$N)


### repeat the steps above but impute your data per locality (hint, you'll have to fix the locality tag for your samples to be equal to differentiate it from the other sites)
  seqResetFilter(genofile)
  subsamp.variantIds <- sample(snp.dt[chr%in%c("2L", "2R", "3L","3R")]$variant.id, 1000)
  seqSetFilter(genofile,
               sample.id=samps$sampleId,
               variant.id=subsamp.variantIds)

# Extract allele frequency data - this will take ~3 mins
  ad <- seqGetData(genofile, "annotation/format/AD") # Alternative allele depth
  dp <- seqGetData(genofile, "annotation/format/DP") # Full depth (coverage) at sites
  dat <- ad$data/dp

# List number of: ( samples & sites )
  dim(dat)

# Add samples to the rownames
  rownames(dat) <- seqGetData(genofile, "sample.id")

# Transpose
  dat_t <- data.table(t(dat))

# We need to clean the data in two ways. First, we need to remove sites that have 100% missing data, then we need to impute the missing data
  # first, covert this file to long format
    dat_t[,variant.id:=seqGetData(genofile, "variant.id")]
    dat.long <- melt(dat_t, id.vars="variant.id", variable.name="sampleId", value.name="af")

  # next, calcualte missing rate per SNP and also calulate average allele frequency per SNP
    setkey(dat.long, sampleId)
    dat.long <- merge(dat.long, samps, by="sampleId")
    dat.long.missing <- dat.long[,list(missing=mean(is.na(af)), af_mean=mean(af, na.rm=T)), list(variant.id, locality)]

  # only keep sites with low missing data
    setkey(dat.long, variant.id, locality)
    setkey(dat.long.missing, variant.id, locality)

    dat.long.filtered <- merge(dat.long, dat.long.missing[missing<1])
    dat.long.filtered

  # impute
    dat.long.filtered[is.na(af), af:=af_mean]

  # now, to conduct PCA we need to turn our long data back into a matrix
      dat.wide <- dcast(dat.long.filtered, variant.id~sampleId, value.var="af")
      dim(dat.wide)

  # Normalize the data using scale() function with centering and scaling
    dat.matrix <- dat.wide[,-"variant.id", with=F]
    dat.matrix <- na.omit(dat.matrix)
    normalized_data <- scale(dat.matrix, center = TRUE, scale = TRUE)

  # Perform PCA using PCA function
    pca.snps <- prcomp(t(normalized_data))

  # variance proportion (%)
    summ.pca <- summary(pca.snps)
    prop.var <- data.table(var=c(summ.pca$importance[2,]*100),
                           PC=colnames(summ.pca$importance))

  # plot the first and second PC
    pca.dt <- data.table(PC1=pca.snps$x[,1], PC2=pca.snps$x[,2], sample=dimnames(pca.snps$x)[[1]])
    pca.dt[,trt:=tstrsplit(sample, "_")[[3]]]
