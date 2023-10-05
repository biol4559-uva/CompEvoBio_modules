# Conduct principal component analyses on SNP data & estimate ancestry
# 10.5.2023

# Load packages
# Libraries
.libPaths(c("/scratch/csm6hg/biol4559-R-packages-newer")); .libPaths()

####################################################################################################
############## Move onto loading their own Drosophila GDS data and preform a PCA ###################
####################################################################################################

# Libraries
library(data.table)
library(ggplot2)
library(patchwork)
library(SeqArray)
library(adegenet)
library(ggforce)
library(zoo)
library(ggrepel)
library(doMC)
registerDoMC(4)

######### Load SNP data for larger metapopulation study of Drosophila melanogaster (DEST) #########

### open GDS file
genofile <- seqOpen("/standard/vol186/bergland-lab/biol4559-aob2x/data/dest.expevo.PoolSNP.001.50.25Sept2023.norep.ann.gds")

### get target populations
samps <- fread("/standard/vol186/bergland-lab/biol4559-aob2x/data/biol4559_sampleMetadata.csv")

# Extract samples from GDS
samps.gds <- seqGetData(genofile, "sample.id")

# Prepare SNP files - Select the samples you worked on!
# Enter your BioProject here:
bioproj = "SRP002024"
samps.new <- rbind(samps[sampleId %like% bioproj])

# Get subsample of data to work on. 
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps.new$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

# Choose biallelic sites - multiallelic positions are difficult to work with and violate assumptions
snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, sample.id=samps.new$sampleId, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

# Select sites 1) on autosomes, 2) common sites (intermediate allele frequencies), and 3) have low missingness.
seqSetFilter(genofile, 
             sample.id=samps.new$sampleId,
             snps.dt[chr%in%c("2L", "2R", "3L","3R")][missing < 0.05][af > 0.2]$variant.id)

# Extract updated filtered snps
snps.filt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

# Extract allele frequency data - this will take ~3 mins
ad <- seqGetData(genofile, "annotation/format/AD") # Alternative allele depth
dp <- seqGetData(genofile, "annotation/format/DP") # Full depth (coverage) at sites
dat <- ad$data/dp

# List number of: ( samples & sites )
dim(dat)

# Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"),
                       seqGetData(genofile, "position"), 
                       paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")

# Add samples to the rownames
rownames(dat) <- seqGetData(genofile, "sample.id")

# Reformat
dat_filt_t <- data.frame(t(dat))

# Add names
Sample_names <- names(dat_filt_t)

# Impute missing loci – as means
dat_filt_t_im = na.aggregate(dat_filt_t)

# Normalize the data using scale() function with centering and scaling
normalized_data <- scale(dat_filt_t_im, center = TRUE, scale = TRUE)

# Perform PCA using PCA function
pca.snps <- prcomp(t(normalized_data))

# variance proportion (%)
summ.pca <- summary(pca.snps)
prop.var <- data.table(var=c(summ.pca$importance[2,]*100),
                       PC=colnames(summ.pca$importance))

# plot the first and second PC
plot(pca.snps$x[,1], 
     pca.snps$x[,2],
     pch=20, cex=4)

# plot the second and third PC
plot(pca.snps$x[,2], 
     pca.snps$x[,3],
     pch=20, cex=4)

# Conditional fix sample names
if(bioproj=="SRP002024") {
  samps.new$sampleId <- gsub(samps.new$sampleId, pattern = "-", replacement = ".")
}

# Now it is your turn to merge your samples with the metadata object "samps" and:
# 1) plot the PCA using data.table and ggplot2 to show your experimental treatment groups!
# Make sure to utilize the geom_mark_eclipse() function for your treament groups, do the groupings make sense?

# Use: geom_label_repel(data=dt, aes(x=PC1, y=PC2, label=treatment)) from ggrepel to add labels!

##############################################################################
################## Move onto Larger species-wide dataset #####################
##############################################################################

# Restrict to your samples
samps.new <- rbind(samps[sampleId %like% bioproj])

# Restrict to well-known populations across D. melanogaster
# Remove D. simulans sample in dgn dataset
samps <- rbind(samps[set=="DrosRTEC"],
               samps[set=="DrosEU"],
               samps[set=="dgn"][!sampleId=="SIM_SIM_w501_1_NA-MM-DD"],
               samps.new)

# Reset and extract full dataset 
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

# Choose biallelic sites 
snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

# Select sites 1) on the 2L chromosome (much faster), 
# 2) common sites (intermediate allele frequencies), and 3) have low missingness.
seqSetFilter(genofile, snps.dt[chr%in%c("2L")][missing < 0.05][af > 0.2]$variant.id)

# Extract filtered snps
snps.filt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile))

# Extract allele frequency data - this will take ~3 mins
ad <- seqGetData(genofile, "annotation/format/AD") # Alternative allele depth
dp <- seqGetData(genofile, "annotation/format/DP") # Full depth (coverage) at sites
dat <- ad$data/dp

# List number of: ( samples & sites )
dim(dat)

# Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"),
                       seqGetData(genofile, "position"), 
                       paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")

# Add samples to the rownames
rownames(dat) <- seqGetData(genofile, "sample.id")

# Reformat for PCA
dat_filt_t <- data.frame(t(dat))

# Add names 
Sample_names <- names(dat_filt_t)

# Impute missing loci – as means
dat_filt_t_im = na.aggregate(dat_filt_t)

# Normalize the data using scale() function with centering and scaling
normalized_data <- scale(dat_filt_t_im, center = TRUE, scale = TRUE)

# Perform PCA using PCA function
pca.snps <- prcomp(t(normalized_data))

# variance proportion (%)
summ.pca <- summary(pca.snps)
prop.var <- data.table(var=c(summ.pca$importance[2,]*100),
                       PC=colnames(summ.pca$importance))

# plot the first and second PC
plot(pca.snps$x[,1], 
     pca.snps$x[,2],
     pch=20, cex=4)

# plot the second and third PC
plot(pca.snps$x[,2], 
     pca.snps$x[,3],
     pch=20, cex=4)

# Merge by Sample - fix weird naming issue
pca.snps1 <- data.table(sample=gsub(rownames(pca.snps$x), 
                                    pattern = "\\.", 
                                    replacement = "-"), 
                        PC1=pca.snps$x[,1],
                        PC2=pca.snps$x[,2],
                        PC3=pca.snps$x[,3],
                        PC4=pca.snps$x[,4],
                        PC5=pca.snps$x[,5])

# Now it is your turn to:
# 1) plot the PC data and get a sense of what variables are driving the samples to group together.

# To accomplish this, you will have to merge your samples with the metadata object "samps" and:
# 1) plot the PCA using data.table, ggplot2, and ggrepel to show your different experimental treatment groups within the species

# This section will find the closest sample in PC space:

# Function to find the closest sample by distance
dist.samp <- function(sampsToKeep, bioproj) {
  # sampsToKeep="ExpEvo_SRP002024_ACO_1_1975-MM-DD"; bioProject="SRP002024"
  print(sampsToKeep)
  
  dt.temp <- data.table(rbind(pca.snps1[!sample %like% bioproj],
                              pca.snps1[sample==sampsToKeep]))
  
  # Calculate PCA scores for all samples
  all_sample_scores <- dt.temp[,-c(1)]  # Assuming the first column is not part of PCA scores
  
  # Calculate pairwise Euclidean distances between all samples
  distances <- as.matrix(dist(all_sample_scores, diag = FALSE))
  diag(distances) <- Inf
  rownames(distances) <- dt.temp$sample
  colnames(distances) <- dt.temp$sample
  
  # Find the closest sample for each sample
  closest_indices <- apply(distances, 1, function(x) which.min(x[-1]))
  
  # Create a data.table to display the results
  result_df <- data.table(SampleIndex = 1:length(closest_indices),
                          ClosestSampleIndex = closest_indices,
                          ogSample=rownames(distances),
                          closestSample=names(closest_indices)[closest_indices])
  
  # Merge w/metadata
  result_df2 <- data.table(merge(result_df[ogSample==sampsToKeep], samps[,c(1,6,7,14)],
                                 by.x="closestSample", by.y="sampleId"))
  
  # progress message
  print(paste("Closest continental population:", result_df2$continent))
  
  # finis
  return(result_df2)
}

# Apply dist.samp to each input value using lapply
results_list <- lapply(samps.new$sampleId, function(sampsToKeep) {
  dist.samp(sampsToKeep, bioproj)
})

# Combine the results into a single data.table
combined_results <- do.call(rbind, results_list)

# 2) Show which group (continental population) is most like each of your samples using euclidean distances.
# Do your samples show concordance with what is the most similar continental population?
