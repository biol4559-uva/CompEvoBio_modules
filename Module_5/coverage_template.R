# ijob -A biol4559-aob2x -c2 -p instructional --mem=12G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(ggplot2)
  library(data.table)
  .libPaths(c("/scratch/aob2x/biol4559/biol4559-R-packages/", .libPaths()))
  library(R.utils)

### read the SYNC file
  sync <- fread("/scratch/aob2x/biol4559/ExpEvo_PRJEB5713_ancestral1_1_2007-MM-DD.tab.gz", sep2=":")
  setnames(sync, names(sync), c("chr", "pos", "ref", "A", "T", "C", "G", "N","del"))
  sync <- sync[,-"del",]
  sync[,id:=1:dim(sync)[1]]


### what is the total read depth?
  sync[,depth:=A+T+C+G+N]

### aggregation is a technique to summarize data
  sync.ag <- sync[,
                  list(mean_rd=mean(depth),
                       median_rd=median(depth)),
                  list(chr)]
  sync.ag

### do different reference nucleotides have different coverages?

### make a box and whisker plot of coverage across the chromosomes
  ### first, we are going to subsample the data
  setkey(sync, id)

  subsamp <- sort(sample(1:137567484, 10000, replace=F))

  sync.small <- sync[J(subsamp)]
  ggplot(data=sync, aes(x=chr, y=depth)) + geom_boxplot()

###########
### brief interlude to go explore UCSC Genome browser; see the instructions
###########

### flag repetitive regions
