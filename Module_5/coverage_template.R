
### libraries
c  library(ggplot2)
  library(data.table)
  library(R.utils)
  library(patchwork)

### read the SYNC file
  sync <- fread("/scratch/COMPUTEID/coverage/ExpEvo_PRJEB5713_ancestral1_1_2007-MM-DD.tab.gz", sep2=":")

  setnames(sync, names(sync), c("chr", "pos", "ref", "A", "T", "C", "G", "N","del"))
  sync <- sync[,-"del",]
  sync[,id:=1:dim(sync)[1]]

### what is the total read depth?
  sync[,depth:=A+T+C+G+N]

### refresher on how to subset this data.table
  sync[chr=="2L"] ### this lets you see only sites on 2L
  sync[chr!="2L"] ### this lets you see all sites not on 2L

### aggregation is a technique to summarize data.
  sync.ag <- sync[,
                  list(mean_rd=mean(depth),
                       median_rd=median(depth)),
                  list(chr)]
  sync.ag

### Your turn:
### do different reference nucleotides have different coverages? (modify the aggregation call above)

### make a box and whisker plot of coverage across the chromosomes
  ### first, we are going to subsample the data
    setkey(sync, id)
    subsamp <- sort(sample(1:137567484, 100000, replace=F))
    sync.small <- sync[J(subsamp)]

  ### make the first plot
    coverage.boxplot <- ggplot(data=sync.small, aes(x=chr, y=depth)) + geom_boxplot()
    coverage.boxplot

### Your turn:
  ### How do you make that plot easier to interpret? Try removing the mitochondria.
  ### (I'll leave it to you to try and figure that out).
  ### Use the code example above and make a new object with your plot called "coverage_nomito.boxplot"

  ### combine the two plots together
    cov_boxplot <- coverage.boxplot + coverage_nomito.boxplot + plot_annotation(tag_levels="A", title="Coverage")
    cov_boxplot

### Poisson simulation. Let's work just on chromosome 2L
    setkey(sync, id)
    subsamp <- sort(sample(1:137567484, 1000000, replace=F))
    sync.small <- sync[J(subsamp)]

    sync.2L <- sync.small[chr=="2L"]
    nSNPs <- dim(sync.2L)[1]
    meanCov <- median(sync.2L$depth)

    sync.2L[,rand_pois:=rpois(nSNPs, meanCov)]

    coverage_density_plot <- ggplot(data=sync.2L) +
    geom_density(aes(x=depth), adjust=2) +
    geom_density(aes(x=rand_pois), color="red", adjust=2)

    coverage_density_plot

###########
### brief interlude to go explore UCSC Genome browser; see the instructions
###########

### flag repetitive regions
    ### first load in the repeat data. This data comes as a bed file, with chromosome, start, stop and the type
      rep <- fread("/scratch/aob2x/coverage/repeats.sort.merge.clean.bed")
      setnames(rep, names(rep), c("chr", "start", "stop", "type"))
      rep[,rep_region:=T]

    ### we need to do a little back-end work to align our sync object with the rep object. To to this alignment, we will use the `foverlaps` function
    ### this function overlaps ranges. Our sync file is in basepairs, so the range is one.
      sync[,start:=pos]
      sync[,stop:=pos]

    ### we need to specify which columns we are merging these two files on
      setkey(rep, chr, start, stop)
      setkey(sync, chr, start, stop)

      sync <- foverlaps(sync, rep[,-"type",], nomatch=NA)
      sync[is.na(rep_region), rep_region:=F]

### Your turn:
### test if sequence depth is different between repetitive and non-repetitive regions.
### use the Aggregate method described above to calculate the mean, median, min, max for each chromosome & rep region classification


### Your turn:
### Make a density plot from above and separates out the rep & non-rep regions.
### To make the final plot, you'll need to fiddle with two parameters in the geom_density call: fill & alpha
### Save that plot as `coverage_density_rep_plot`
### the code below will get you in the right direction.
      setkey(sync, id)
      subsamp <- sort(sample(1:137567484, 1000000, replace=F))
      sync.small <- sync[J(subsamp)]

      sync.2L <- sync.small[chr=="2L"]
      nSNPs <- dim(sync.2L)[1]
      meanCov <- median(sync.2L$depth)

      sync.2L[,rand_pois:=rpois(nSNPs, meanCov)]

### Make a composite plot
  layout <- "
  AC
  BC"

  coverage.boxplot +
    coverage_nomito.boxplot +
    coverage_density_rep_plot +
    plot_layout(design=layout) +
    plot_annotation(tag_levels="A", title="Coverage")
