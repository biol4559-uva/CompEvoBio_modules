### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)

### random allele frequencies at manyloci for two populations
  p1 <- runif(1000, 0, 1)
  p2 <- runif(1000, 0, 1)

### your turn: write a function to calculate Fst between two samples at many loci simultaneously. have your function output fst, average subpop heterozygosity, total heterozygosity, allele freq in p1, allele freq in p2
  FstFunction <- function(p1, p2) {

  }

### calculate
  fst.dt <- FstFunction(p1=p1, p2=p2)

### plot
  ggplot(data=fst.dt, aes(x=ht, y=fst)) + geom_point()
  ggplot(data=fst.dt, aes(x=p1, y=p2, color=fst)) + geom_point(size=.5)
  ggplot(data=fst.dt, aes(x=abs(p1-p2),y=fst,  color=fst)) + geom_point(size=.5)
