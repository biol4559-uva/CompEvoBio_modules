### libraries
  .libPaths(c("/project/biol4559-aob2x/biol4559-R-packages/", .libPaths()))
  library(ggplot2)
  library(data.table)
  library(foreach)

### A first view of bean-bag geneticcs
### let's say that we have a population of 5 diploid individuals.
### That means that there are 10 alleles that differ by origin.
### if there is a locus with a polymorphism at 50% frequency, and everyone is heterozygote:

gen1 <- c(1,0, 1,0, 1,0, 1,0, 1,0)
gen2 <- sample(gen1, replace=T)

table(gen1)
table(gen2)

### 
