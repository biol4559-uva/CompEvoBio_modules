### libraries
  .libPaths(c("~/USERNAME", .libPaths())) ### <- you'll need to change this to your username
  library(ggplot2)
  library(data.table)
  library(foreach)  ### you will probably need to install this package
  library(doMC)     ### you will probably need to install this package
  registerDoMC(2)

## A first attempt.
  gen1 <- c(1,0, 1,0, 1,0, 1,0, 1,0) ### initialize our population
  gen2 <- sample(gen1, replace=T) ### The sample function draws samples from the input vector.

  table(gen1) ## this function tabulates the allele counts for gen1..
  table(gen2) ### and gen2.

  table(gen1)/length(gen1) ## calculate frequencies
  table(gen2)/length(gen2)

## Streamlining things a bit with `rbinom`
  rbinom(n=1, size=10, prob=.5) ### this returns the number of successes.
  rbinom(n=1, size=10, prob=.5)/10 ### this returns the frequency of the 'A' allele.

## One generation of drift using `rbinom`
  gen1 <- .5
  gen2 <- rbinom(n=1, size=10, prob=gen1)/10
  gen2

## Multiple generations
  gen1 <- .5
  gen2 <- rbinom(n=1, size=10, prob=gen1)/10
  gen3 <- rbinom(n=1, size=10, prob=gen2)/10
  gen4 <- rbinom(n=1, size=10, prob=gen3)/10

### First we define our parameters
    nGens <- 10
    popSize=10
    startingAlleleFreq=.5

### next, we generate a data.table that will hold our output.
  tmp <- data.table(gen=c(1:nGens), af=-1, popSize=popSize) ### we set `af` equal to -1 because it is just a placeholder. We include popSize for bookkeeping
  tmp[gen==1]$af <- startingAlleleFreq ### we need to initialize the very first generation at the specified allele frequency

### we run the for loop. We start at generation two because we need to use generation 1 as the starting generation
  for(i in 2:nGens) {
    tmp[gen==i]$af <- rbinom(1, popSize, tmp[gen==(i-1)]$af)/popSize
  }
  tmp

  ggplot(data=tmp, aes(x=gen, y=af)) + geom_line() + ylim(0,1)


## Multiple alleles
  test_fun <- function(popSize) {
    return(data.table(initialPopSize=popSize, newPopSize=popSize * 2))
  }

  tmp <- foreach(popSize.i=c(100, 500), .combine="rbind")%do%{
    test_fun(popSize=popSize.i)
  }
