# **A Wright-Fisher simulation**

## Overview
In this exercise we will learn how to model drift. Drift is the random change in allele frequencies due to finite population size. Drift in small populations can be a strong evolutionary force, but in large populations is a weak force. When drift is a strong force (and thus when population sizes are small) allele frequencies will change quickly. Eventually drift will cause allele frequencies to go to 0% (lost) or 100% (fixed), and in small populations this happens quickly. In large populations, this happens very slowly. This is why populations that have more individuals have more genetic diversity, and is how we can estimate historical population size in a species by studying levels of genetic diversity. By simulating the drift process using a Wright-Fisher simulation, we can begin to gain insight into the role of population size on genetic diversity.

Goals:
1. Practice more with R & Github & functions
2. Learn about one of the basic evolutionary forces, drift
3. Learn about mulit-core computing

Objective:
1. Make the classic 'drift' figure
<p align="center">
  <img src="/Module_2/images/wf.png" width="1000"/>
</p>


## Assumptions
1. Neutrality: One assumption of a basic drift model is that all alleles (that differ by state) at a locus are functionally equivalent. That is, they are neutral. <br>

2. Reproduction: Although every allele is functionally equivalent, individuals can leave 0, 1, 2, or more offspring. <br>

3. Population is in Hardy-Weinberg equilibrium. One implication of this assumption is that we only have to model allele frequencies, and can ignore genotypes.


## Start up an RStudio on Demand Job
Choose 2 cores. Load the libraries in the HWE_template.R script.

## A first attempt.
This structure of the WF simulation resembles the "bean-bag genetics" picture. We create the first generation ("gen1") of 2N gametes. In this example the population size is 5.

```
gen1 <- c(1,0, 1,0, 1,0, 1,0, 1,0) ### initialize our population. What do the 1's and 0's mean?
gen2 <- sample(gen1, replace=T) ### The sample function draws samples from the input vector.

table(gen1) ## this function tabulates the allele counts for gen1..
table(gen2) ### and gen2.

table(gen1)/length(gen1) ## calculate frequencies
table(gen2)/length(gen2)
```

If you run these lines of code multiple times, you will see that the frequency shifts in gen2.

## Streamlining things a bit with `rbinom`
The code above requires that we initialize our population with a certain number of individuals. If we wanted to model a population of size 1e6, this would be very tedious. We can simplify things by using the binomial random number generator function. <br>

Why do we use this function? In this model, we are assuming that polymorphisms are bi-allelic. This means that there are two alleles, different by state. This assumption fits for most SNP markers. It also fits for coins, where there are heads and tails (two alleles).

A binomial random number generator takes two basic parameters: the true frequency and the sample size. For the purposes of the WF simulation, the true frequency represents the frequency in the current generation; the sample size represents the population size of the next generation. Like before, the frequency in gen1 is 50% and there is a population size of 10.

```
###simulating 10 individual draws takes just about the same amount of time as 1e6 individuals
rbinom(n=1, size=10, prob=.5) ### this returns the number of successes.
rbinom(n=1, size=10, prob=.5)/10 ### this returns the frequency of the 'A' allele.

rbinom(n=1, size=1e6, prob=.5) ### this returns the number of successes.
rbinom(n=1, size=1e6, prob=.5)/1e6 ### this returns the frequency of the 'A' allele.

```

## One generation of drift using `rbinom`
All we really care about is allele frequencies...
```
gen1 <- .5
gen2 <- rbinom(n=1, size=10, prob=gen1)/10
gen2
```

## Multiple generations
How do we run this for multiple generations? One solution would be to do it manually:
```
gen1 <- .5
gen2 <- rbinom(n=1, size=10, prob=gen1)/10
gen3 <- rbinom(n=1, size=10, prob=gen2)/10
gen4 <- rbinom(n=1, size=10, prob=gen3)/10

c(gen1, gen2, gen3, gen4)
```

Would this solution scale well if we wanted to simulate many generations? Probably not and it is very tedious to write out. Instead, we can use a `for` loop and a pre-constructed matrix to store the output. Each iteration of the `for` loop is a generation.

```
library(ggplot2)
library(data.table)

### First we define our parameters
  nGens <- 10
  popSize=10
  startingAlleleFreq=.5

### next, we generate a data.table that will hold our output.
  tmp <- data.table(gen=c(1:nGens), af=-1, popSize=popSize) ### we set `af` equal to -1 because it is just a placeholder. We include popSize for bookkeeping
  tmp[gen==1]$af <- startingAlleleFreq ### we need to initialize the very first generation at the specified allele frequency
  tmp

### we run the for loop. We start at generation two because we need to use generation 1 as the stargin
  for(i in 2:nGens) {
    tmp[gen==i]$af <- rbinom(1, popSize, tmp[gen==(i-1)]$af)/popSize
  }
  tmp

  ggplot(data=tmp, aes(x=gen, y=af)) + geom_line() + ylim(0,1)

```
<p align="center">
  <img src="/Module_2/images/single.wf.jpg" width="500"/>
</p>


## make a WF function
Now, re-write your WF procedure to make it a function. Your function will need to take several parameters: nGens, popSize, startingAlleleFrequency. What sections of the code above should be included in your funciton?

## Multiple alleles
The specific trajectory of a single neutral allele only gives us a glimpse into the dynamics of drift. What we are really interested in is the distribution of allele frequencies at different loci through time. How do we generate this? One solution is to use a `foreach` loop. A `foreach` loop is like a `for` loop except that each iteration of the `foreach` loop is independent of the others. This is in contrast to a `for` loop where you can easily access information about previous iterations. The benefit of `foreach` loops are twofold: (1) they make bookkeeping and data-output much easier, and (2) they can be spead up on multicore computers, where each core runs the independent interaction of the foreach loop.

```
install.packages("foreach", "doMC")
library(foreach)
library(doMC)
registerDoMC(2)

test_fun <- function(nGens=10, popSize=500, startingAlleleFrequency=.5, locus) {
  return(data.table(nGens=nGens, locus=locus))
}

### runs each element of the for-each loop sequentially
tmp <- foreach(locus.i=c(1:10), .combine="rbind")%do%{
  test_fun(locus=locus.i)
}
tmp

### runs each element of the for-each loop in parallel
tmp <- foreach(locus.i=c(1:10), .combine="rbind")%dopar%{
  test_fun(locus=locus.i)
}
tmp
```

The `.combine` parameter tells `foreach` how to combine the different iterations.

Using all the tools above, your goal is now to build the objective figure below. This figure incorporates multiple alleles & multiple population sizes. If you are feeling especially curious, you might also modify starting allele frequency and plot that as separate rows of facets.


## Your objective:
Population size (the sample size) will determine the speed of drift. You can see this property by observing that allele frequencies change more, per generation when population size is small compared to when it is large. Using the information provided above, write a script that generates this basic figure. This figure includes multiple loci, and various population sizes. You are free to add your own artistic flair. 

To complete this project, user GitHub to upload your script and figure to a new folder in your repo and submit the link.

Note: the separate panels of the plot are called "facets" in ggplot jargon. To add facets to your plot use the `facet_grid()` function. Read the help file for that function to learn more about it.

You can also chagne the color scheme of your figure. Try installing `ggtheme` and explore some possibilities [here](https://ggplot2.tidyverse.org/reference/ggtheme.html)

<p align="center">
  <img src="/Module_2/images/wf.png" width="1000"/>
</p>
