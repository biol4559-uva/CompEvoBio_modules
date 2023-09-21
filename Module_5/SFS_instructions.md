# **Site Frequency Spectrum**

## Overview
In class so far, we have been thinking about coverages. The other important feature of your data is allele frequency. In this exercise, we will be identifying polymorphic sites, that is sites in the genome where there is more than one allele in the population. The [site-frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum) is basically a histogram that shows the number of sites in the genome that have allele frequency of 1%, 2%, 3% ... 99%. In most sexually reproducing species, the SFS has an expected shape. Most SNPs in the genome have low allele frequency, and very few are common. The image below depicts a standard SFS for a sample of ~30 individuals in a population. The SFS contains a ton of information about the population, believe it or not. Populations that have recently undergone an expansion ("Growth") tend to have a lot of very rare alleles, for instance. As an aside, in humans there are a ton of very-very rare alleles that are present because human population size has expanded at a [super-expontential rate](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3586590/) in the last several thousand years.


Goals:
1. Practice writing a function in R
2. Gain some insight into the consequences of not correcting for effective coverage

Objective:
1. How much does effective coverage change between the X and Autosomes if we have a pool of male flies.

## Instructions
1. Spin up an Rstudio job, copy the contents of [effective_coverage_template.R](/Module_5/effective_coverage_template.R)

2. Write a function that takes two parameters, the number of chromosomes in a sample and the read depth and returns effective coverage. Have it spit out a data.table with two three columns: rd, nChr, and effective_coverage.

3. Use your function to calculate effective coverage for an observed read depth of 40X and 20 flies. On an autosome, what will effective coverage be? On the X-chromosome, what will effective coverage be?

4. To complete this assignment, look at the Effective Coverage assignment on Canvas and help your "friend" solve their problem.
