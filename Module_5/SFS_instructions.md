# **Site Frequency Spectrum**

## Overview
In class so far, we have been thinking about coverages. The other important feature of your data is allele frequency. In this exercise, we will be identifying polymorphic sites, that is sites in the genome where there is more than one allele in the population. The [site-frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum) is basically a histogram that shows the number of sites in the genome that have allele frequency of 1%, 2%, 3% ... 99%. In most sexually reproducing species, the SFS has an expected shape. Most SNPs in the genome have low allele frequency, and very few are common. The image below depicts a standard SFS for a sample of ~30 individuals in a population. The SFS contains a ton of information about the population, believe it or not. Populations that have recently undergone an expansion ("Growth") tend to have a lot of very rare alleles, for instance. As an aside, in humans there are a ton of very-very rare alleles that are present because human population size has expanded at a [super-expontential rate](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3586590/) in the last several thousand years.

<p align="center">
  <img src="/Module_5/images/SFS.jpeg" width="1000"/>
</p>

Goals:
1. Explore the SFS of a Pool-Seq sample
2. Figure out how to identify polymorphic sites from your SYNC file.

Objective:
1. Determine if the SFS is different between genic regions (https://en.wikipedia.org/wiki/Coding_region) and non-genic regions. For the purposes of this assignment, we are considering genic regions to be coding sequence, UTR, and introns. Non-coding is everything else.
<p align="center">
  <img src="/Module_5/images/SFS_2.jpeg" width="1000"/>
</p>

## Instructions
1. Spin up an Rstudio job, and copy the contents of [SFS_template.R](/Module_5/SFS_template.R)
2. Work your way through the assignment and fill in the missing pieces.
