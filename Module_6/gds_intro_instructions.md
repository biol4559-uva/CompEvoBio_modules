# **Intro to GDS files & metadata**

## Overview
Phew. Rivanna is at least operational, and we've been able to pull most of your samples off of the problematic drive. It is always a drag when tools fail in the middle of a project. But, we're on the right track and have a strategy for moving forward!

At this point, you have mapped your short reads and compiled your metadata. I've taken the output of your mapping jobs, and combined them all into what is called a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file. A VCF file is a text file that contains information about the read depth and allele frequencies for SNPs that are polymorphic in any particular sample. VCF files were originally built to handle individual sequencing data, where each entry represents a single fly (or person or cat or carp or whatever). The basic structure of a VCF file is as follows:

There are any number of lines that start with a "#". Those are header lines that define various features of the data. They all start with "##" at the begining of the line. Then there is the line that starts with a single "#". That defines the column names. There are ~10 columns that contain standardized information such as chromosome, position, reference allele, alternative allele, some information about the quality and effect of the paticular mutation, and so on. Then, starting around the 10th column, you have the sample information.

The sample information contains all of the important stuff, and the various bits of information are separated by a ":". The first field contains the genotype, encoded as 0/0, 0/1, 1/1. Because your data is pool-seq, this field is irrelevant. But, if you were dealing with individual based sequencing you would care about it the most. The second field contains the counts of the reference allele, the third the counts of the alternative allele, the fourth the total read depth, and the fifth the frequency of the alternate allele.

Our VCF file contains ~4.5 million SNPs. The reason that it is so large is that it contains your samples from the experimental evolution studies plus estimates of allele frequencies from fly populations worldwide. Fortunately, we do not need deal with this large file and instead we can use another version of it called a GDS file. The GDS file is a compressed version of the VCF file that we can work with using R. The nice thing about the GDS file is that we can access various pieces of information inside the VCF file with very little text processing.


Goals:
1. Learn how to interface with the GDS object using R

Objective:
1. Make a plot that shows the amount of "missing data" in your sample and test speed improvements with parallelization.
<p align="center">
  <img src="/Module_6/images/gds_plot.jpeg" width="1000"/>
</p>


## Instructions
1. Start an Rstudio job. Request 10 cores and 20Gb of RAM for 4 hours on the Standard partition.
2. Copy the [template](/Module_5/gds_intro_template.R) and work through the script. You'll have to draw on techniques that we have developed over the class to make the figure.
