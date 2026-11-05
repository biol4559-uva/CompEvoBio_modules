# **Fst outliers & sliding window**

## Overview
The goal of class today is to examine our full Fst output. Are there specific SNPs that have high Fst? What about genes? 

Goals:
1. Identify outlier windows, genes, SNPs
2. Determine if your outliers are potentially connected to your traits

Objective:
1. Combine your SNP-based FST analysis from last class. Identify outliers, and ask whether these outliers are enriched for certain types of features. 

## Instructions
1. Make sure that you have run your [Fst_gather.R](/Module_9/Fst_gather.R) script from last class, and that you have saved your `big_fstOutput.Rdata` object.
2. Using the script [Fst_summarization_version2.R](/Module_10/Fst_summarization_version2.R), examine various ways of collapsing your SNP based Fst estimates to sliding window based estimates. There is also a template for asking whether there is enrichment for different types of SNPs. This is a script that will help you identify outliers.
3. Using the example in [Fst_outlier_allele_freq_plot.R](/Module_10/Fst_outlier_allele_freq_plot.R), examine the patterns of change in allele frequency at your outliers. Are you convinced? At the bottom of that script, there is an example of how to study allele frequency change in natural popualation. This is a simple example, and you can explore other aspects of allele frequency change of your SNPs in natural popualations.
