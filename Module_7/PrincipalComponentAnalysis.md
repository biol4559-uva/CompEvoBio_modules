# **Principal component analysis**

## Overview
Our focal species, <i>D. melanogaster</i>, is a member of a clade of drosophilids that originated in central Africa. Intriguingly, these tiny flies have embarked on an epic journey, accompanying humans in their migrations across the globe. Approximately 15,000 years ago, they ventured out of Africa into Eurasia, adapting to new environments and evolving alongside human populations. In more recent history, within the last 150 years, they spread to new frontiers in Australia and the Americas, further diversifying their populations. Now, armed with the power of genomic data, we have an opportunity to trace this remarkable journey. The <i>D. melanogaster</i> genome is a treasure trove of information, holding clues about the shared and fixed mutations that have arisen during this species' evolution.

## Challenge
However, the <i>D. melanogaster</i> genome is complex, with over five million characterized mutations. Analyzing this data is a formidable task. This is where Principal Component Analysis (PCA) comes into play. PCA is a versatile and powerful technique in the realm of data analysis and dimensionality reduction. It enables us to distill complex and multidimensional genomic data into a more understandable and informative form. The goal of this class will be to learn how to apply PCA to transform high-dimensional data into a format that retains essential information while simplifying the complexity.

Goals:
1. Use the Iris flower dataset in R to learn PCA.
2. See how your mapped <i>D. melanogaster</i> samples group according to their experimental treatments.
3. Find out which continental populations are most similar to your mapped <i>D. melanogaster</i> samples.

Objective:
1. Use PCA on the Iris flower dataset, and make a version of this figure from the example dataset. In this figure, the points represent an individual and the filled ellipses are group boundaries for the species.

<p align="center">
  <img src="/Module_7/images/Iris.comp.PCA.pdf" width="1000"/>
</p>

## Instructions
1. Open a new Rstudio job. 

2. Copy the contents of [1.IrisPrincipalComponentAnalysis.R](/Module_7/1.IrisPrincipalComponentAnalysis.R) into a new file in your GitHub folder.

3. Run the code to generate the base figures of the floral phenotypes and run the PCA.

4. Add to your version of [2.PCADrosophila.R](/Module_7/2.PCADrosophila.R) to use PCA on the genome-wide SNPs within your samples. You will have to merge the PCA results with the metadata and use ggplot to ultimately graph how your samples group according to experimental treatment. 

5. To complete this assignment, upload your script and figures to GitHub. In the Canvas assignment, upload your figure and the path to your github folder.
