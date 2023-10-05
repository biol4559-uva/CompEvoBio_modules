# **Principal component analysis**

## Overview
Our focal species, <i>D. melanogaster</i>, is a member of a clade of drosophilids that originated in central Africa. Intriguingly, these flies have accompanied humans in their migrations across the globe. Approximately 15,000 years ago, they ventured out of Africa into Eurasia, adapting to new environments and evolving alongside human populations. In more recent history, within the last 150 years, they spread to both Australia and the Americas, further diversifying. Now, armed with large genomic datasets, we have an opportunity to trace this journey. The <i>D. melanogaster</i> genome holds the clues about the shared and fixed mutations that have arisen during this [species' evolution](https://academic.oup.com/mbe/article/38/12/5782/6361628).

## Challenge
However, the <i>D. melanogaster</i> genome is complex, with over five million characterized mutations [(see DGRP paper)](https://www.nature.com/articles/nature10811). Analyzing this data is a formidable task. This is where Principal Component Analysis (PCA) comes into play. PCA is a versatile and powerful technique in data analysis and is used to reduce dimensionality. PCA enables us to distill complex and multidimensional genomic data into a more understandable and informative form, usually in a two-dimensional space. The goal of this class will be to learn how to apply PCA to transform high-dimensional data into a format that retains essential information while simplifying complexity.

Goals:
1. Use the Iris flower dataset in R to learn PCA.
2. See how your mapped <i>D. melanogaster</i> samples group according to their experimental treatments.
3. Find out which continental populations are most similar to your mapped <i>D. melanogaster</i> samples.

Objective:
1. Use PCA on the [Iris flower dataset](https://en.wikipedia.org/wiki/Iris_flower_data_set), and make a version of this figure from the example dataset. In this figure, the points represent an individual and the filled ellipses are group boundaries for the species.

<p align="center">
  <img src="/Module_7/images/Iris.comp.PCA.png" width="600"/>
</p>

## Instructions
1. Open a new Rstudio job. 

2. Copy the contents of [1.IrisPrincipalComponentAnalysis.R](/Module_7/1.IrisPrincipalComponentAnalysis.R) into a new file in your GitHub folder.

3. Run the code to generate the base figures of the Iris floral phenotypes and run the PCA.

4. Add to your version of [2.PCADrosophila.R](/Module_7/2.PCADrosophila.R) to use PCA on the genome-wide SNPs within your samples. You will have to merge the PCA results with the metadata and use ggplot to graph how your samples group according to experimental treatment. 

5. To complete this assignment, upload your script and figures to GitHub. In the Canvas assignment, upload your figure and the path to your github folder.
