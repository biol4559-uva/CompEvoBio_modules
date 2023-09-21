# **Coverage instructions**

## Overview
In our last class, you used the number of reads that map back to each chromosome to test if your samples were contaminated by the sister species, D. simulans and to test the sex ratio of your samples. Both of these calculations "coverage" as a metric. Coverage referrs to the number of reads that map to a region of the genome. You were calculating the average coverage across a chromosome. In this exercise, we will be calculating the actual coverage at every position in the genome.

Goals:
1. Practice data-table and string manipulation.
2. Learn about the difference in effective coverage vs nominal coverage
2. Practice making plots and learn how to make a multi-panel plot

Objective:
1. Make this composite figure based on the test data. Panel A shows a box-and-whisker plot of coverage across all chromosomes. Panel B shows the same, but removes the mitochdrion genome. Panel C shows the coverages for regions annotated as repetitive (blue) or not (red filled). The red line shows the expectated distribution based on a Poisson distribution.
<p align="center">
  <img src="/Module_5/images/coverage_composite.jpeg" width="1000"/>
</p>

## Instructions
1. Open a terminal window and run the following commands:
```
mkdir /scratch/COMPUTEID/coverage
cp -R /standard/vol186/bergland-lab/biol4559-aob2x/ /scratch/COMPUTEID/coverage/.
```

2. Open an Rstudio job, make sure to request 4 cores and 30Gb of RAM.

3. Make a new script, and copy the contents of [coverage_template.R](/Module_5/coverate_template.R) to a new script.

4. Work through lines 45 of the script. You'll have to come up with your own solution at one stage of the script to generate this figure:

<p align="center">
  <img src="/Module_5/images/coverage_boxwhisker.jpeg" width="1000"/>
</p>

5. In a frictionless universe where everything works according to simple models, the number of reads that cover any give site in the genome shoudl follow a [Poisson Distribution](https://en.wikipedia.org/wiki/Poisson_distribution). The Poisson distribution governs the number of events that occur over a given unit of time or space. One example of a Poisson distribution would be the number of students to sign into the AFC per hour across the course of a week. There are 7*24 observations (1 hour each), and the counts of students in each window across the week form our distribution. The Poisson distribution is determined by a single parameter, the average. <br>

In theory, another example of a Poisson distribution is the number of reads that map to the genome, per position in the genome. Read depth will only follow a Poisson distribution if each read is uniquely maps back to the genome, all reads have the same probability of mapping back to the genome, and all regions of the genome are equally likely to be sequenced. This is a tall order! <br>

We will assess whether our mapping data conform to a Poisson distribution by simulating a Poisson and comparing our real data to the observed data.

6. Examine the code in the `Poisson simulation` section of the script. How do you interpret the output?

7. How can we interpret the patterns that we observe? What aspects of the genome are contributing to this signal? One hypothesis to test is that [repetitive sequences](https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)) are causing us issues. Repetitive sequences are regions of the genome that have small, repeating sequences such as a microsattelite, a tandem repeat, etc. These are regions that are very hard to map short reads uniquely to. Microsat motifs occur throughout the genome and individuals will have variation in how many repeats they have. This can cause extrememly high or low read depths. Let's confirm manually that some of our very high or low read depth regions are likely to come from these repetitive sequences

8. Find a site with very high coverage and record its chromosome/position. Navigate to the [UCSC genome-browser](https://genome.ucsc.edu). Hover over "Genomes", and click on "Other"
<p align="center">
  <img src="/Module_5/images/UCSC_1.jpeg" width="1000"/>
</p>

9. Next, click on "Fruitfly", and make sure that "Aug 2014...dm6" is selected from the dropdown menu.
<p align="center">
  <img src="/Module_5/images/UCSC_2.jpeg" width="1000"/>
</p>

10. Enter your chosen site in the coordinates window using the format: "chr2L:23332290..23332290". On the next screen, make sure that all of the Repeat tracks are selected to "dense", and then hit "refresh". The Repeat tracks include "Interrupted Rpts", "Microsatellite", "Simple Repeats", "WM+SDust", "RepeatMasker". All of these tracks are are the output of different programs that predict repeat sequence from a genome sequence. Does your wierd site land in one of these repeat regions? Mine does, as indicated by the black line.
<p align="center">
  <img src="/Module_5/images/UCSC_3.jpeg" width="1000"/>
</p>

11. You can try zooming out. In my case, this particular region is gene-poor, has low conservation, and has multiple repeat annotations. We'd want to filter this site, plus all others like it from our dataset.
<p align="center">
  <img src="/Module_5/images/UCSC_5.jpeg" width="1000"/>
</p>

12. By comparison, let's look at our Adh F/S polymorphism. The location is "chr2L:14617051..14617051". If you zoom out a bit, you can see the structure of the gene and the high degree of conservation between species. The F/S polymorphism is a known polymorphism so it shows up in the EVA track.
<p align="center">
  <img src="/Module_5/images/UCSC_6.jpeg" width="1000"/>
</p>

13. The UCSC genome browswer is a great way to explore the Drosophila genome. The other excellent resource is [FlyBase](https://flybase.org), and specifically the [Jbrowse broswer](https://flybase.org/jbrowse/?data=data%2Fjson%2Fdmel&loc=2L%3A14617051..14617051&tracks=Gene_span%2CRNA&highlight=), which can display more information about functional characterization of the gene (explore the Gene Expression tracks)

14. Back to our data.... I've pulled together a file that defines the boundaries of all known repetitive regions in the reference genome. Work through the R-script starting at the "Flag Repetitive Regions"

15. To complete this assignment, upload your final composite figure. Feel free to fiddle with colors, layout, etc.
