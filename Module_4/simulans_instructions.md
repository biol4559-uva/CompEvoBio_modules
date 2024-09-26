# **D.simulans contamination**

## Overview
The Drosophila genus contains hundreds of species that live all around the world. Our focal species, <i>D. melanogaster</i>, is a member of a clade of drosophilids that originated in central Africa. For tens of thousands of years, it has lived amongst humans and human food waste. <i>D. melanogaster</i> and its close sister species <i>D. simulans</i> have colonized the world in the wake of human migration. They are thought to have migrated out of Africa into Eurasia about 15K years ago, and both are throught to have migrated to Australia and the Americas within the last 150 years. These two species live in similar habitats, and are morphologically nearly identical. Despite their similarity, they are also different in many ways including behaviors, physiology, environmental tolerances, plus (subtle) morphological differences. When sampling from the field, it can be challenging to identify individuals to these species perfectly. Even experienced fly researchers have been known to make mistakes.

The data we are using comes from "pool-seq" where many individual flies are combined, and DNA is extracted en-masse from the group of flies. If a <i>D. simulans</i> individual snuck in somehow, we can identify its presence by identifying short-read sequences that are more likely to come from <i>D. simulans</i> than <i>D. melanogaster</i>. To accomplish this, we have mapped our pools to a reference genome that contains the genome of both species. We can use the numnber of reads that map to each genome to estimate contamination.

Goals:
1. Practice using R.
2. Learn how to install packages from Bioconductor.
3. Learn how to interact with BAM files using R.

Objective:
1. Estimate contamination, and make a version of this figure for your data. In this figure, the points represent an relationship between the estimated contamination rate from the mapping approach and the true contamination rate (more on that below). The horizontal colored lines are the contamination estiamtes from samples in my study. In general, the contamination rates are estimated at ~1-2%. Given that there are 25 flies in my pool, this represents <1 fly and therefore these samples should not be considered suspect for contamination.

<p align="center">
  <img src="/Module_4/images/simulans_output.jpeg" width="1000"/>
</p>

## Instructions
1. Open a new Rstudio job. Select 5 cores. We will utilize a little bit of multicore computing in this example.

2. Copy the contents of [simulans_template.R](/Module_4/simulans_template.R) into a new file in your GitHub folder.

3. Run the code to generate the base figure with the black points. The data for this figure come from a "experiment" where we took short read sequence data from samples we were very certain were either D.melanogaster or D. simulans and mixed them together in various proportions. We generated pools of 40 diploid individuals with varying numbers and ran them through our mapping pipeline. The estimated proportions reflect the number of total reads that map to the D. simulans genome divided by the number of reads that map to both genomes. The correlation is pretty good!

4. Add to your version of [simulans_template.R](/Module_4/simulans_template.R) to estimate levels of contamination in your samples. You will have to take the parts of the script that load in the test data and fiddle with a few things to get it to work. One is figuring out how to get the list of bam files that you want. The other is some details inside the `%dopar%` call.

5. To complete this assignment, upload your script and figure to GitHub. In the Canvas assignment, upload your figure and the path to your github folder.
