# **Inversions**

## Overview

Inversions are a type of structural mutation where the order of genes along the chromosome are rearranged. Inversions can be quite large, and in Drosophila span 5Mb-15Mb (Mb=million basepairs). Recalling that each chromosomal arm of the Drosophila genome is ~25Mb, these inversions can rearrange large fractions of the genome. In Drosophila melanogaster, there are large inversions on all of the autosomes and they are referred to as In(2L)t, In(2R)Ns, In(3L)Payne, and so on. These particular inversions are said to be "cosmopolitan" because they are found in populations around the world. Individuals who have tow copies of a particular inversion are inversion homozygotes; individuals that have no copies of a particular inversion are called "standard" homozygotes; individuals can also be heterozygous for an inverted and standard allele.

<p align="center">
  <img src="/Module_8/images/inversion_order.jpeg" width="500"/>
</p>


Inversions also hold a special place in evolutionary biology, especially for Drosophila. This is because inversions were among the first identifiable mutations, and were studied before DNA was even shown to be the hereditary material or and before anyone had actually sequenced DNA let alone find the order of genes. Inversions can be identified by studying the banding pattern of polytene chromosomes extracted from Drosophila larval salivary glands. Polytene chromosomes are when bundles of hundreds or thousands of replicate copies of the genome. Basically, the cell is 1000N not 2N. Polytenization occurrs when a cell needs to produce a ton of an enzyme and does so by copying the chromsoome to increase the rate of transcription. Polytene chromosomes display distinct banding patterns that are very stereotypic. When individuals are heterozygous for inverted and standard alleles their polytene chromosomes form loops when the chromosomes pair. This looping pattern enabled early Drosophila scientists to track physical changes in the chromosomes before specific knowledge of DNA was available.

<p align="center">
  <img src="/Module_8/images/Inversioncartoon.jpeg" width="500"/>
</p>

Most species harbor inversion polymorphisms, and these inversions have been implicated in a variety of important functional, phenotypic variation. This is because inversions link together genetic variation by "supressing" recombination. During meiosis in an inverted/standard heterozygote, the chromosomes will pair and form loops (just like the polytene chromosomes). When cross-overs (recombination) happen between the inverted and standard chromosomes, the recombinant produces can have large deletions, large duplications, loose centromeres, or pick up a second cetromere. These events almost always lead to the selective death of the recombinant chromosomes and prevent them from being passed onto the next generation. Thus, inversions do not supress recombination in the physical sense, they supress recombination in the population sense.

The consequence of supressed recombination between inverted and standard chromosomes means that the genes between the inversion breakpoints co-evolve. Say, for instance, that for some special new adaptation to occur mutations must happen in two genes and that those genes are near the inversion breakpoints. Because the inversion supresses recombination, those mutations are always trasmitted together during meisosis and the special new adaptation will always manifest. Sometimes, inversions are referred to as 'supergenes' or 'co-adapted gene complexes' because of this feature.

Detecting inversion polymorphisms from short-read data can be a challenging endeavor. Fortunately, we know specific mutations that are tightly linked to the inversion. We can calculate the frequency of these mutations from the GDS file to get a sense of their frequency. You can find that file [here](/Module_8/inversion_markers_v6.txt).


Goals:
1. Learn a bit about inversions
2. Learn one way to test for differences in allele frequency
3. Start thinking about figure legends

Objective:
1. Make a plot showing the frequency of the inversions. Include the p-value from an anova for each inverison.
<p align="center">
  <img src="/Module_8/images/inversion_freq.jpeg" width="1000"/>
</p>

<p align="center">
  <img src="/Module_8/images/inversion_stats.jpeg" width="1000"/>
</p>
