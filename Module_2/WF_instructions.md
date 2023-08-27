# **A Wrigh-Fisher simulation**

## Overview
In this exercise we will learn how to model drift. Drift is the random change in allele frequencies due to finite population size. Drift in small populations can be a strong evolutionary force, but in large populations is a weak force. When drift is a strong force (and thus when population sizes are small) allele frequencies will change quickly. Eventually drift will cause allele frequencies to go to 0% (lost) or 100% (fixed), and in small populations this happens faster. In large populations, this happens very slowly. This is why populations that have more individuals have more genetic diversity. By simulating the drift process using a Wrigh-Fisher simulation, we can begin to gain insight into the role of population size on genetic diversity

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

2. Reproduction: Although every allele is functionally equivalent, individuals can leave 0, 1, 2, or more offspring.
