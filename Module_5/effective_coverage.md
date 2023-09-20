# **Effective overage instructions**

## Overview
We have just looked at the statistics regarding the nominal coverage from a sample. The coverage of a sample is important because it reflects the precision of allele frequency estimtes. Precision is similar to variance (it is the inverse of variance!). Larger sample sizes give more precise esimates of allele frequency.

To understand the importance of read depth, consider the analogy of a coin. If we want to estimate the frequency that a coin flip gives a "tail", we would flip that coin multiple times. The mulitple flips of a coin is analagous to the multiple reads at a position in the genome. If you flip the coin ten times, maybe you get 3 heads and 7 tails. That does not mean that the true frequency of tails is 70%! The confidence intervals for our frequency estimate are calculated as: $\sqrt{pq}+(1+x)^2$

$$ SE = {sqrt{pq / n} $$

This sentence uses `$` delimiters to show math inline:  $\sqrt{3x-1}+(1+x)^2$


Goals:
1. Practice data-table and string manipulation.
2. Learn about the difference in effective coverage vs nominal coverage
2. Practice making plots and learn how to make a multi-panel plot

Objective:
1. Make this composite figure based on the test data
<p align="center">
  <img src="/Module_5/images/coverage_composite.jpeg" width="1000"/>
</p>
