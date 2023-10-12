# **Sliding windows**

## Overview

Sliding window plots are good ways to look at average signal across the genome. We will continue to develop these in the next few weeks but I want to give you a flavor of them now.

In a sliding window analysis, you take a window of the genome (say 10,000bp) and then perform some summarization of that region. The summarization could be, for instance, the number of polymorphisms. It could be average read depth. Or, it could be the a measure of how different allele frequencies are between treatment groups. Then, you slide over a certain number of basepairs and take another window.

Goals:
1. Implement a SW analysis. I'll give you the template, and then you need to come up with a creative application. Work with your group.

Objective:
1. Make at least three sliding window plots. Describe the y-axis, and provide a brief interpretation of your figure. Here is one example:

This figure shows the average change in allele frequency (y-axis) between 'contsys' and 'virsys' treatments. The y-axis represents the average absolute difference in allele frequency for all SNPs wihtin a 10,000 bp window. There is one standout region, on chromosome 3L and this region contains XYZ genes.
<p align="center">
  <img src="/Module_8/images/sliding_window.jpeg" width="1000"/>
</p>
