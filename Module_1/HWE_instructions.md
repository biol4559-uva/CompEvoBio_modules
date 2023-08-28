# **Instructions for HWE in-class assignment**

## Overview
In this exercise, you will examine the features of Hardy-Weinberg equilibrium to gain a better understanding of some aspects of allele frequency change, and to learn technical skills used throughout this class.

Objectives:
1. Use Rstudio on Rivanna for basic scripting.
2. Learn about basic object types in R, and make a function in R
3. Plot data that you generate using your function

Goals:
1.	Increased familiarity with R.
2.	Increased understanding of allele frequencies, precision, and statistical tests
3.  Understand the difference between long and wide data
---

## Steps
### 1.	Log into [OpenOneDemand](https://rivanna-portal.hpc.virginia.edu/pun/sys/dashboard) and start an interactive Rstudio job.
<p align="center">
  <img src="/Module_1/images/OOD1.png" width="750"/>
</p>

### 2. Select the options below. Why do these options matter? (we’ll discuss but here is the notes)

>a.	Version = which libraries to use<br>
>b.	Partition = what sort of computer your job gets given to, constraints, SU<br>
>c.	Hours = Your Rstudio session will terminate after this amount of time<br>
>d. Cores = billable SU (service unit)<br>
>e. Memory = for small jobs like this lesson we only need >Gb<br>
>f. Allocation = the account to bill SUs to. The instructional allocation prioritizes our jobs during our class-time so that we don’t have to wait. bio4559-aob2x<br>

<p align="center">
<img src="/Module_1/images/OOD2.png" width="750"/>
</p>

### 3. Create a new script in Rstudio
1. Copy the contents of the [HWE_template.R](/Module_1/HWE_template.R) to a new script in Rstudio. Try running the script. Play around for a minute modifying the functions. Can you write a new function that takes input parameters, performs calculations using two different sets of equations and returns the output from those two equations?

### 4. Write your own script/function to answer question 1.4 from Gillespie
1. The main goal is to recreate a version of this figure using a function that you write yourself in R and plot using ggplot.
<p align="center">
<img src="/Module_1/images/hwe_plot.png" width="750"/>
