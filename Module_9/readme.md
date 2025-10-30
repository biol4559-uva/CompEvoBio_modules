# **Fst**

## Overview
The goal of class today is to implement a Fst analysis. First, we will do a small warm-up exercise to implement our own simple Fst test. Then, we will use a more sophisticated method that will take into account the double binomial sampling of the Pool-Seq data (first round of binomial sampling = flies from the cage; second round is sequence reads from the flies).

Goals:
1. Learn about Fst
2. Learn about scatter-gather approaches

Objective:
1. Make your own Fst plot. You'll have to make some decisions about what you want to keep

## instructions
1. First, open the [Fst function script](/CompEvoBio_modules/Module_9/Fst_function_template.R). We'll talk in class about Fst and you will be tasked with writing a fuction that calculates it. This is just a warm-up exercise.

Remember the equation for Fst:  Fst=(Htot-Hwith)/Htot

2. Calculating Fst between two populations using a naive estimator (like in step 1) is fast and easy. Doing a more proper analysis when there are many more samples, and pool-seq is a bit more complicated. It is also slower. We are going to break the task into windows and then submit an array job where each job processes a window. Then, we'll collect the jobs together. That way, hopefully, things will run more quickly. Use your scatter gather approach from last week, and try scattering across all windows (~2000 jobs). [Start here](/CompEvoBio_modules/Module_9/Fst_scatter.sh), and then look [here](/CompEvoBio_modules/Module_9/Fst_template.R)

3. Now we need to collect our results and make sense of them. Use

4. You might also want to try generalized linear model or Fisher's exact tests or CMH tests on a Per-SNP basis:
  ```
  summary(glm(af_nEff~TREATMENT, family=binomial(), data=DATA, weights=DATA$nEff))

  for a FET, you'll need to make a 2x2 table with ref and alt counts for control and treatment

  for a CMH. you'll need to make a multidimensional table with ref and alt counts for control and treatment
  ```

5. You might also want to look at the "population branch statistic" which is based on three populations and Fst values between them. It calculates the branch lengths for population from the collective "center"

<p align="center">
  <img src="/Module_9/images/pbs.jpeg" width="500"/>
</p>
