# **Fst**

## Overview
The goal of class today is to implement a Fst analysis. First, we will do a small warm-up exercise to implement our own simple Fst test. Then, we will use a more sophisticated method that will enable us to identify Fst outliers. The more sophisticated method is called [outFLANK](http://www.jstor.org/stable/10.1086/682949).

Goals:
1. Learn about Fst
2. Learn about scatter-gather approaches

Objective:
1. Make your own Fst plot. You'll have to make some decisions about what you want to keep

## instructions
1. First, open the [Fst function script](/CompEvoBio_modules/Module_9/Fst_function_template.R). We'll talk in class about Fst and you will be tasked with writing a fuction that calculates it. This is just a warm-up exercise.

2. Calculating Fst between two populations using a naive estimator (like in step 1) is fast and easy. Doing a more proper analysis when there are many more samples is a bit more complicated. It is also slower. We are going to break the task into windows and then submit an array job where each job processes a window. Then, we'll collect the jobs together. That way, hopefully, things will run more quickly.

3. Open the [outFLANK_Fst.R script](/CompEvoBio_modules/Module_9/outFLANK_Fst.R). You'll need to walk thorugh this script (much of it will be familiar) and modify to point to your specific directory. Is there any additional information that you want? Do you want to calculate Fst between every sample? Or, do you want to target some samples? How do you want to deal with replicates?

4. Any operation that you perform on a single SNP is likely going to be time consuming. While we are running through your data, consider:
  • the Population branch statistic (this is a 3-population test)
  <p align="center">
    <img src="/Module_9/pbs.jpeg" width="1000"/>
  </p>

• generalized linear model:
  ```
  summary(glm(af_nEff~TREATMENT, family=binomial(), data=DATA, weights=DATA$nEff))
  ```
