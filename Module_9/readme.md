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

2. Calculating Fst between two populations using a naive estimator (like in step 1) is fast and easy. Doing a more proper analysis when there are many more samples is a bit more complicated. It is also slower. We are going to break the task into windows and then submit an array job where each job processes a window. Then, we'll collect the jobs together. That way, hopefully, things will run more quickly. Use your scatter gather approach from last week, and try scattering across all windows (~2000 jobs). [Start here](/CompEvoBio_modules/Module_9/Fst_scatter.sh)

3. How do you know if your Fst values are larger than you'd expect by chance? That is a thorny problem. There are multiple ways of identifying Fst outliers. We will try one, called OutFLANK. The OutFLANK approach assumes that most of your genome is not under selection and then builds a null model based on that assumption. Then, it can take the full distribution of Fst values and identify if the outliers are lager than you'd expect based on that model. Open the [outFLANK_Window.R script](/CompEvoBio_modules/Module_9/outFLANK_Window.R). You'll need to walk thorugh this script (much of it will be familiar) and modify to point to your specific directory. Is there any additional information that you want? Do you want to calculate Fst between every sample? Or, do you want to target some samples? How do you want to deal with replicates?


4. You might also want to try generalized linear model or Fisher's exact tests or CMH tests on a Per-SNP basis:
  ```
  summary(glm(af_nEff~TREATMENT, family=binomial(), data=DATA, weights=DATA$nEff))

  for a FET, you'll need to make a 2x2 table with ref and alt counts for control and treatment

  for a CMH. you'll need to make a multidimensional table with ref and alt counts for control and treatment
  ```
